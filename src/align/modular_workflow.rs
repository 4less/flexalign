use std::cmp::min;

use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{
    align::{
        common::{Print},
    },
    database::common::FlexalignDatabase,
    flexalign::time,
    options::Options,
    GOLDSTD_EVAL,
};

use super::{
    common::{
        Align, AnchorExtractor, AnchorPair, Heuristic,
        KmerExtractor, Or, PAFOutput, PairedAnchorExtender, PairedAnchorExtractor,
        PairedAnchorMAPQ, RangeExtractor, SAMOutput, SeedExtractor, StdPairedAnchorMAPQ,
    },
    process::{
        alignment::ani_abort_score,
        evaluate::{self, get_id_from_header},
    },
    stats::Stats,
};
use itertools::Itertools;

#[derive(Clone)]
pub struct Modular<
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor<C>,
    RE: RangeExtractor<C, F>,
    SE: SeedExtractor<F>,
    AE: AnchorExtractor,
    PO: PAFOutput,
    SO: SAMOutput,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,
    pub db: &'a D,
    pub kmer_extractor: KE,
    pub range_extractor: RE,
    pub seed_extractor: SE,
    pub anchor_extractor: AE,

    pub rec_rev: OwnedFastqRecord,
    pub(crate) output: Or<PO, SO>,
}

impl<
        'a,
        const C: usize,
        const F: usize,
        KE: KmerExtractor<C>,
        RE: RangeExtractor<C, F>,
        SE: SeedExtractor<F>,
        AE: AnchorExtractor,
        PO: PAFOutput,
        SO: SAMOutput,
        D: FlexalignDatabase,
    > Modular<'a, C, F, KE, RE, SE, AE, PO, SO, D>
{
    //RE, SE,
    pub fn run(&mut self, rec: &RefFastqRecord, stats: &mut Stats) -> () {
        stats.reads_processed += 1;

        let (duration, kmers) = time(|| self.kmer_extractor.generate(rec, stats));
        stats.time_get_kmers += duration;

        let (duration, ranges) = time(|| self.range_extractor.generate(kmers, rec, stats));
        stats.time_get_ranges += duration;

        let (duration, seeds) = time(|| self.seed_extractor.generate(ranges, stats));
        stats.time_range_header += duration;
        stats.seeds += seeds.len();

        let (duration, anchors) = time(|| {
            self.anchor_extractor
                .generate(seeds, rec.seq().len(), stats)
        });
        stats.time_range_header += duration;
        stats.anchors += anchors.len();

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return;
        }

        let (_duration, _) = time(|| {
            anchors.sort_unstable_by_key(|a| {
                -((a.core_matches() - a.mismatches as usize - a.indels()) as i64)
            });
        });

        let (duration, _) = time(|| {
            rec.reverse_complement(&mut self.rec_rev);
        });
        stats.time_reverse_complement += duration;
        stats.time_anchor_sorting += duration;

        let best = anchors.first().unwrap();
        let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
        let reference = &self.db.get_reference(best.reference as usize).unwrap();

        let best_corelen = best.core_matches() - best.mismatches as usize - best.indels();
        let second_best_corelen = if anchors.len() > 1 {
            let second_best = anchors.get(1).unwrap();
            second_best.core_matches() - second_best.mismatches as usize - second_best.indels()
        } else {
            0
        };

        let pseudo_mapq = best_corelen - second_best_corelen;

        // Compile time switch
        if GOLDSTD_EVAL {
            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec.head().len())]
                == &rec.head()[..min(ref_string.len(), rec.head().len())];
            stats
                .gold_std_evaluation
                .as_mut()
                .unwrap()
                .add(correct, pseudo_mapq as u64);
        }

        if self.output.has_a() {
            self.output.a.as_mut().unwrap().write(
                &String::from_utf8_lossy(rec.head()),
                rec.seq().len(),
                best.seeds.first().unwrap().qbegin() as i32,
                best.seeds.last().unwrap().qend() as i32,
                best.forward,
                ref_string,
                reference.len(),
                best.seeds.first().unwrap().rbegin() as i32,
                best.seeds.last().unwrap().rend() as i32,
                best.seed_count,
                0,
                pseudo_mapq as u8,
            );
        }
    }
}

#[derive(Clone)]
pub struct ModularPE<
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor<C>,
    RE: RangeExtractor<C, F>,
    SE: SeedExtractor<F>,
    AE: PairedAnchorExtractor,
    AS: PairedAnchorExtender,
    PO: PAFOutput,
    SO: SAMOutput,
    A: Align + Heuristic + Send,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,
    pub db: &'a D,
    pub kmer_extractor_fwd: KE, // Extracts k-mers of relevance (e.g. minimizer)
    pub kmer_extractor_rev: KE,
    pub range_extractor_fwd: RE, // Extracts range where positions for a core-mer are stored
    pub range_extractor_rev: RE,
    pub seed_extractor_fwd: SE, // Extracts seeds from ranges
    pub seed_extractor_rev: SE,
    pub anchor_extractor: AE, // Generates anchors from seedlist
    pub anchor_extender: AS,  // Extends anchors with reference

    pub align: A,

    pub output: Or<PO, SO>,

    pub rec_fwd_revc: OwnedFastqRecord,
    pub rec_rev_revc: OwnedFastqRecord,
}

impl<
        'a,
        const C: usize,
        const F: usize,
        KE: KmerExtractor<C>,
        RE: RangeExtractor<C, F>,
        SE: SeedExtractor<F>,
        AE: PairedAnchorExtractor,
        AS: PairedAnchorExtender,
        PO: PAFOutput,
        SO: SAMOutput,
        A: Align + Heuristic + Send,
        D: FlexalignDatabase,
    > ModularPE<'a, C, F, KE, RE, SE, AE, AS, PO, SO, A, D>
{
    //RE, SE,
    pub fn run(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        stats: &mut Stats,
    ) -> () {
        stats.reads_processed += 2;

        log::trace!("\n####################################################################\n#### New pair \n####################################################################");

        if GOLDSTD_EVAL {
            let header_str = String::from_utf8_lossy(rec_fwd.head());
            let id = get_id_from_header(&header_str, self.db);
            log::trace!("True ID: {}\n", id);
        }

        log::trace!("{}\n{}", rec_fwd.to_string(), rec_rev.to_string());

        // Extract minimizer
        let (duration, kmers_fwd) = time(|| self.kmer_extractor_fwd.generate(rec_fwd, stats));
        stats.time_get_kmers += duration;
        let (duration, kmers_rev) = time(|| self.kmer_extractor_rev.generate(rec_rev, stats));
        stats.time_get_kmers += duration;

        // Get ranges from minimizers
        let (duration, ranges_fwd) =
            time(|| self.range_extractor_fwd.generate(kmers_fwd, rec_fwd, stats));

        stats.time_get_ranges += duration;
        let (duration, ranges_rev) =
            time(|| self.range_extractor_rev.generate(kmers_rev, rec_rev, stats));
        stats.time_get_ranges += duration;

        log::trace!(
            "\n##### Ranges FWD \n{}\n",
            ranges_fwd
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    format!(
                        "{} - {} {} \nVRANGE:{}\n {}",
                        index,
                        range.qpos,
                        range.fmer.to_string().unwrap(),
                        "", //range.vrange.to_string(),
                        range.range_size
                    )
                })
                .collect_vec()
                .join("\n")
        );

        log::trace!(
            "\n##### Ranges REV \n{}\n",
            ranges_rev
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    format!(
                        "{} - {} {} \nVRANGE:{}\n {}",
                        index,
                        range.qpos,
                        range.fmer.to_string().unwrap(),
                        "", //range.vrange.to_string(),
                        range.range_size
                    )
                })
                .collect_vec()
                .join("\n")
        );

        if false && GOLDSTD_EVAL {
            let mut _index = 0;
            let mut correct = false;
            let mut _correct_id = 0;
            for (idx, range) in ranges_rev.iter().enumerate() {
                let kmer_len = if range.vrange.header.is_some() {
                    31
                } else {
                    15
                };
                range
                    .vrange
                    .best_flex_match(&range.fmer, |rpos, rval, distopt| {
                        let refstr = &self.db.get_rname(rval as usize).unwrap();
                        correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())]
                            == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                        if correct && idx > 0 {
                            eprintln!(
                                "{}/{} index of first correct with {:?} at refpos {}",
                                idx,
                                ranges_fwd.len(),
                                distopt,
                                rpos
                            );
                        }
                    });

                if !correct {
                    let quals = rec_fwd.qual();
                    let _test = &quals[range.qpos..range.qpos + kmer_len];
                    let kmer_quals = &quals[range.qpos..range.qpos + kmer_len]
                        .iter()
                        .map(|q| q - 33);
                    let minq = kmer_quals.clone().into_iter().min();
                    eprintln!(
                        "--fwd-> {}: {:?} {:?}",
                        idx,
                        minq,
                        kmer_quals.clone().into_iter().join(",")
                    );
                }

                // vrange.all_matches(|rpos, rval| {
                //     let refstr = &self.db.get_rname(rval as usize).unwrap();
                //     correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                //     if correct && idx > 5 {
                //         correct_id = rval;
                //         eprintln!("{}/{} index of first correct with --- at refpos {} vrange size {}", idx, ranges_fwd.len(), rpos, vrange.len());
                //     }
                // });

                if correct {
                    _index = idx;
                    break;
                }
            }
            // if index > 0 {
            //     eprintln!("-----------------------------------");
            //     eprintln!("Index {}", index);
            //     eprintln!("{}\n", rec_fwd.to_string());
            //     eprintln!("-----------------------------------");
            //     eprintln!("CorrectID: {}", correct_id);
            //     for (ridx, (qpos, flex, range, len)) in ranges_fwd.iter().enumerate() {
            //         eprintln!("{}----\n{}", ridx, range.to_verbose_string());
            //     }
            //     exit(9);
            // }

            let mut _index = 0;
            let mut correct = false;
            let mut _correct_id = 0;
            //(qpos, flex, vrange, len)
            for (idx, range) in ranges_rev.iter().enumerate() {
                let kmer_len = if range.vrange.header.is_some() {
                    31
                } else {
                    15
                };
                range
                    .vrange
                    .best_flex_match(&range.fmer, |rpos, rval, distopt| {
                        let refstr = &self.db.get_rname(rval as usize).unwrap();
                        correct |= &refstr.as_bytes()[..min(refstr.len(), rec_rev.head().len())]
                            == &rec_rev.head()[..min(refstr.len(), rec_rev.head().len())];
                        if correct && idx > 0 {
                            eprintln!(
                                "rev {}/{} index of first correct with {:?} at refpos {}",
                                idx,
                                ranges_rev.len(),
                                distopt,
                                rpos
                            );
                        }
                    });

                if !correct {
                    let quals = rec_rev.qual();
                    let kmer_quals = &quals[range.qpos..range.qpos + kmer_len]
                        .iter()
                        .map(|q| q - 33);
                    let minq = kmer_quals.clone().into_iter().min();
                    eprintln!(
                        "--rev-> {}: {:?} {:?}",
                        idx,
                        minq,
                        kmer_quals.clone().into_iter().join(",")
                    );
                }
                // vrange.all_matches(|rpos, rval| {
                //     let refstr = &self.db.get_rname(rval as usize).unwrap();
                //     correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                //     if correct && idx > 5 {
                //         correct_id = rval;
                //         eprintln!("{}/{} index of first correct with --- at refpos {} vrange size {}", idx, ranges_fwd.len(), rpos, vrange.len());
                //     }
                // });

                if correct {
                    _index = idx;
                    break;
                }
            }
            // if index > 0 {
            //     eprintln!("-----------------------------------");
            //     eprintln!("Index {}", index);
            //     eprintln!("{}\n", rec_rev.to_string());
            //     eprintln!("-----------------------------------");
            //     eprintln!("CorrectID: {}", correct_id);
            //     for (ridx, (qpos, flex, range, len)) in ranges_rev.iter().enumerate() {
            //         eprintln!("{}----\n{}", ridx, range.to_verbose_string());
            //     }
            //     exit(9);
            // }
        }

        // ranges_fwd.iter().for_each(|x| {
        //     let qend = x.qpos + 31;
        //     if qend as usize > rec_fwd.seq().len() {
        //         log::debug!("->> {:?} {}/{}", x.qpos, qend, rec_fwd.seq().len());
        //         log::debug!("{}", x.vrange);
        //     }
        //     assert!(qend as usize <= rec_fwd.seq().len());
        // });

        // ranges_rev.iter().for_each(|x| {
        //     let qend = x.qpos + 31;
        //     if qend as usize > rec_rev.seq().len() {
        //         log::debug!("->> {:?} {}/{}", x.qpos, qend, rec_rev.seq().len());
        //         log::debug!("{}", x.vrange);
        //     }
        //     assert!(qend as usize <= rec_rev.seq().len());
        // });

        // Get Seeds from ranges
        let (duration, seeds_fwd) = time(|| self.seed_extractor_fwd.generate(ranges_fwd, stats));
        stats.time_range_header += duration;
        stats.seeds += seeds_fwd.len();
        let (duration, seeds_rev) = time(|| self.seed_extractor_rev.generate(ranges_rev, stats));
        stats.time_range_header += duration;
        stats.seeds += seeds_rev.len();

        // seeds_fwd.iter().for_each(|s| {
        //     let qend = s.qpos + s.length as u32;
        //     if qend as usize > rec_fwd.seq().len() {
        //         println!("->> {:?} {}/{}", s, qend, rec_fwd.seq().len());
        //     }
        //     assert!(qend as usize <= rec_fwd.seq().len());
        // });
        // seeds_rev.iter().for_each(|s| {
        //     let qend = s.qpos + s.length as u32;
        //     if qend as usize > rec_rev.seq().len() {
        //         println!("->> {:?} {}/{}", s, qend, rec_rev.seq().len());
        //     }
        //     assert!(qend as usize <= rec_rev.seq().len());
        // });

        log::trace!(
            "\n##### Seeds FWD \n{}\n\n",
            seeds_fwd
                .iter()
                .enumerate()
                .map(|(index, seed)| {
                    format!(
                        "{} - ref: {} qpos: {} rpos: {}, RNAME: {:?}",
                        index,
                        seed.rval,
                        seed.qpos,
                        seed.rpos,
                        self.db.get_rname(seed.rval as usize)
                    )
                })
                .collect_vec()
                .join("\n")
        );

        log::trace!(
            "\n##### Seeds REV \n{}\n\n",
            seeds_rev
                .iter()
                .enumerate()
                .map(|(index, seed)| {
                    format!(
                        "{} - ref: {} qpos: {} rpos: {}, RNAME: {:?}",
                        index,
                        seed.rval,
                        seed.qpos,
                        seed.rpos,
                        self.db.get_rname(seed.rval as usize)
                    )
                })
                .collect_vec()
                .join("\n")
        );

        // eprintln!("Header {} ... \nID {}", String::from_utf8_lossy(rec_fwd.head()), get_id_from_header(&String::from_utf8_lossy(rec_fwd.head()), self.db));
        // eprintln!("FWD: {}\nREV: {}", rec_fwd.to_string(), rec_rev.to_string());
        // eprintln!("Ref 1: {}", String::from_utf8_lossy(self.db.get_reference(1).unwrap()));

        // eprintln!("Seed number: {}", seeds_rev.len());
        // for seed in seeds_rev {
        //     eprintln!("----- Seed ({})\n{}", self.db.get_rname(seed.rval as usize).unwrap(), seed);
        //     let reference = self.db.get_reference(seed.rval as usize).unwrap();
        //     eprintln!("Q: {}", String::from_utf8_lossy(&rec_fwd.seq()[seed.query_range()]));
        //     eprintln!("R: {}", String::from_utf8_lossy(&reference[seed.reference_range()]));
        // }

        let (duration, anchors) = time(|| {
            self.anchor_extractor.generate(
                seeds_fwd,
                seeds_rev,
                rec_fwd.seq().len(),
                rec_rev.seq().len(),
                stats,
            )
        });
        stats.time_get_anchors += duration;
        stats.anchors += anchors.len();

        log::trace!(
            "\n#################################################\n##### Anchors\n#################################################\n{}",
            anchors
                .iter()
                .enumerate()
                .map(|(index, AnchorPair(r1, r2))| {
                    format!(
                        "\n####### Anchor Index: {} \n------Read /1: \n{}\n------Read /2: \n{}",
                        index,
                        if let Some(r1) = r1 {
                            r1.to_string()
                        } else {
                            "No anchor Read /1".to_string()
                        },
                        if let Some(r2) = r2 {
                            r2.to_string()
                        } else {
                            "No anchor Read /2".to_string()
                        },
                    )
                })
                .collect_vec()
                .join("\n")
        );

        // for AnchorPair(a1, a2) in anchors.iter() {
        //     match a1 {
        //         Some(a) => {
        //             let first: &AnchorSeed = &a.seeds[0];
        //             assert!(first.qend() <= rec_fwd.seq().len());
        //         }
        //         None => (),
        //     }
        //     match a2 {
        //         Some(a) => {
        //             // println!("Querylen: {},\nAnchor... {}", rec_rev.seq().len(), a);
        //             let first: &AnchorSeed = &a.seeds[0];
        //             assert!(first.qend() <= rec_rev.seq().len());
        //         }
        //         None => (),
        //     }
        // }

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return;
        }

        // eprintln!("Read: {}", String::from_utf8_lossy(rec_fwd.head()));
        let best_before = anchors.first().as_mut().unwrap().clone();

        // ######################################################################################################
        // Now here starts the reference-based portion of the algorithm. Before, no sequence comparison
        // Between query and reference is done
        // ######################################################################################################

        let (duration, _) = time(|| {
            rec_fwd.reverse_complement(&mut self.rec_fwd_revc);
            rec_rev.reverse_complement(&mut self.rec_rev_revc);
        });
        stats.time_reverse_complement += duration;

        let anchors_len = anchors.len();
        let max_hamming = 10;

        // assert!(extension_anchors.iter().all(|pair| pair.validate().is_ok()));

        // Assumes valid anchor seeds!!
        let (duration, extended_count) = time(|| {
            let extended_count = self.anchor_extender.extend(
                anchors,
                rec_fwd,
                &self.rec_fwd_revc,
                rec_rev,
                &self.rec_rev_revc,
                stats,
            );

            glidesort::sort_by_key(&mut anchors[0..extended_count], |AnchorPair(a1, a2)| {
                -(a1.as_ref().map_or(0, |a| a.score) + a2.as_ref().map_or(0, |a| a.score))
            });

            // assert!(extension_anchors.iter().all(|pair| pair.validate().is_ok()));
            extended_count
        });

        // Assumes sorted anchors !!
        let mut extension_anchors = &mut anchors[0..extended_count];

        stats.time_extend_anchors += duration;

        if GOLDSTD_EVAL {
            let mut any_correct = false;
            for AnchorPair(a1, a2) in extension_anchors.iter() {
                let rval = a1.as_ref().map_or_else(|| a2.as_ref().unwrap().reference, |a| a.reference);
                let ref_string = &self.db.get_rname(rval as usize).unwrap();
                let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())]
                    == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];
                any_correct |= correct;
            }
            log::trace!("Any correct? {}", any_correct);
        }

        log::trace!(
            "\n##### Sorted and Extended Anchors{}\n\n",
            extension_anchors
                .iter()
                .enumerate()
                .map(|(index, AnchorPair(r1, r2))| {
                    format!(
                        "\n####### Anchor Index: {} \n------Read /1: \n{}\n------Read /2: \n{}",
                        index,
                        if let Some(r1) = r1 {
                            r1.to_string()
                        } else {
                            "No anchor Read /1".to_string()
                        },
                        if let Some(r2) = r2 {
                            r2.to_string()
                        } else {
                            "No anchor Read /2".to_string()
                        },
                    )
                })
                .collect_vec()
                .join("\n")
        );

        //#####################################################################################
        //# At this point enforce valid seeds

        // extension_anchors
        //     .iter_mut()
        //     .enumerate()
        //     .for_each(|(i, (AnchorPair(a1, a2)))| {
        //         match a1 {
        //             Some(a) => {
        //                 if a.seed_status.is_valid()
        //                     && a.seeds.len() > 1
        //                     && a.seeds[0].qbegin() > a.seeds[1].qbegin()
        //                 {
        //                     eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
        //                     panic!("Z 1  {}", a);
        //                 }
        //             }
        //             _ => {}
        //         }

        //         match a2 {
        //             Some(a) => {
        //                 if a.seed_status.is_valid()
        //                     && a.seeds.len() > 1
        //                     && a.seeds[0].qbegin() > a.seeds[1].qbegin()
        //                 {
        //                     eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
        //                     panic!("Z 2  {}", a);
        //                 }
        //             }
        //             _ => {}
        //         }
        //     });

        let pseudo_mapq = StdPairedAnchorMAPQ::anchor_mapq(extension_anchors);
        let old_score = &extension_anchors[0].0.as_ref().map_or(0, |a| a.score)
            + &extension_anchors[0].1.as_ref().map_or(0, |a| a.score);


        // Assumes sorted anchors !!
        let anchors_len: usize = extension_anchors.len();
        let alignment_anchors =
            &mut extension_anchors[0..min(self.options.args.align_top_y, anchors_len)];

        let (duration, _) = time(|| {
            let mut min_score_1 = None;
            let mut min_score_2 = None;

            alignment_anchors
                .iter_mut()
                .enumerate()
                .for_each(|(_i, (AnchorPair(a1, a2)))| {
                    let reference = match a1 {
                        Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
                        None => &self
                            .db
                            .get_reference(a2.as_ref().unwrap().reference as usize)
                            .unwrap(),
                    };

                    match a1 {
                        Some(a) => {
                            let query = if a.forward {
                                rec_fwd.seq()
                            } else {
                                self.rec_fwd_revc.seq()
                            };
                            if query.len() == 0 {
                                a.score = 0i32;
                            } else {
                                if min_score_1.is_none() {
                                    min_score_1 =
                                        Some(ani_abort_score(0.5, 4, query.len() as i32).abs());
                                }
                                self.align.set_max_alignment_score(min_score_1.unwrap());
                                // eprintln!("Align max score: {}", min_score_1.unwrap());

                                if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                    log::debug!("1  {}", a);
                                }

                                // let status = a.smart_align(&mut self.align, query, reference, 10, min_score_1.unwrap());
                                // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_1.unwrap());

                                let status = a.align(
                                    &mut self.align,
                                    query,
                                    reference,
                                    10,
                                    min_score_1.unwrap(),
                                );

                                // let (qr, rr) = a.whole(query.len(), reference.len());
                                // let (duration, (score, cigar, status)) = time(|| self.align.align(&query[qr], &reference[rr]));

                                match status {
                                    super::common::Status::OK => stats.alignments_successful += 1,
                                    super::common::Status::Dropped => stats.alignments_dropped += 1,
                                    super::common::Status::Partial => stats.alignments_partial += 1,
                                }

                                let score = a.score;
                                // stats.time_offset += duration;
                                // stats.alignments += 1;
                                // a.score = score / -4;

                                let _ani = (1.0 - a.score as f64 / a.cigar().0.len() as f64);
                                // let ani: f64 = (1.0 - a.score as f64/cigar.0.len() as f64);
                                // let ani: f64 = (1.0 - a.score as f64/a.cigar().0.len() as f64);
                                // if score < -50 && score != std::i32::MIN {
                                //     eprintln!("{}/1: {} ANI: {}", i, score, ani);
                                // }

                                if score != std::i32::MIN && -score < min_score_1.unwrap() {
                                    // eprintln!("Set {} -> {}", min_score_1.unwrap(), -score);
                                    min_score_1 = Some(-score);
                                }
                                // eprintln!("{} (asize: {}) Set score {} {} {} {}", i, a.seeds.len(), score, a.score, (1.0 - a.score as f64/cigar.0.len() as f64),  String::from_utf8_lossy(&cigar.0));
                            }
                            // eprintln!("{}", query.len());
                        }
                        None => (),
                    };
                    match a2 {
                        Some(a) => {
                            let query = if a.forward {
                                rec_rev.seq()
                            } else {
                                self.rec_rev_revc.seq()
                            };
                            if query.len() == 0 {
                                a.score = 0i32;
                            } else {
                                if min_score_2.is_none() {
                                    min_score_2 =
                                        Some(ani_abort_score(0.5, 4, query.len() as i32).abs());
                                }

                                if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                    log::debug!("2  {}", a);
                                }

                                // self.align.set_max_alignment_score(min_score_2.unwrap());
                                // let status = a.smart_align(&mut self.align, query, reference, 10, min_score_2.unwrap());
                                // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_2.unwrap());

                                let status = a.align(
                                    &mut self.align,
                                    query,
                                    reference,
                                    10,
                                    min_score_2.unwrap(),
                                );

                                // let (qr, rr) = a.whole(query.len(), reference.len());
                                // let (duration, (score, cigar, status)) = time(|| self.align.align(&query[qr], &reference[rr]));

                                match status {
                                    super::common::Status::OK => stats.alignments_successful += 1,
                                    super::common::Status::Dropped => stats.alignments_dropped += 1,
                                    super::common::Status::Partial => stats.alignments_partial += 1,
                                }

                                // match status {
                                //     super::common::Status::OK => {
                                //         if a.reference_cigar_range.len() == 0 {
                                //             eprintln!("Invalid range... {:?}", a.reference_cigar_range);
                                //         }
                                //         if is_alignment_valid(&query, &reference[a.reference_cigar_range.clone()], &a.cigar().0) {
                                //             // print_alignment(&query, &reference[a.reference_cigar_range.clone()], &a.cigar().0);
                                //         } else {
                                //             // eprintln!("------------------------");
                                //             // eprintln!("Valid ? {}", a.validate_seeds(query, reference));
                                //             // eprintln!("{}", a);
                                //             // eprintln!("{}\n{}\n{}",
                                //             //     String::from_utf8_lossy(query),
                                //             //     String::from_utf8_lossy(&reference[a.reference_cigar_range.clone()]),
                                //             //     String::from_utf8_lossy(&a.cigar().0));
                                //             // panic!("Issue {:?}", a.reference_cigar_range.clone());
                                //             eprintln!("Flag issue {:?}", a.reference_cigar_range.clone());
                                //         }
                                //     },
                                //     _ => ()
                                // }

                                let score = a.score;
                                // stats.time_offset += duration;
                                // stats.alignments += 1;
                                // a.score = score / -4;

                                let _ani = (1.0 - a.score as f64 / a.cigar().0.len() as f64);
                                // let ani = (1.0 - a.score as f64/cigar.0.len() as f64);
                                // if score < -50 && score != std::i32::MIN {
                                //     eprintln!("{}/2: {} ANI: {}", i, score, ani);
                                // }

                                if score != std::i32::MIN && -score < min_score_2.unwrap() {
                                    // eprintln!("Set {} -> {}", min_score_2.unwrap(), -score);
                                    min_score_2 = Some(-score);
                                }
                                // eprintln!("{} (asize: {}) Set score {} {} {} {}", i, a.seeds.len(), score, a.score, (1.0 - a.score as f64/cigar.0.len() as f64),  String::from_utf8_lossy(&cigar.0));
                            }
                            // eprintln!("{}", query.len());
                        }
                        None => (),
                    };
                });

            // glidesort::sort_by_key(&mut extension_anchors,|AnchorPair(a1, a2)| {
            //     let s1 = match a1 {
            //         Some(a) => a.score,
            //         None => 0,
            //     };
            //     let s2 = match a2 {
            //         Some(a) => a.score,
            //         None => 0,
            //     };

            //     - ((s1 + s2) as i64)
            // });
        });
        stats.time_alignment += duration;

        let new_score = &extension_anchors[0].0.as_ref().map_or(0, |a| a.score)
            + &extension_anchors[0].1.as_ref().map_or(0, |a| a.score);

        // eprintln!("SCORE  {} -> {}", old_score, new_score);
        // eprintln!(
        //     "MAPQ   {} -> {}",
        //     pseudo_mapq,
        //     StdPairedAnchorMAPQ::anchor_mapq(extension_anchors)
        // );

        //##########

        let best_extended_anchor_pair = extension_anchors.first().unwrap();

        let reference_id = if best_extended_anchor_pair.0.is_some() {
            &best_extended_anchor_pair.0.as_ref().unwrap().reference
        } else {
            &best_extended_anchor_pair.1.as_ref().unwrap().reference
        };

        let reference = &self.db.get_reference(*reference_id as usize).unwrap();

        // if pseudo_mapq == 0 && extension_anchors.len() > 1 {
        //     let one = &extension_anchors[0];
        //     let two = &extension_anchors[1];

        //     let rone = if let Some(a) = &one.0 {
        //         a.reference
        //     } else {
        //         one.1.as_ref().unwrap().reference
        //     };
        //     let rtwo = if let Some(a) = &two.0 {
        //         a.reference
        //     } else {
        //         two.1.as_ref().unwrap().reference
        //     };

        //     if rone == rtwo {
        //         eprintln!("Read at fault: {}", String::from_utf8_lossy(rec_fwd.head()));
        //         eprintln!(
        //             "{}\n\nBest Alignment >>\n{}\n{}\n\nSecond Best Alignment>>\n{}{}\n",
        //             pseudo_mapq,
        //             one.0.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             one.1.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             two.0.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             two.1.as_ref().map_or("None".to_string(), |a| a.to_string())
        //         );

        //         if extension_anchors.len() > 2 {
        //             let three = &extension_anchors[2];

        //             let rthree = if let Some(a) = &three.0 {
        //                 a.reference
        //             } else {
        //                 three.1.as_ref().unwrap().reference
        //             };
        //             eprintln!(
        //                 "\nThird Best Alignment >>\n{:?}\n{:?}\n\n",
        //                 three.0, three.1
        //             );
        //         }

        //         let mut input = String::new();
        //         std::io::stdin()
        //             .read_line(&mut input)
        //             .expect("error: unable to read user input");
        //     }

        //     // assert!(rone != rtwo);
        // }

        // ---- Reactivate or put in different piece of code vv

        // let valid_fwd = best_extended_anchor_pair.0.as_ref().map(|a| {
        //     a.validate_seeds(
        //         if a.forward {
        //             rec_fwd.seq()
        //         } else {
        //             self.rec_fwd_revc.seq()
        //         },
        //         reference,
        //     )
        // });

        // log::trace!("MAPQ: {}", pseudo_mapq);

        // let valid_rev = best_extended_anchor_pair.1.as_ref().map(|a| {
        //     a.validate_seeds(
        //         if a.forward {
        //             rec_rev.seq()
        //         } else {
        //             self.rec_rev_revc.seq()
        //         },
        //         reference,
        //     )
        // });
        // let valid = valid_fwd.unwrap_or(true) && valid_rev.unwrap_or(true);

        // ---- Reactivate or put in different piece of code ^^

        // if anchor_pair.0.is_some() {
        //     let a = anchor_pair.0.as_ref().unwrap();
        //     let query = if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
        //     let (qr, rr) = a.whole(query.len(), reference.len());
        //     a.visualize_alignment(query, reference);
        //     eprintln!("Hammingo: {}\n{}\n{}", a.hamming(query, reference), (&query[qr]).ts(),  (&reference[rr]).ts());

        // }
        // if anchor_pair.1.is_some() {
        //     let a = anchor_pair.1.as_ref().unwrap();
        //     let query = if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };
        //     a.visualize_alignment(query, reference);
        //     let (qr, rr) = a.whole(query.len(), reference.len());
        //     eprintln!("Hammingo: {}\n{}\n{}", a.hamming(query, reference), (&query[qr]).ts(),  (&reference[rr]).ts());
        // }

        // if !valid {
        //     eprintln!("Incidence\n{:?} -> {:?}\n{:?} -> {:?}", valid_fwd, anchor_pair.0, valid_rev, anchor_pair.1)
        // }

        // if anchor_pair.0.is_some() && anchor_pair.1.is_some() {
        //     let normal = anchor_pair.0.as_ref().unwrap().forward;
        //     if normal {
        //         eprintln!("{}\n{}", String::from_utf8_lossy(rec_fwd.seq()), String::from_utf8_lossy(self.rec_rev_revc.seq()));
        //     } else {
        //         eprintln!("{}\n{}", String::from_utf8_lossy(self.rec_fwd_revc.seq()), String::from_utf8_lossy(rec_rev.seq()));
        //     }
        // }

        // let before_ref = best_before.reference();
        // let after_ref = best_after.reference();

        // let before_correct = correct(rec_fwd.head(), before_ref, self.db);
        // let after_correct = correct(rec_fwd.head(), after_ref, self.db);

        // if before_correct && !after_correct {
        //     // eprintln!("Anchors... {}", anchors.len());

        //     // eprintln!("Best before: {:?}", best_before);
        //     // eprintln!("Best after:  {:?}", best_after);
        //     anchors.iter().for_each(|AnchorPair(a1, a2)| {
        //         let reference = match a1 {
        //             Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
        //             None => &self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
        //         };

        //         let hamming1 = match a1 {
        //             Some(a) => {
        //                 let query = if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
        //                 if query.len() == 0 { 0 } else {
        //                     query.len() as u64 - a.hamming(query, reference)
        //                 }
        //             },
        //             None => 0,
        //         };
        //         let hamming2 = match a2 {
        //             Some(a) => {
        //                 let query = if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };
        //                 if query.len() == 0 { 0 } else {
        //                     query.len() as u64 - a.hamming(query, reference)
        //                 }
        //             },
        //             None => 0,
        //         };
        //         let score1 = match a1 {
        //             Some(a) => StdAnchorScore::score(a),
        //             None => 0,
        //         };
        //         let score2 = match a2 {
        //             Some(a) => StdAnchorScore::score(a),
        //             None => 0,
        //         };
        //         eprintln!("{}", - ((hamming1 + hamming2) as i64));
        //         eprintln!("anchor {} -- {} {} .....   {} -> {} _____ {} ->{}", AnchorPair(a1.clone(), a2.clone()).reference(),  score2+score1, hamming1+hamming2, score1, hamming1, score2, hamming2);
        //     });

        //     let mut name = String::new();
        //     std::io::stdin().read_line(&mut name).expect("Read line failed.");
        // }

        let mut print_reads = false;
        if best_extended_anchor_pair.0.is_some() {
            let best = best_extended_anchor_pair.0.as_ref().unwrap();
            let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
            let reference = &self.db.get_reference(best.reference as usize).unwrap();
            let query = if best.forward {
                rec_fwd.seq()
            } else {
                self.rec_fwd_revc.seq()
            };
            let hamming = best.hamming(query, reference);

            // let (qr, rr) = best.whole(query.len(), reference.len());

            // let (duration, (score, cigar)) = time(|| self.align.align(&query[qr], &reference[rr]));
            // stats.time_alignment += duration;

            // let (qr, rr) = best.whole(query.len(), reference.len());

            // let hamming = score / -4;

            // eprintln!("----{:?}\n Score {}, Hamming {}, cigar: {}", valid, (score / -4), hamming, cigar);
            // eprintln!("{:?} {:?}", best.seeds.first(), best.whole(query.len(), reference.len()));
            // eprintln!("{} {}", ref_string, String::from_utf8_lossy(rec_fwd.head()));
            // eprintln!("{}\n{}", String::from_utf8_lossy(&query[qr]), String::from_utf8_lossy(&reference[rr]));

            // let best_corelen = best.core_matches() - best.mismatches as usize - best.indels();
            // let second_best_corelen = if anchors.len() > 1 {
            //     let second_best = anchors.get(1).unwrap();
            //     second_best.core_matches() - second_best.mismatches as usize - second_best.indels()
            // } else { 0 };

            if GOLDSTD_EVAL {
                let correct = evaluate::evaluate(
                    stats.gold_std_evaluation.as_mut().unwrap(),
                    ref_string,
                    pseudo_mapq as u64,
                    &rec_fwd,
                    self.db,
                    true,
                );
                print_reads |= !correct
            }

            if self.options.args.debug {
                let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())]
                    == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];

                if !correct {
                    eprintln!("\n\nIncorrect fwd:");
                    eprintln!("{}", String::from_utf8_lossy(rec_fwd.head()));
                    extension_anchors.print();
                    eprintln!("\nFrom seeds:");
                    eprintln!("\nForward Seeds {}", seeds_fwd.len());
                    for seed in seeds_fwd {
                        eprintln!("\t{}", seed);
                    }
                    eprintln!("\nReverse Seeds {}", seeds_rev.len());
                    for seed in seeds_rev {
                        let seed_ref = self.db.get_rname(seed.rval as usize).unwrap();
                        let seed_correct = &seed_ref.as_bytes()
                            [..min(seed_ref.len(), rec_fwd.head().len())]
                            == &rec_fwd.head()[..min(seed_ref.len(), rec_fwd.head().len())];

                        eprintln!(
                            "\t{} -- {} -- {}",
                            seed,
                            self.db.get_rname(seed.rval as usize).unwrap(),
                            seed_correct
                        );
                    }
                }
            }

            if self.output.has_a() {
                self.output.a.as_mut().unwrap().write(
                    &String::from_utf8_lossy(rec_fwd.head()),
                    rec_fwd.seq().len(),
                    best.seeds.first().unwrap().qbegin() as i32,
                    best.seeds.last().unwrap().qend() as i32,
                    best.forward,
                    ref_string,
                    reference.len(),
                    best.seeds.first().unwrap().rbegin() as i32,
                    best.seeds.last().unwrap().rend() as i32,
                    (query.len() - hamming as usize) as u32,
                    0,
                    pseudo_mapq,
                );
            }
        } else {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().count_fn();
            }
        }

        if best_extended_anchor_pair.1.is_some() {
            let best = best_extended_anchor_pair.1.as_ref().unwrap();
            let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
            let reference = &self.db.get_reference(best.reference as usize).unwrap();
            let query = if best.forward {
                rec_rev.seq()
            } else {
                self.rec_rev_revc.seq()
            };

            let hamming = best.hamming(query, reference);

            // let (qr, rr) = best.whole(query.len(), reference.len());

            // let (duration, (score, cigar)) = time(|| self.align.align(&query[qr], &reference[rr]));
            // stats.time_alignment += duration;

            // let (qr, rr) = best.whole(query.len(), reference.len());

            // let hamming = score / -4;

            if GOLDSTD_EVAL {
                let correct = evaluate::evaluate(
                    stats.gold_std_evaluation.as_mut().unwrap(),
                    ref_string,
                    pseudo_mapq as u64,
                    &rec_rev,
                    self.db,
                    false,
                );
                print_reads |= !correct;
            }

            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())]
                == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];

            if self.options.args.debug {
                if !correct {
                    eprintln!("\n\nIncorrect Rev:");
                    eprintln!("{}", String::from_utf8_lossy(rec_rev.head()));
                    extension_anchors.print();
                    eprintln!("\nFrom seeds:");
                    eprintln!("\nForward Seeds {}", seeds_fwd.len());
                    for seed in seeds_fwd {
                        let seed_ref = self.db.get_rname(seed.rval as usize).unwrap();
                        let seed_correct = &seed_ref.as_bytes()
                            [..min(seed_ref.len(), rec_rev.head().len())]
                            == &rec_rev.head()[..min(seed_ref.len(), rec_rev.head().len())];

                        eprintln!(
                            "\t{} -- {} -- {}",
                            seed,
                            self.db.get_rname(seed.rval as usize).unwrap(),
                            seed_correct
                        );
                    }
                    eprintln!("\nReverse Seeds {}", seeds_rev.len());
                    for seed in seeds_rev {
                        eprintln!("\t{}", seed);
                    }
                }
            }

            if self.output.has_a() {
                self.output.a.as_mut().unwrap().write(
                    &String::from_utf8_lossy(rec_rev.head()),
                    rec_rev.seq().len(),
                    best.seeds.first().unwrap().qbegin() as i32,
                    best.seeds.last().unwrap().qend() as i32,
                    best.forward,
                    ref_string,
                    reference.len(),
                    best.seeds.first().unwrap().rbegin() as i32,
                    best.seeds.last().unwrap().rend() as i32,
                    (query.len() - hamming as usize) as u32,
                    0,
                    pseudo_mapq,
                );
            }
        } else {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().count_fn();
            }
        }

        if GOLDSTD_EVAL {
            if print_reads && pseudo_mapq != 0 {
                let eval = stats.gold_std_evaluation.as_mut().unwrap();
                if let Some(or1) = eval.output_read1.as_mut() {
                    assert!(rec_fwd.head().ends_with(b"1"));
                    or1.write(format!("@{}\n", rec_fwd.to_string()));
                }
                if let Some(or2) = eval.output_read2.as_mut() {
                    assert!(rec_rev.head().ends_with(b"2"));
                    or2.write(format!("@{}\n", rec_rev.to_string()));
                }
            }
        }

        // stats.time_reverse_complement += duration;
        // stats.time_anchor_sorting += duration;
        // let (duration, _) = time(|| {
        //     rec_rev.reverse_complement(&mut self.rec_rev_revc);
        // });
        // stats.time_reverse_complement += duration;
        // stats.time_anchor_sorting += duration;

        // let best = anchors.first().unwrap();
        // let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
        // let reference = &self.db.get_reference(best.reference as usize).unwrap();

        // let best_corelen = best.core_matches() - best.mismatches as usize - best.indels();
        // let second_best_corelen = if anchors.len() > 1 {
        //     let second_best = anchors.get(1).unwrap();
        //     second_best.core_matches() - second_best.mismatches as usize - second_best.indels()
        // } else { 0 };

        // let pseudo_mapq = best_corelen - second_best_corelen;

        // Compile time switch
        // if  {
        // GOLDSTD_EVAL
        //     // @NC_009436.1_4088855_4089351_1:2:0_1:5:2_2/1

        //     let header_str = String::from_utf8_lossy(rec_fwd.head());
        //     let first_part_a = header_str.split('-').next().unwrap_or("");
        //     let first_part_b = header_str.splitn(3, '_').take(2).collect::<Vec<&str>>().join("_");
        //     let mut true_id = *self.db.get_rid(first_part_a).unwrap_or(&0);

        //     if true_id == 0 {
        //         true_id = *self.db.get_rid(&first_part_b).unwrap_or(&0);
        //     }

        //     if true_id == 0 {
        //         panic!("True id is {}", true_id);
        //     }

        //     let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];
        //     // eprintln!("{}\t{}\t{}\t{}", ref_string, header_str, correct, pseudo_mapq);

        //     stats.gold_std_evaluation.as_mut().unwrap().add(correct, pseudo_mapq as u64);

        // }

        // self.output_paf.as_mut().unwrap().write(
        //     &String::from_utf8_lossy(rec_fwd.head()),
        //     rec_fwd.seq().len(),
        //     best.seeds.first().unwrap().qbegin() as i32,
        //     best.seeds.last().unwrap().qend() as i32,
        //     best.forward,
        //     ref_string,
        //     reference.len(),
        //     best.seeds.first().unwrap().rbegin() as i32,
        //     best.seeds.last().unwrap().rend() as i32,
        //     best.seed_count,
        //     0,
        //     pseudo_mapq as u8);
    }
}
