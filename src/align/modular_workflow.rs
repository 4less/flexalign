use std::{cmp::min, os::linux::raw::stat};

use bioreader::sequence::fastq_record::{print_color_qualities, OwnedFastqRecord, RefFastqRecord};
use colored::Colorize;
use kmerrs::{consecutive::kmer::{Kmer, KmerIter}, minimizer::context_free::Minimizer};

use crate::{align::{common::{AnchorScore, Print, StdAnchorScore}, data_structures::ToString}, database::common::FlexalignDatabase, flexalign::time, options::Options, GOLDSTD_EVAL};

use super::{common::{is_alignment_valid, print_alignment, Align, AnchorExtractor, AnchorPair, Heuristic, KmerExtractor, Or, PAFOutput, PairedAnchorExtractor, PairedAnchorMAPQ, PairedAnchorSorter, RangeExtractor, SAMOutput, SeedExtractor, StdPairedAnchorMAPQ}, process::{alignment::ani_abort_score, evaluate::{self, correct, get_id_from_header}, output::StdPAFOutput}, stats::Stats};


#[derive(Clone)]
pub struct Modular<
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor<C>,
    RE: RangeExtractor<C, F>,
    SE: SeedExtractor::<F>,
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
    KE: KmerExtractor::<C>,
    RE: RangeExtractor::<C, F>,
    SE: SeedExtractor::<F>,
    AE: AnchorExtractor,
    PO: PAFOutput,
    SO: SAMOutput,
    D: FlexalignDatabase
    > Modular<'a, C, F, KE, RE, SE, AE, PO, SO, D> { //RE, SE, 
    pub fn run(
        &mut self,
        rec: &RefFastqRecord,
        stats: &mut Stats) -> ()
    {
        stats.reads_processed += 1;

        let (duration, kmers) = time(|| {
            self.kmer_extractor.generate(rec, stats)
        });
        stats.time_get_kmers += duration;

        let (duration, ranges) = time(|| {
            self.range_extractor.generate(kmers, stats)
        });
        stats.time_get_ranges += duration;

        let (duration, seeds) = time(|| {
            self.seed_extractor.generate(ranges, stats)
        });
        stats.time_range_header += duration;
        stats.seeds += seeds.len();

        let (duration, anchors) = time(|| {
            self.anchor_extractor.generate(seeds, rec.seq().len(), stats)
        });
        stats.time_range_header += duration;
        stats.anchors += anchors.len();

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return
        }

        let (duration, _) = time(|| {
            anchors.sort_unstable_by_key(|a| {
                - ((a.core_matches() - a.mismatches as usize - a.indels()) as i64)
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
        } else { 0 };

        let pseudo_mapq = best_corelen - second_best_corelen;

        // Compile time switch
        if GOLDSTD_EVAL {

            // @NC_009436.1_4088855_4089351_1:2:0_1:5:2_2/1

            // let header_str = String::from_utf8_lossy(rec.head());
            // let first_part_a = header_str.split('-').next().unwrap_or("");
            // let first_part_b = header_str.splitn(3, '_').take(2).collect::<Vec<&str>>().join("_");
            // let mut true_id = *self.db.get_rid(first_part_a).unwrap_or(&0);

            // if true_id == 0 {
            //     true_id = *self.db.get_rid(&first_part_b).unwrap_or(&0);
            // }

            // if true_id == 0 {
            //     panic!("True id is {}", true_id);
            // }


            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec.head().len())] == &rec.head()[..min(ref_string.len(), rec.head().len())];
            // eprintln!("{}\t{}\t{}\t{}", ref_string, header_str, correct, pseudo_mapq);


            stats.gold_std_evaluation.as_mut().unwrap().add(correct, pseudo_mapq as u64);
            
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
                pseudo_mapq as u8);
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
    SE: SeedExtractor::<F>,
    AE: PairedAnchorExtractor,
    AS: PairedAnchorSorter,
    PO: PAFOutput,
    SO: SAMOutput,
    A: Align + Heuristic + Send,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,
    pub db: &'a D,
    pub kmer_extractor_fwd: KE,
    pub kmer_extractor_rev: KE,
    pub range_extractor_fwd: RE,
    pub range_extractor_rev: RE,
    pub seed_extractor_fwd: SE,
    pub seed_extractor_rev: SE,
    pub anchor_extractor: AE,
    pub anchor_sorter: AS,

    pub align: A,

    pub output: Or<PO, SO>,

    pub rec_fwd_revc: OwnedFastqRecord,
    pub rec_rev_revc: OwnedFastqRecord,
}

impl<   
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor::<C>,
    RE: RangeExtractor::<C, F>,
    SE: SeedExtractor::<F>,
    AE: PairedAnchorExtractor,
    AS: PairedAnchorSorter,
    PO: PAFOutput,
    SO: SAMOutput,
    A: Align + Heuristic + Send,
    D: FlexalignDatabase
    > ModularPE<'a, C, F, KE, RE, SE, AE, AS, PO, SO, A, D> { //RE, SE, 
    pub fn run(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        stats: &mut Stats) -> ()
    {
        stats.reads_processed += 2;

        // Extract minimizer
        let (duration, kmers_fwd) = time(|| {
            self.kmer_extractor_fwd.generate(rec_fwd, stats)
        });
        stats.time_get_kmers += duration;
        let (duration, kmers_rev) = time(|| {
            self.kmer_extractor_rev.generate(rec_rev, stats)
        });
        stats.time_get_kmers += duration;


        // Get ranges from minimizers
        let (duration, ranges_fwd) = time(|| {
            self.range_extractor_fwd.generate(kmers_fwd, stats)
        });
        stats.time_get_ranges += duration;
        let (duration, ranges_rev) = time(|| {
            self.range_extractor_rev.generate(kmers_rev, stats)
        });
        stats.time_get_ranges += duration;


        // Get Seeds from ranges
        let (duration, seeds_fwd) = time(|| {
            self.seed_extractor_fwd.generate(ranges_fwd, stats)
        });
        stats.time_range_header += duration;
        stats.seeds += seeds_fwd.len();
        let (duration, seeds_rev) = time(|| {
            self.seed_extractor_rev.generate(ranges_rev, stats)
        });
        stats.time_range_header += duration;
        stats.seeds += seeds_rev.len();

        // eprintln!("Header {} ... \nID {}", String::from_utf8_lossy(rec_fwd.head()), get_id_from_header(&String::from_utf8_lossy(rec_fwd.head()), self.db));
        let (duration, mut anchors) = time(|| {
            self.anchor_extractor.generate(seeds_fwd, seeds_rev, rec_fwd.seq().len(), rec_rev.seq().len(), stats)
        });
        stats.time_get_anchors += duration;
        stats.anchors += anchors.len();

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return
        }

        // eprintln!("Read: {}", String::from_utf8_lossy(rec_fwd.head()));
        let best_before = anchors.first().as_mut().unwrap().clone();

        // Now here starts the reference-based portion of the algorithm. Before, no sequence comparison
        // Between query and reference is done
        let (duration, _) = time(|| {
            rec_fwd.reverse_complement(&mut self.rec_fwd_revc);
            rec_rev.reverse_complement(&mut self.rec_rev_revc);
        });
        stats.time_reverse_complement += duration;


        let anchors_len = anchors.len();
        let max_hamming = 10;



        // Assumes sorted anchors !!
        let mut extension_anchors = &mut anchors[0..min(self.options.args.extend_top_x, anchors_len)];

        
        extension_anchors.iter_mut().enumerate().for_each(|(i, (AnchorPair(a1, a2)))| {
            match a1 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("Y 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("Y 2  {}", a);
                    }
                }, _ => {},
            }
        });

        // Assumes valid anchor seeds!!
        let (duration, _) = time(|| {
            self.anchor_sorter.sort(extension_anchors, rec_fwd, &self.rec_fwd_revc, rec_rev, &self.rec_rev_revc, stats);
        });
        stats.time_extend_anchors += duration;

        
        extension_anchors.iter_mut().enumerate().for_each(|(i, (AnchorPair(a1, a2)))| {
            match a1 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("Z 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("Z 2  {}", a);
                    }
                }, _ => {},
            }
        });


        // Assumes sorted anchors !!
        let anchors_len: usize = extension_anchors.len();
        let alignment_anchors = &mut extension_anchors[0..min(self.options.args.align_top_y, anchors_len)];

        let (duration, _) = time(|| {
            let mut min_score_1 = None;
            let mut min_score_2 = None;

            alignment_anchors.iter_mut().enumerate().for_each(|(i, (AnchorPair(a1, a2)))| {
                let reference = match a1 {
                    Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
                    None => &self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
                };

                match a1 {
                    Some(a) => {
                        let query = if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
                        if query.len() == 0 { 
                            a.score = 0i32;
                        } else {

                            if min_score_1.is_none() {
                                min_score_1 = Some(ani_abort_score(0.5, 4, query.len() as i32).abs());
                            }
                            self.align.set_max_alignment_score(min_score_1.unwrap());
                            // eprintln!("Align max score: {}", min_score_1.unwrap());

                            if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                eprintln!("1  {}", a);
                            }

                            let status = a.smart_align(&mut self.align, query, reference, 10, min_score_1.unwrap());
                            // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_1.unwrap());
                            

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

                            let ani = (1.0 - a.score as f64/a.cigar().0.len() as f64);
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
                    },
                    None => (),
                };
                match a2 {
                    Some(a) => {
                        let query = if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };
                        if query.len() == 0 { 
                            a.score = 0i32;
                        } else {
                            if min_score_2.is_none() {
                                min_score_2 = Some(ani_abort_score(0.5, 4, query.len() as i32).abs());
                            }

                            if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                eprintln!("2  {}", a);
                            }

                            self.align.set_max_alignment_score(min_score_2.unwrap());
                            let status = a.smart_align(&mut self.align, query, reference, 10, min_score_2.unwrap());
                            // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_2.unwrap());
                            
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

                            let ani = (1.0 - a.score as f64/a.cigar().0.len() as f64);
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
                    },
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

//#######################


        // if self.output.has_a() {

        // }


        // for AnchorPair(a1, a2) in extension_anchors.iter() {
        //     let s1 = match a1 {
        //         Some(a) => a.score,
        //         None => 0,
        //     };
        //     let s2 = match a2 {
        //         Some(a) => a.score,
        //         None => 0,
        //     };
        //     eprintln!("Score: {} {}", s1, s2);
        // }


        let best_after = extension_anchors.first().unwrap().clone();


        let pseudo_mapq = StdPairedAnchorMAPQ::anchor_mapq(extension_anchors);
        let anchor_pair = extension_anchors.first().unwrap();
        
        let reference_id = if anchor_pair.0.is_some() { &anchor_pair.0.as_ref().unwrap().reference } else { &anchor_pair.1.as_ref().unwrap().reference };

        let reference = &self.db.get_reference(*reference_id as usize).unwrap();
        
        
        let valid_fwd = anchor_pair.0.as_ref().map(|a| a.validate_seeds(if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() }, reference));
        let valid_rev = anchor_pair.1.as_ref().map(|a| a.validate_seeds(if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() }, reference));
        let valid = valid_fwd.unwrap_or(true) && valid_rev.unwrap_or(true);
        

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

        let valid = valid_fwd.unwrap_or(true) && valid_rev.unwrap_or(true);
        

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

        if anchor_pair.0.is_some() {

            let best = anchor_pair.0.as_ref().unwrap();
            let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
            let reference = &self.db.get_reference(best.reference as usize).unwrap();
            let query = if best.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
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
                evaluate::evaluate(stats.gold_std_evaluation.as_mut().unwrap(), ref_string, pseudo_mapq as u64, &rec_fwd, self.db);
            }

            if self.options.args.debug {
                let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];

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
                        let seed_correct = &seed_ref.as_bytes()[..min(seed_ref.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(seed_ref.len(), rec_fwd.head().len())];
            
                        eprintln!("\t{} -- {} -- {}", seed, self.db.get_rname(seed.rval as usize).unwrap(), seed_correct);
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
                    pseudo_mapq);
            }

        }

                
        if anchor_pair.1.is_some() {
            let best = anchor_pair.1.as_ref().unwrap();
            let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
            let reference = &self.db.get_reference(best.reference as usize).unwrap();
            let query = if best.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };

            let hamming = best.hamming(query, reference);

            // let (qr, rr) = best.whole(query.len(), reference.len());
            
            // let (duration, (score, cigar)) = time(|| self.align.align(&query[qr], &reference[rr]));
            // stats.time_alignment += duration;
            
            // let (qr, rr) = best.whole(query.len(), reference.len());
            
            // let hamming = score / -4;

            if GOLDSTD_EVAL {
                evaluate::evaluate(stats.gold_std_evaluation.as_mut().unwrap(), ref_string, pseudo_mapq as u64, &rec_fwd, self.db);
            }
            
            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];

            if self.options.args.debug {
                if !correct {
                    eprintln!("\n\nIncorrect Rev:");
                    eprintln!("{}", String::from_utf8_lossy(rec_rev.head()));
                    extension_anchors.print();
                    eprintln!("\nFrom seeds:");
                    eprintln!("\nForward Seeds {}", seeds_fwd.len());
                    for seed in seeds_fwd {
                        let seed_ref = self.db.get_rname(seed.rval as usize).unwrap();
                        let seed_correct = &seed_ref.as_bytes()[..min(seed_ref.len(), rec_rev.head().len())] == &rec_rev.head()[..min(seed_ref.len(), rec_rev.head().len())];
            
                        eprintln!("\t{} -- {} -- {}", seed, self.db.get_rname(seed.rval as usize).unwrap(), seed_correct);
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
                    pseudo_mapq);
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