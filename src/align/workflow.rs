
use std::{cmp::{max, min}, io::stdin, mem::swap};

use super::{super::GOLDSTD_EVAL, common::{KmerExtractor}, stats::Stats};
use bioreader::sequence::fastq_record::{print_color_qualities, OwnedFastqRecord, RefFastqRecord};
use flexmap::{values::{VData, VRange}, VD};
use kmerrs::{consecutive::kmer::{Kmer, KmerIter}, minimizer::context_free::Minimizer};
use colored::Colorize;

use crate::{database::common::FlexalignDatabase, flexalign::time, io::output_buffer::OutputBuffer, options::Options};

use super::data_structures::{Anchor, AnchorSeed, Seed};


#[derive(Clone)]
pub struct Standard<
    'a,
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize,
    T: Minimizer,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,

    pub ranges: Vec<(usize, Kmer<F>, VRange<'a>, usize)>,
    pub seeds: Vec<Seed>,
    pub anchors: Vec<Anchor>,
    pub indices: Vec<usize>,
    pub other_indices: Vec<usize>,
    pub minimizer: T,
    pub db: &'a D,

    pub output_buffer: Vec<u8>,
    pub output_buffer_threshold: usize,

    pub ob: OutputBuffer,

    pub rev_rec: OwnedFastqRecord,

    pub count_difficult_anchors: usize,

    // pub range_iteration_order: Vec<(usize, usize)>,
    // pub seeds: Vec<(usize, VCell)>,
}

impl<
        'a,
        const K: usize,
        const C: usize,
        const F: usize,
        const S: usize,
        const L: usize,
        const HEADER_THRESHOLD: usize,
        T: Minimizer,
        D: FlexalignDatabase
    > Standard<'a, K, C, F, S, L, HEADER_THRESHOLD, T, D>
{
    pub fn new(db: &'a D, minimizer: T, options: &'a Options, output: OutputBuffer) -> Self {
        Self {
            ranges: Vec::new(),
            // range_iteration_order: Vec::new(),
            seeds: Vec::new(),
            anchors: Vec::new(),

            indices: Vec::new(),
            other_indices: Vec::new(),

            ob: output,
            output_buffer: Vec::new(),
            output_buffer_threshold: 1_000_000,
            
            minimizer,
            db: db,

            options: options,

            rev_rec: OwnedFastqRecord::new(),

            count_difficult_anchors: 0,
        }
    }

    pub fn run(
        &mut self,
        rec: &RefFastqRecord,
        stats: &mut Stats) -> ()
    {
        stats.reads_processed += 1;

        self.ranges.clear();
        self.seeds.clear();
        self.anchors.clear();
        
        self.get_ranges(rec, stats);
        self.get_seeds(stats);
        // self.get_anchors(stats);

        // eprintln!("--Seeds--");
        // for seed in &self.seeds {
        //     eprintln!("{}", seed.to_string());
        // }
        // eprintln!("---------");

        let (duration, _) = time(|| {
            self.seeds_to_anchors(stats, rec.seq().len());
        });
        stats.time_get_anchors += duration;

        let (duration, _) = time(|| {
            self.anchors.sort_unstable_by_key(|a| {
                - ((a.core_matches() - a.mismatches as usize - a.indels()) as i64)
            });
        });
        stats.time_anchor_sorting += duration;

        
        if self.anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return
        }

        let (duration, _) = time(|| {
            rec.reverse_complement(&mut self.rev_rec);
        });
        stats.time_reverse_complement += duration;


        let extend_anchors = true;
        if extend_anchors {
            let (duration, _) = time(|| {
                self.extend_anchors(rec);
            });
            stats.time_extend_anchors += duration;
        }

        // let (duration, _) = time(|| {
        //     for anchor in &self.anchors {
        //         let sane = self.sanity_check_anchor(&anchor, rec, &self.rev_rec);
        //         if !sane { 
        //             self.count_difficult_anchors += 1;
        //             // eprintln!("Difficult: {}", self.count_difficult_anchors);
        //         };
        //     }
        // });
        // stats.time_checking_anchors += duration;
        // eprintln!("{}", self.count_difficult_anchors as f64 / stats.anchors as f64);


        let best = self.anchors.first().unwrap();
        let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
        let reference = &self.db.get_reference(best.reference as usize).unwrap();

        let best_corelen = best.core_matches() - best.mismatches as usize - best.indels();
        let second_best_corelen = if self.anchors.len() > 1 {
            let second_best = self.anchors.get(1).unwrap();
            second_best.core_matches() - second_best.mismatches as usize - second_best.indels()
        } else { 0 };

        let pseudo_mapq = best_corelen - second_best_corelen;

        // Compile time switch
        if GOLDSTD_EVAL {

            // @NC_009436.1_4088855_4089351_1:2:0_1:5:2_2/1

            let header_str = String::from_utf8_lossy(rec.head());
            let first_part_a = header_str.split('-').next().unwrap_or("");
            let first_part_b = header_str.splitn(3, '_').take(2).collect::<Vec<&str>>().join("_");
            let mut true_id = *self.db.get_rid(first_part_a).unwrap_or(&0);

            if true_id == 0 {
                true_id = *self.db.get_rid(&first_part_b).unwrap_or(&0);
            }

            if true_id == 0 {
                panic!("True id is {}", true_id);
            }


            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec.head().len())] == &rec.head()[..min(ref_string.len(), rec.head().len())];
            // eprintln!("{}\t{}\t{}\t{}", ref_string, header_str, correct, pseudo_mapq);


            stats.gold_std_evaluation.as_mut().unwrap().add(correct, pseudo_mapq as u64);

            if !correct {
                let any_seed_match = self.anchors.iter().any(|a| a.reference == true_id as u64);
                let any_anchor_match = self.anchors.iter().any(|a| a.reference == true_id as u64);
                
                if any_anchor_match {

                }

                if any_seed_match && !any_anchor_match {

                    eprintln!("\n\n_______{}\t\t{}\t{}\t{}\t{}\t{}\t{}", true_id, any_seed_match, any_anchor_match, ref_string, header_str, correct, pseudo_mapq);

                    eprintln!("\n--------------------------------------- Anchor {}", self.anchors.len());
                    for (i,anchor) in self.anchors.iter().enumerate() {
                        let anchor_ref = &self.db.get_rname(anchor.reference as usize).unwrap();
                        let correct = &anchor_ref.as_bytes()[..min(anchor_ref.len(), rec.head().len())] == &rec.head()[..min(anchor_ref.len(), rec.head().len())];
                        
                        eprintln!("\n{}  {}  ---  {}   /   {} ___________________sane? {}",
                            if correct { ">>>>>".green().bold() } else { "_____".red() },
                                i, anchor_ref, String::from_utf8_lossy(rec.head()), self.sanity_check_anchor(anchor, rec, &self.rev_rec));
                        self.print_debug(rec, &self.rev_rec, anchor);
                    }

                    eprintln!("\n--------------------------------------- Seeds {} (True: {})", self.seeds.len(), true_id);
                    for seed in &self.seeds {
                        eprintln!("{}", seed);
                    }

                    eprintln!("QUALITIES-------------");
                    print_color_qualities(rec.qual(), Some(33));

                    let mut input: String = String::default();
                    stdin().read_line(&mut input).expect("Did not enter a correct string");
                }
            }
        }

        self.ob.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            String::from_utf8_lossy(rec.head()), 
            rec.seq().len(),
            best.seeds.first().unwrap().qbegin(),
            best.seeds.last().unwrap().qend(),
            if best.forward { '+' } else { '-' },
            ref_string,
            reference.len(),
            best.seeds.first().unwrap().rbegin(),
            best.seeds.last().unwrap().rend(),
            best.seed_count, 
            pseudo_mapq));

        
    }


    pub fn print_debug(&self, rec_fwd: &RefFastqRecord, rec_rev: &OwnedFastqRecord, anchor: &Anchor) {
        eprintln!("{}\n{}", anchor.to_string(), self.sanity_check_anchor(anchor, rec_fwd, rec_rev));
    }

    pub fn sanity_check_anchor(&self, anchor: &Anchor, rec_fwd: &RefFastqRecord, rec_rev: &OwnedFastqRecord) -> bool {
        if !anchor.orientation_set {
            return true
        }
        
        let reference = &self.db.get_reference(anchor.reference as usize).unwrap();

        if anchor.seeds.first().unwrap().rpos as usize > reference.len() {
            eprintln!("Seed and anchor are invalid\n{}", anchor.to_string());
            panic!("Not good");
        }

        let mut _sane = true;

        let mut difficult_anchor = false;
        for seed in &anchor.seeds {
            let mut seed_fwd = seed.clone();
            let mut seed_rev = seed.clone();

            if !anchor.forward { 
                let _ = &seed_fwd.reverse(rec_fwd.seq().len()); 
            } else {
                let _ = &seed_rev.reverse(rec_fwd.seq().len());
            };

            let qseed_fwd = &rec_fwd.seq()[seed_fwd.qpos as usize..(seed_fwd.qpos + seed_fwd.length) as usize];
            let qseed_rev = &rec_rev.seq()[seed_rev.qpos as usize..(seed_rev.qpos + seed_rev.length) as usize];
            let qseed_fwd2 = &rec_rev.seq()[seed_fwd.qpos as usize..(seed_fwd.qpos + seed_fwd.length) as usize];
            let qseed_rev2  = &rec_fwd.seq()[seed_rev.qpos as usize..(seed_rev.qpos + seed_rev.length) as usize];
            let rseed = &reference[seed.rpos as usize..(seed.rpos+seed.length as u64) as usize];

            if anchor.forward && qseed_rev == rseed || !anchor.forward && qseed_fwd == rseed {
                difficult_anchor = true;
            }
            
            if qseed_fwd != rseed && qseed_rev != rseed {
                _sane = false;
                panic!("No seed is perfect match with reference. Problem\n{}\n{}\n{}\n{} fwd\n{} ref\n{} rev\n{}\n{}", 
                    seed.to_string(),
                    seed_fwd.to_string(),
                    seed_rev.to_string(),
                    String::from_utf8_lossy(qseed_fwd), 
                    String::from_utf8_lossy(rseed), 
                    String::from_utf8_lossy(qseed_rev), 
                    String::from_utf8_lossy(qseed_fwd2), 
                    String::from_utf8_lossy(qseed_rev2));
            }
        }
        !difficult_anchor
    }

    pub fn extend_anchors(&mut self, rec: &RefFastqRecord) {
        let best = self.anchors.last_mut().unwrap();
        let best_first = best.seeds.first().unwrap();
        let _ref_string = &self.db.get_rname(best.reference as usize).unwrap();

        let reference = &self.db.get_reference(best.reference as usize).unwrap();
        let seq = *reference;
        // let start = if best.ref_pos < 100 { 0 } else {best.ref_pos - 100} as usize;
        // let end: usize = min(seq.len(), best.ref_pos as usize + 100) as usize;
        let start = best_first.rpos as usize + F/2 as usize;
        let end: usize = start + C;
        // let read_slice = &rec.seq()[best.read_pos as usize + F/2..best.read_pos as usize + F/2 as usize + C];
        let ref_seed = &seq[start..end];

        let _read_seed = if best.seed_count == 1 && !best.orientation_set {

            let start_fwd = best_first.qpos as usize;
            let seed_fwd = &rec.seq()[start_fwd..start_fwd + best_first.length as usize];
            let start_rev = rec.seq().len() - best_first.length as usize - best_first.qpos as usize;
            let seed_rev = &self.rev_rec.seq()[start_rev..start_rev + best_first.length as usize];
            if ref_seed == seed_fwd {
                best.forward = true;
                seed_fwd
            } else {
                best.forward = false;
                seed_rev
            }
            
        } else if best.forward {
            let start = best_first.qpos as usize;
            &rec.seq()[start..start + best_first.length as usize]
        } else {
            let start = best_first.qpos as usize;
            &self.rev_rec.seq()[start..start + best_first.length as usize]
        };
        
    }

    pub fn seed_group_indices(&self) -> (Vec<(u32, u32)>, usize) {
        let mut last_idx = 0;
        let mut groups = Vec::new();
        let mut max_size = 0;

        for i in 1..self.seeds.len() {
            let prev = &self.seeds[i-1];
            let next = &self.seeds[i];
            if prev.rval != next.rval {
                groups.push((last_idx as u32, i as u32));
                max_size = max(max_size, i - last_idx);
                last_idx = i;
            }
        }

        groups.push((last_idx as u32, self.seeds.len() as u32));
        max_size = max(max_size, self.seeds.len() - last_idx);

        
        let acc = groups.iter().fold(0, |acc, (start, end)| { acc + (end-start) });
        if self.seeds.len() != acc as usize {
            panic!("{} {}", acc, self.seeds.len());
        }


        (groups, max_size)
    }

    pub fn get_best(&self) -> Option<Anchor> {
        self.anchors.last().cloned()
    }

    pub fn seeds_to_anchors(
        &mut self, stats: &mut Stats, read_length: usize
    ) -> () {
        
        let (groups, max_size) = self.seed_group_indices();
        // eprintln!(" ------------------ BARRIER --------------------- {}", groups.len());

        // for (start, end) in &groups {
        //     eprintln!("{} {}", start, end);
        // }

        // eprintln!("Max size {}", max_size);

        let skip_threshold = max_size as i32 - 10;

        stats.anchors += groups.len();
        for (start, end) in groups {
            // eprintln!("{} < {} == {} ({})", end-start, max_size - 5, (end - start) < (max_size as u32 - 5), max_size);
            if ((end - start) as i32) < skip_threshold && (end-start) <= self.ranges.len() as u32 { 
                // eprintln!("Skip {} {}, {}, {},  {}", start, end, end-start, self.options.args.ranges, skip_threshold);
                continue 
            };

            self.group_into_anchor(start as usize, end as usize, read_length);
        }
    }

    pub fn group_into_anchor(&mut self, start: usize, end: usize, read_length: usize) {
        let seeds = &mut self.seeds[start..end];
        
        // eprintln!("-- Seed group -- {} ... {} - {}", seeds.first().unwrap().rval, start, end);
        // for seed in seeds.iter() {
        //     eprintln!("{} --- {}", seed, seed.reverse(read_length));
        // }

        let _anchor_group_index = self.anchors.len();

        if !seeds.is_empty() {
            // Group by exact offset. If there are seeds left, distribute them onto 

            self.indices.clear();
            self.indices.extend(0..seeds.len());
            self.other_indices.clear();

            while !self.indices.is_empty() {
                let first = seeds.get(*self.indices.first().unwrap()).unwrap();
                
                // Set first seed as an anchor.
                let mut a = Anchor::from_seed(first);

                // eprintln!("Starting with anchor: {}", first.to_string());
                
                let mut offset = None;
                let mut forward = None;

                // Find seeds that have 0 indels with respect to the first anchor
                for index in self.indices.iter().skip(1) {
                    let next: &Seed = &seeds[*index];

                    let (offset_first, fwd, indel_first) = first.closest_offset(&next, read_length);

                    if indel_first == 0 && (forward.is_none() || forward.unwrap() == fwd) {
                        if offset.is_none() {
                            // eprintln!("Establish offset: {} {} {}", offset_first, fwd, indel_first);
                            // eprintln!("Between: {} {}", );
                            offset = Some(offset_first);
                            a.set_forward(fwd, read_length);
                            forward = Some(fwd);
                        }

                        // eprintln!("Add seed: {} ... {} == {} && {} == 0", next.to_string(), offset.unwrap(), offset_first, indel_first);
                        let read_dummy = String::from_utf8(vec![b'X'; read_length as usize]).unwrap();
                        // eprintln!("{}", next.to_visual_string(read_dummy.as_bytes()));
                        a.add_seed(&next, read_length as u32);

                    } else {
                        self.other_indices.push(*index);
                    }
                }

                // eprintln!("Final anchor: {}", a);
                self.anchors.push(a);

                self.indices.clear();
                swap(&mut self.other_indices, &mut self.indices);
            }

            // let mut anchors = &self.anchors[anchor_group_index..];
            
        }
        // eprintln!("END -- Seed group -- {}", seeds.first().unwrap().rval );


    }

    pub fn retrieve_seeds(
            &mut self, 
            max_best_flex: usize, 
            ranges: usize, 
            stats: &mut Stats) -> (usize, usize) {

        let mut matches = 0;
        let mut discarded_max_flex_count = 0;
        for (qpos, flex, range, _range_size) in &self.ranges {
            match range.header {
                Some(headers) => {
                    let mut min_dist = u32::MAX;
                    let mut count = 0;
                    for header in headers {
                        let dist = header.dist(flex.0 as u32);
                        if dist < min_dist { min_dist = dist; count = 0; }
                        if dist == min_dist { count += 1};
                    }
                    
                    let take = count <= max_best_flex;
                    // eprintln!("{} Range count = {}/{} < {}", if take { "X".green() } else { "O".red() }, count, range.positions.len(), self.options.args.max_best_flex);
                    if !take {
                        discarded_max_flex_count += 1;
                        continue;
                    }

                    // for header in headers {
                    //     let dist = header.dist(flex.0 as u32);
                    //     if dist < min_dist { min_dist = dist }
                    // }


                    // eprintln!("Header------");
                    for (index, header) in headers.iter().enumerate() {
                        let dist = header.dist(flex.0 as u32);
                        if dist == min_dist {
                            // self.seeds.push((*pos, range.positions[index].clone()))
                            let (value, rpos) = VD::get(range.positions[index].0);
                            self.seeds.push(Seed::from_flexmer::<K,C,F>(*qpos, rpos, value, dist));
                            // eprintln!("H Push: {} {} {}", self.seeds.last().unwrap().to_string(), *pos, rpos);
                        }
                    }
                },
                None => {
                    for cell in range.positions {
                        // self.seeds.push((*pos, cell.clone()));
                        let (value, rpos) = VD::get(cell.0);
                        self.seeds.push(Seed::from_coremer::<K,C,F>(*qpos, rpos, value));
                        // eprintln!("N Push: {} {} {}", self.seeds.last().unwrap().to_string(), *pos, rpos);
                    }
                },
            };
            matches += 1;
            if matches >= ranges { break }
        }
        (matches, discarded_max_flex_count)
    }

    pub fn get_seeds(
        &mut self, stats: &mut Stats
    ) -> () {
        let (duration, _) = time(|| {
            // eprintln!("----------------- Get Ranges....");
            let (ranges, discarded_max_flex_count) = self.retrieve_seeds(
                self.options.args.max_best_flex,
                self.options.args.ranges as usize,
                stats
            );
            // eprintln!("Ranges: {}/{}, Discarded: {} ({})", ranges, self.ranges.len(), discarded_max_flex_count, self.options.args.max_best_flex);
    
            if ranges < self.options.args.min_ranges && discarded_max_flex_count > 0  {
                // eprintln!("----------------- Recover Ranges....");
                let old_ranges = ranges;
                let (ranges, discarded_max_flex_count) = self.retrieve_seeds(
                    128,
                    self.options.args.ranges as usize,
                    stats
                );
                // eprintln!("{} -> {} (Still discarded: {})", old_ranges, ranges, discarded_max_flex_count);
            }
    
        });
        
        stats.time_range_header += duration;
        stats.seeds += self.seeds.len();

        let (duration, _) = time(|| {
            glidesort::sort_by_key(&mut self.seeds, |seed: &Seed| {
                (seed.rval, seed.rpos)//, seed.length) // seed.offset
            });

            // self.seeds.sort_unstable_by_key(|seed| {
            //     (seed.rval, seed.rpos)//, seed.length) // seed.offset
            // })
        });
        stats.time_seed_sorting += duration;

        
    }




    pub fn get_ranges(
        &mut self,
        rec: &RefFastqRecord,
        stats: &mut Stats,
    ) -> () {
        let iter = KmerIter::<K, true>::new(rec.seq());
        self.ranges.clear();
        for (pos, kmer_fwd, kmer_rev) in iter {
            stats.kmers_processed += 1;

            let cmer_fwd = kmer_fwd.middle::<C>();
            let cmer_rev = kmer_rev.middle::<C>();
            let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
            let cmer = min(cmer_fwd, cmer_rev);

            if !self.minimizer.is_minimizer(cmer.0) {
                continue;
            };

            stats.minimizer += 1;

            let fmer = kmer.flanks::<F>();

            let (duration, result) = time(|| self.db.get_vrange(cmer.0));
            stats.time_get_ranges += duration;

            let range: VRange<'a> = match result {
                Some(range) => range,
                None => continue,
            };
            let range_len = (&range).positions.len();
            self.ranges.push((pos, fmer, range, range_len));

            // eprintln!("{}: {} / {} -> {} .. Is own rc? {}", 
            //     pos, cmer.to_string().unwrap(), cmer.rc().to_string().unwrap(), range_len, cmer.is_own_rc());
        }


        let (duration, _) = time(|| self.sort_ranges(stats));
        stats.ranges += self.ranges.len();
        stats.time_range_sorting += duration;
        



        

        // let mut last_taken_pos = self.ranges.first().unwrap().0;
        // let mut iterations = 0;
        // let mut penalty = 0;
        // self.range_iteration_order.push((0, 0));

        // self.range_iteration_order.clear();
        // for (i, (rpos, _, _, range_len)) in self.ranges[1..].iter().enumerate() {
        //     let penalty = K - (rpos - last_taken_pos);
        //     let absdiff = rpos.abs_diff(last_taken_pos);
        //     if absdiff > K {
        //         println!("-----> Take {} absdiff {} last pos {} current pos {}", i, absdiff, last_taken_pos, rpos);
        //         last_taken_pos = *rpos;
        //         iterations = 0;
        //     }
        //     self.range_iteration_order.push((penalty, i));
        //     iterations += 1;

        //     println!("index {} rpos {} penalty {} .... {}", i, rpos, penalty, last_taken_pos);
        // }

        // for pair in &self.range_iteration_order {
        //     dbg!(pair);
        // }

        // let mut s= String::new();
        // stdin().read_line(&mut s).expect("Did not enter a correct string");

    }

    pub fn sort_ranges(&mut self, _stats: &mut Stats) {
        // Faster than glidesort in this example
        self.ranges.sort_unstable_by_key(|r| r.2.positions.len());
    }

    pub fn get_anchors(&mut self, _stats: &mut Stats) {
        if self.seeds.len() < 10 { return }

        let mut stop = false;
        for seed in self.seeds.iter() {
            if seed.mismatch > 0 { stop = true };
        }

        if stop {
            for seed in self.seeds.iter() {
                println!("{}", seed);
            }
            let mut s= String::new();
            stdin().read_line(&mut s).expect("Did not enter a correct string");
        }

    }
}

