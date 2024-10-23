
use flexmap::VD;

use crate::{align::{common::SeedExtractor, data_structures::Seed, stats::Stats}, flexalign::time};

use super::range_extractor::Range;

#[derive(Clone)]
pub struct StdSeedExtractor<const K: usize, const C: usize, const F: usize> {
    pub seeds: Vec<Seed>,
    pub max_best_flex: usize,
    pub max_ranges: usize,
    pub min_ranges: usize,
}

impl<const K: usize, const C: usize, const F: usize> StdSeedExtractor<K, C, F> {
    pub fn new(max_best_flex: usize, max_ranges: usize, min_ranges: usize) -> Self {
        Self {
            seeds: Vec::new(),
            max_best_flex,
            max_ranges,
            min_ranges,
        }
    }

    pub fn retrieve_seeds(
        &mut self,
        ranges: &[Range<F>],
        max_best_flex: usize, 
        max_ranges: usize, 
        stats: &mut Stats) -> (usize, usize) {

        let mut matches = 0;
        let mut discarded_max_flex_count = 0;
        for (qpos, flex, range, _range_size) in ranges {
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
                            let (value, rpos) = VD::get(range.positions[index].0);
                            self.seeds.push(Seed::from_flexmer::<K,C,F>(*qpos, rpos, value, dist));
                        }
                    }
                },
                None => {
                    for cell in range.positions {
                        // self.seeds.push((*pos, cell.clone()));
                        let (value, rpos) = VD::get(cell.0);
                        self.seeds.push(Seed::from_coremer::<K,C,F>(*qpos, rpos, value));
                    }
                },
            };
            matches += 1;
            if matches >= max_ranges { break }
        }
        (matches, discarded_max_flex_count)
    }

}


impl<const K: usize, const C: usize, const F: usize> SeedExtractor<F> for StdSeedExtractor<K, C, F> {
    fn generate(&mut self, ranges: &[Range<F>], stats: &mut crate::align::stats::Stats) -> &[Seed] {
        self.seeds.clear();

        let (retrieved_ranges, discarded_max_flex_count) = self.retrieve_seeds(
            ranges, 
            self.max_best_flex,
            self.max_ranges,
            stats
        );
        
        if retrieved_ranges < self.min_ranges && discarded_max_flex_count > 0  {
            // eprintln!("----------------- Recover Ranges....");
            let old_ranges = ranges;
            let (ranges, discarded_max_flex_count) = self.retrieve_seeds(
                ranges,
                128,
                self.max_ranges,
                stats
            );
            // eprintln!("{} -> {} (Still discarded: {})", old_ranges, ranges, discarded_max_flex_count);
        }
        
        // stats.time_range_header += duration;
        stats.seeds += self.seeds.len();


        let (duration, _) = time(|| {
            glidesort::sort_by_key(&mut self.seeds, |seed: &Seed| {
                (seed.rval, seed.rpos)//, seed.length) // seed.offset
            });
        });
        stats.time_seed_sorting += duration;

        &self.seeds
    }

    fn retrieve(&self) -> &[Seed] {
        &self.seeds
    }
}