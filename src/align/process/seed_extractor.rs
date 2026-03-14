
use flexmap::VD;

use crate::{align::{common::SeedExtractor, data_structures::common::Seed, stats::Stats}, flexalign::time};
use itertools::Itertools;

use super::range_extractor::Range;

#[derive(Clone)]
pub struct StdSeedExtractor<const K: usize, const C: usize, const F: usize> {
    pub seeds: Vec<Seed>,
    pub max_best_flex: usize,
    pub max_range_size: usize,
    pub min_ranges: usize,
    pub max_ranges: usize,
}

impl<const K: usize, const C: usize, const F: usize> StdSeedExtractor<K, C, F> {
    pub fn new(max_best_flex: usize, max_range_size: usize, min_ranges: usize, max_ranges: usize) -> Self {
        Self {
            seeds: Vec::new(),
            max_best_flex,
            max_range_size,
            min_ranges,
            max_ranges,
        }
    }

    pub fn retrieve_seeds(
        &mut self,
        ranges: &[Range<F>],
        max_best_flex: usize, 
        max_ranges: usize, 
        stats: &mut Stats) -> (usize, usize) {

        let mut retrieved_matches = 0;
        let mut discarded_matches = 0;

        let mut discarded_max_flex_count = 0;

        // eprintln!("{}", ranges.iter().map(|(_,_,_,size)| size.to_string())
        //     .collect_vec()
        //     .join(","));

        // Iterate VRanges, each vrange is one 
        // for (qpos, flex, range, _range_size) in ranges {
        for range in ranges {
            match range.vrange.header {
                Some(headers) => {
                    let mut min_dist = u32::MAX;
                    let mut count = 0;
                    for header in headers {
                        let dist = header.dist(range.fmer.0 as u32);
                        if dist < min_dist { min_dist = dist; count = 0; }
                        if dist == min_dist { count += 1};
                    }
                    
                    let take = count <= max_best_flex;
                    // eprintln!("{} Range count = {}/{} < {}", if take { "X".green() } else { "O".red() }, count, range.positions.len(), self.options.args.max_best_flex);
                    if !take {
                        discarded_max_flex_count += 1;
                        continue;
                    }
                    
                    // eprintln!("Header------");
                    for (index, header) in headers.iter().enumerate() {
                        let dist = header.dist(range.fmer.0 as u32);
                        if dist == min_dist {    
                            let (value, rpos) = VD::get(range.vrange.positions[index].0);
                            self.seeds.push(Seed::from_flexmer::<K,C,F>(range.qpos, rpos, value, dist));
                        }
                    }
                    retrieved_matches += 1;
                },
                None => {
                    for cell in range.vrange.positions {
                        // self.seeds.push((*pos, cell.clone()));
                        let (value, rpos) = VD::get(cell.0);
                        self.seeds.push(Seed::from_coremer::<K,C,F>(range.qpos, rpos, value));
                    }
                },
            };
            // if retrieved_matches >= max_ranges { break }
            // eprintln!("RetrievedMatches: {} >= {}", retrieved_matches, self.max_ranges);
            if retrieved_matches >= self.max_ranges { 
                // eprintln!("{}/{}/{}", retrieved_matches, self.max_ranges, ranges.len()); 
                break 
            }
        }
        stats.retrieved_ranges += retrieved_matches;

        (retrieved_matches, discarded_max_flex_count)
    }

}


impl<const K: usize, const C: usize, const F: usize> SeedExtractor<F> for StdSeedExtractor<K, C, F> {
    fn generate(&mut self, ranges: &[Range<F>], stats: &mut crate::align::stats::Stats) -> &[Seed] {
        self.seeds.clear();

        let (retrieved_ranges, discarded_max_flex_count) = self.retrieve_seeds(
            ranges, 
            self.max_best_flex,
            self.max_range_size,
            stats
        );
        
        if retrieved_ranges < self.min_ranges && discarded_max_flex_count > 0  {
            // eprintln!("----------------- Recover Ranges....");
            let old_ranges = ranges;
            let (ranges, discarded_max_flex_count) = self.retrieve_seeds(
                ranges,
                128,
                self.max_range_size,
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