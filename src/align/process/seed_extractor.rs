
use flexmap::values::HeaderSeq;
use flexmap::VD;

use crate::{align::{common::SeedExtractor, data_structures::common::Seed, stats::Stats}, flexalign::time};

use super::range_extractor::Range;

/// What one range did, which is what drives the `--ranges` budget.
///
/// The distinction between `Headerless` and `Headered` is not cosmetic: headerless ranges push
/// seeds but **do not consume the budget** (ProjectShard.md §3.2), which is why the walk counts
/// them separately rather than treating "emitted seeds" as "retrieved".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RangeVerdict {
    /// The flex filter rejected it: too many flank headers tie at the minimum distance.
    Discarded,
    /// A block at or below `HEADER_THRESHOLD` positions, carrying no flank headers.
    Headerless,
    /// A headered block that survived the flex filter. Consumes one unit of the budget.
    Headered,
}

/// One range's evidence, whatever its provenance.
///
/// The unsharded pipeline resolves this from a live `VRange` (scanning flank headers for the
/// min-distance set); the rejoin reads it off the wire, because the shard pass already did that
/// scan and wrote the surviving cells down. The walk below must not be able to tell the
/// difference -- that is the point of the trait.
///
/// ProjectShard.md hurdle 11: two copies of this walk would drift the first time anyone tuned the
/// heuristic, and since sharded output is not bit-compared against unsharded (§4 relaxes
/// exactness), nothing would catch it. So it exists once.
pub trait SeedRange<const K: usize, const C: usize, const F: usize> {
    /// Pushes this range's seeds onto `seeds` and reports what it did.
    fn emit_seeds(&self, max_best_flex: usize, seeds: &mut Vec<Seed>) -> RangeVerdict;
}

/// Walks ranges in order, emitting seeds until the `--ranges` budget is spent.
///
/// `ranges` must already be ordered rarest-first: `StdRangeExtractor` sorts by
/// `positions.len()` ascending, and the rejoin re-establishes that order by `(n_positions, qpos)`
/// after merging shards (§3.3). The cutoff is a **global rank**, so this order is load-bearing --
/// it is the one real cross-key dependency in the design (§5).
///
/// Returns `(retrieved, discarded)`.
pub fn walk_ranges<const K: usize, const C: usize, const F: usize, R>(
    ranges: &[R],
    max_best_flex: usize,
    max_ranges: usize,
    seeds: &mut Vec<Seed>,
) -> (usize, usize)
where
    R: SeedRange<K, C, F>,
{
    let mut retrieved_matches = 0;
    let mut discarded_max_flex_count = 0;

    for range in ranges {
        match range.emit_seeds(max_best_flex, seeds) {
            RangeVerdict::Discarded => {
                discarded_max_flex_count += 1;
                // Note: skips the cutoff check below, matching the original `continue`.
                continue;
            }
            // Emits seeds, spends no budget (§3.2). It can still be cut off by a range that
            // trips the budget before it.
            RangeVerdict::Headerless => {}
            RangeVerdict::Headered => retrieved_matches += 1,
        }
        if retrieved_matches >= max_ranges {
            break;
        }
    }

    (retrieved_matches, discarded_max_flex_count)
}

impl<'a, const F: usize> Range<'a, F> {
    /// The distance to the closest flank header, and how many headers tie at it.
    ///
    /// `count` -- not `positions.len()` -- is what the flex filter tests: a block can hold many
    /// positions while only a few flanks actually match the read.
    pub fn min_flank_dist(&self, headers: &[HeaderSeq]) -> (u32, usize) {
        let mut min_dist = u32::MAX;
        let mut count = 0;
        for header in headers {
            let dist = header.dist(self.fmer.0 as u32);
            if dist < min_dist {
                min_dist = dist;
                count = 0;
            }
            if dist == min_dist {
                count += 1
            }
        }
        (min_dist, count)
    }

    /// The cells whose flank header ties at `min_dist`, paired with that header's index.
    ///
    /// Shared with the shard emitter, which writes exactly these cells to the wire instead of
    /// turning them into seeds -- so the two modes select the same cells by construction.
    pub fn cells_at_dist<'r>(
        &'r self,
        headers: &'r [HeaderSeq],
        min_dist: u32,
    ) -> impl Iterator<Item = u64> + 'r {
        headers
            .iter()
            .enumerate()
            .filter(move |(_, header)| header.dist(self.fmer.0 as u32) == min_dist)
            .map(move |(index, _)| self.vrange.positions[index].0)
    }
}

impl<'a, const K: usize, const C: usize, const F: usize> SeedRange<K, C, F> for Range<'a, F> {
    fn emit_seeds(&self, max_best_flex: usize, seeds: &mut Vec<Seed>) -> RangeVerdict {
        match self.vrange.header {
            Some(headers) => {
                let (min_dist, count) = self.min_flank_dist(headers);
                if count > max_best_flex {
                    return RangeVerdict::Discarded;
                }
                for cell in self.cells_at_dist(headers, min_dist) {
                    let (value, rpos) = VD::get(cell);
                    seeds.push(Seed::from_flexmer::<K, C, F>(self.qpos, rpos, value, min_dist));
                }
                RangeVerdict::Headered
            }
            None => {
                for cell in self.vrange.positions {
                    let (value, rpos) = VD::get(cell.0);
                    seeds.push(Seed::from_coremer::<K, C, F>(self.qpos, rpos, value));
                }
                RangeVerdict::Headerless
            }
        }
    }
}

#[derive(Clone)]
pub struct StdSeedExtractor<const K: usize, const C: usize, const F: usize> {
    pub seeds: Vec<Seed>,
    pub max_best_flex: usize,
    pub max_range_size: usize,
    pub min_ranges: usize,
    pub max_ranges: usize,
    /// Hard ceiling on the min-distance flank tie count (`--mask-flank-mult`); 0 = disabled.
    /// Applied on top of `max_best_flex` in every pass, so — unlike `max_best_flex` — the
    /// recovery pass below cannot relax it.
    pub mask_flank_mult: usize,
}

impl<const K: usize, const C: usize, const F: usize> StdSeedExtractor<K, C, F> {
    pub fn new(max_best_flex: usize, max_range_size: usize, min_ranges: usize, max_ranges: usize, mask_flank_mult: usize) -> Self {
        Self {
            seeds: Vec::new(),
            max_best_flex,
            max_range_size,
            min_ranges,
            max_ranges,
            mask_flank_mult,
        }
    }

    /// The mask as an inclusive upper bound on the tie count: 0 means "no ceiling".
    fn flank_cap(&self) -> usize {
        if self.mask_flank_mult == 0 { usize::MAX } else { self.mask_flank_mult }
    }

    pub fn retrieve_seeds(
        &mut self,
        ranges: &[Range<F>],
        max_best_flex: usize,
        _max_ranges: usize,
        stats: &mut Stats) -> (usize, usize) {

        let (retrieved_matches, discarded_max_flex_count) = walk_ranges::<K, C, F, _>(
            ranges,
            max_best_flex,
            self.max_ranges,
            &mut self.seeds,
        );
        stats.retrieved_ranges += retrieved_matches;

        (retrieved_matches, discarded_max_flex_count)
    }

}


impl<const K: usize, const C: usize, const F: usize> SeedExtractor<F> for StdSeedExtractor<K, C, F> {
    fn generate(&mut self, ranges: &[Range<F>], stats: &mut crate::align::stats::Stats) -> &[Seed] {
        self.seeds.clear();

        // `--mask-flank-mult` is a hard ceiling: cap the flank-tie cutoff in *both* passes so a
        // repetitive flank stays masked even when the recovery pass relaxes `max_best_flex`.
        let cap = self.flank_cap();

        let (retrieved_ranges, discarded_max_flex_count) = self.retrieve_seeds(
            ranges,
            self.max_best_flex.min(cap),
            self.max_range_size,
            stats
        );

        if retrieved_ranges < self.min_ranges && discarded_max_flex_count > 0  {
            // eprintln!("----------------- Recover Ranges....");
            let _old_ranges = ranges;
            let (_ranges, _discarded_max_flex_count) = self.retrieve_seeds(
                ranges,
                128.min(cap),
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
