
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
/// Upper bound on the flank ties a *surviving* range can carry: the recovery pass relaxes
/// `max_best_flex` to 128, and nothing above the active `max_best_flex` is ever kept, so a
/// stack buffer of this size holds every kept range's survivor indices without allocating.
pub const MAX_KEPT_TIES: usize = 128;

pub trait SeedRange<const K: usize, const C: usize, const F: usize> {
    /// Pushes this range's seeds onto `seeds` and reports what it did.
    ///
    /// `flank_slack` widens the kept flank set to every cell within that many mismatches of the
    /// best (minimum) distance (0 = best-distance ties only). `dist_buf` is a reused scratch the
    /// header-scanning implementation caches per-flank distances in; implementations that carry
    /// pre-resolved cells ignore both.
    fn emit_seeds(
        &self,
        max_best_flex: usize,
        flank_slack: u32,
        dist_buf: &mut Vec<u8>,
        seeds: &mut Vec<Seed>,
    ) -> RangeVerdict;
}

/// Walks ranges in order, emitting seeds until the `--ranges` budget is spent.
///
/// `ranges` must already be ordered rarest-first: `StdRangeExtractor` sorts by
/// `positions.len()` ascending, and the rejoin re-establishes that order by `(n_positions, qpos)`
/// after merging shards (§3.3). The cutoff is a **global rank**, so this order is load-bearing --
/// it is the one real cross-key dependency in the design (§5).
///
/// Returns `(retrieved, discarded, headers_scanned)`.
pub fn walk_ranges<const K: usize, const C: usize, const F: usize, R>(
    ranges: &[R],
    max_best_flex: usize,
    flank_slack: u32,
    max_ranges: usize,
    dist_buf: &mut Vec<u8>,
    seeds: &mut Vec<Seed>,
) -> (usize, usize, usize)
where
    R: SeedRange<K, C, F>,
{
    let mut retrieved_matches = 0;
    let mut discarded_max_flex_count = 0;
    let mut headers_scanned = 0usize;

    for range in ranges {
        let verdict = range.emit_seeds(max_best_flex, flank_slack, dist_buf, seeds);
        // `resolve_flex` leaves one cached distance per scanned flank in `dist_buf`, so its length
        // is exactly this range's header-scan count (0 for headerless ranges and for the replay,
        // which pre-resolves cells and never touches `dist_buf`).
        headers_scanned += dist_buf.len();
        match verdict {
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

    (retrieved_matches, discarded_max_flex_count, headers_scanned)
}

impl<'a, const F: usize> Range<'a, F> {
    /// One-pass flex resolution: computes each flank distance **exactly once** (cached into
    /// `dist_buf`), finds the minimum, then selects the header indices whose distance is within
    /// `flank_slack` of it into `out[..n]`.
    ///
    /// Returns `Some((min_dist, n))`, or `None` if the surviving set exceeds `max_best_flex` (the
    /// range is discarded -- and, crucially, its `positions` are never read, so a repetitive range
    /// costs no access to the cold blob pages the positions live in). After the call `dist_buf[i]`
    /// holds header `i`'s distance, so the caller tags each surviving seed with its *own* distance
    /// (which matters: `Seed::from_flexmer` shifts and lengths the seed by it).
    ///
    /// Replaces the old two-pass `min_flank_dist` + `cells_at_dist` (which walked the flanks twice
    /// and recomputed every distance). Shared by the unsharded walk and the shard emitter, so both
    /// still select identical cells by construction (hurdle 11). `out.len()` must be
    /// >= `max_best_flex`; [`MAX_KEPT_TIES`] is the standard size.
    pub fn resolve_flex(
        &self,
        headers: &[HeaderSeq],
        max_best_flex: usize,
        flank_slack: u32,
        dist_buf: &mut Vec<u8>,
        out: &mut [u32],
    ) -> Option<(u32, usize)> {
        debug_assert!(max_best_flex <= out.len(), "out[] must hold a full kept set");
        dist_buf.clear();
        if headers.is_empty() {
            return Some((0, 0));
        }
        dist_buf.reserve(headers.len());
        let flex = self.fmer.0 as u32;
        let mut min_dist = u32::MAX;
        for header in headers {
            let d = header.dist(flex);
            dist_buf.push(d as u8); // dist is a hamming distance over F<=16 flank bases, fits u8
            if d < min_dist {
                min_dist = d;
            }
        }
        let threshold = min_dist + flank_slack;
        let mut count = 0usize;
        let mut stored = 0usize;
        for (i, &d) in dist_buf.iter().enumerate() {
            if d as u32 <= threshold {
                count += 1;
                if stored < out.len() {
                    out[stored] = i as u32;
                    stored += 1;
                }
            }
        }
        if count > max_best_flex {
            None
        } else {
            debug_assert_eq!(stored, count, "a kept set must fit entirely in out[]");
            Some((min_dist, count))
        }
    }
}

impl<'a, const K: usize, const C: usize, const F: usize> SeedRange<K, C, F> for Range<'a, F> {
    fn emit_seeds(
        &self,
        max_best_flex: usize,
        flank_slack: u32,
        dist_buf: &mut Vec<u8>,
        seeds: &mut Vec<Seed>,
    ) -> RangeVerdict {
        match self.vrange.header {
            Some(headers) => {
                let mut idx = [0u32; MAX_KEPT_TIES];
                match self.resolve_flex(headers, max_best_flex, flank_slack, dist_buf, &mut idx) {
                    None => RangeVerdict::Discarded,
                    Some((_min_dist, n)) => {
                        for &i in &idx[..n] {
                            let (value, rpos) = VD::get(self.vrange.positions[i as usize].0);
                            // Each seed carries its own flank distance, not the range minimum.
                            let dist = dist_buf[i as usize] as u32;
                            seeds.push(Seed::from_flexmer::<K, C, F>(self.qpos, rpos, value, dist));
                        }
                        RangeVerdict::Headered
                    }
                }
            }
            None => {
                dist_buf.clear(); // no flanks scanned; keep dist_buf.len() an honest scan count
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
    /// `--flank-slack`: keep flank cells within this many mismatches of the best distance.
    pub flank_slack: u32,
    /// Reused per-flank distance scratch for the single-pass flex resolution (avoids reallocating
    /// per range/read).
    flex_scratch: Vec<u8>,
}

impl<const K: usize, const C: usize, const F: usize> StdSeedExtractor<K, C, F> {
    pub fn new(max_best_flex: usize, max_range_size: usize, min_ranges: usize, max_ranges: usize, mask_flank_mult: usize, flank_slack: u32) -> Self {
        Self {
            seeds: Vec::new(),
            max_best_flex,
            max_range_size,
            min_ranges,
            max_ranges,
            mask_flank_mult,
            flank_slack,
            flex_scratch: Vec::new(),
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

        let (retrieved_matches, discarded_max_flex_count, headers_scanned) = walk_ranges::<K, C, F, _>(
            ranges,
            max_best_flex,
            self.flank_slack,
            self.max_ranges,
            &mut self.flex_scratch,
            &mut self.seeds,
        );
        stats.retrieved_ranges += retrieved_matches;
        stats.header_elements += headers_scanned;

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
