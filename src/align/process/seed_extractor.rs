

use crate::{align::{common::SeedExtractor, data_structures::common::Seed, stats::Stats}, flexalign::stage_time};

use flexmap::range::AnyRange;

use super::range_extractor::Range;

/// `AnyRange::flank_distances`: each flank cell's distance to `flex` into `out`, and the minimum
/// (`u32::MAX` for none). A packed range's cells go 8 to an AVX2 step.
#[inline]
pub fn flank_distances(range: &AnyRange, flex: u32, out: &mut Vec<u8>) -> u32 {
    match range {
        AnyRange::Packed(p) => flank_distances_bytes(p.flex_bytes(), flex, out),
        _ => range.flank_distances(flex, out),
    }
}

/// [`flank_distances`] over packed flank cells (4 bytes each, little-endian).
#[inline]
fn flank_distances_bytes(cells: &[u8], flex: u32, out: &mut Vec<u8>) -> u32 {
    let n = cells.len() / 4;
    out.clear();
    out.reserve(n);
    let mut min = u32::MAX;
    let mut i = 0usize;
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    unsafe {
        use core::arch::x86_64::*;
        let want = _mm256_set1_epi32(flex as i32);
        let m55 = _mm256_set1_epi32(0x5555_5555);
        let m33 = _mm256_set1_epi32(0x3333_3333);
        let m0f = _mm256_set1_epi32(0x0F0F_0F0F);
        let m3f = _mm256_set1_epi32(0x3F);
        // The low byte of each 32-bit lane to the front of its 128-bit half, then both halves'
        // fronts together.
        let low_bytes = _mm256_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        let fronts = _mm256_setr_epi32(0, 4, 0, 0, 0, 0, 0, 0);
        let mut vmin = _mm256_set1_epi32(-1);
        let dst = out.as_mut_ptr();
        while i + 8 <= n {
            let x = _mm256_xor_si256(_mm256_loadu_si256(cells.as_ptr().add(4 * i) as *const __m256i), want);
            // A base differs where either of its two bits does (flank_distance).
            let v = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 1)), m55);
            let v = _mm256_add_epi32(_mm256_and_si256(v, m33), _mm256_and_si256(_mm256_srli_epi32(v, 2), m33));
            let v = _mm256_and_si256(_mm256_add_epi32(v, _mm256_srli_epi32(v, 4)), m0f);
            let v = _mm256_add_epi32(v, _mm256_srli_epi32(v, 8));
            let d = _mm256_and_si256(_mm256_add_epi32(v, _mm256_srli_epi32(v, 16)), m3f);
            vmin = _mm256_min_epu32(vmin, d);
            let packed = _mm256_permutevar8x32_epi32(_mm256_shuffle_epi8(d, low_bytes), fronts);
            _mm_storel_epi64(dst.add(i) as *mut __m128i, _mm256_castsi256_si128(packed));
            i += 8;
        }
        let mut lanes = [0u32; 8];
        _mm256_storeu_si256(lanes.as_mut_ptr() as *mut __m256i, vmin);
        min = lanes.iter().copied().fold(min, u32::min);
        out.set_len(i);
    }
    for c in cells[4 * i..].chunks_exact(4) {
        let d = flexmap::packed::flank_distance(u32::from_le_bytes([c[0], c[1], c[2], c[3]]), flex);
        out.push(d as u8);
        min = min.min(d);
    }
    min
}

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
        if retrieved_matches >= max_ranges {
            continue;
        }
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
        max_best_flex: usize,
        flank_slack: u32,
        dist_buf: &mut Vec<u8>,
        out: &mut [u32],
    ) -> Option<(u32, usize)> {
        debug_assert!(max_best_flex <= out.len(), "out[] must hold a full kept set");
        // dist is a hamming distance over F<=16 flank bases, fits u8
        let min_dist = flank_distances(&self.vrange, self.fmer.0 as u32, dist_buf);
        if dist_buf.is_empty() {
            return Some((0, 0));
        }
        let threshold = min_dist + flank_slack;
        let mut count = 0usize;
        let mut stored = 0usize;
        // Branch-free on the distance (which is data, and mispredicted): every index is written to
        // the next free slot and kept only if within the threshold.
        for (i, &d) in dist_buf.iter().enumerate() {
            let keep = (d as u32 <= threshold) as usize;
            count += keep;
            if stored < out.len() {
                out[stored] = i as u32;
                stored += keep;
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
        self.emit_seeds_inner::<K, C>(max_best_flex, flank_slack, dist_buf, seeds)
    }
}

impl<'a, const F: usize> Range<'a, F> {
    fn emit_seeds_inner<const K: usize, const C: usize>(
        &self,
        max_best_flex: usize,
        flank_slack: u32,
        dist_buf: &mut Vec<u8>,
        seeds: &mut Vec<Seed>,
    ) -> RangeVerdict {
        match self.vrange.has_header() {
            true => {
                let mut small = [0u32; MAX_KEPT_TIES];
                // `--adaptive-ties` keeps more than the stack buffer holds; rare, so it allocates.
                let mut big = Vec::new();
                let idx: &mut [u32] = if max_best_flex <= MAX_KEPT_TIES {
                    &mut small
                } else {
                    big.resize(max_best_flex, 0);
                    &mut big
                };
                match self.resolve_flex(max_best_flex, flank_slack, dist_buf, idx) {
                    None => RangeVerdict::Discarded,
                    Some((_min_dist, n)) => {
                        for &i in &idx[..n] {
                            let (value, rpos) = self.vrange.value(i as usize);
                            // Each seed carries its own flank distance, not the range minimum.
                            let dist = dist_buf[i as usize] as u32;
                            seeds.push(Seed::from_flexmer::<K, C, F>(self.qpos, rpos, value, dist));
                        }
                        RangeVerdict::Headered
                    }
                }
            }
            false => {
                dist_buf.clear(); // no flanks scanned; keep dist_buf.len() an honest scan count
                for i in 0..self.vrange.len() {
                    let (value, rpos) = self.vrange.value(i);
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
    /// The last `generate`'s (ranges seeded, ranges dropped for too many flank ties).
    last_ties: (usize, usize),
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
            last_ties: (0, 0),
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
        self.last_ties = (retrieved_ranges, discarded_max_flex_count);

        if retrieved_ranges < self.min_ranges && discarded_max_flex_count > 0  {
            // eprintln!("----------------- Recover Ranges....");
            let _old_ranges = ranges;
            let (_ranges, _discarded_max_flex_count) = self.retrieve_seeds(
                ranges,
                128.min(cap),
                self.max_range_size,
                stats
            );
            self.last_ties = (_ranges, _discarded_max_flex_count);
            // eprintln!("{} -> {} (Still discarded: {})", old_ranges, ranges, discarded_max_flex_count);
        }
        
        // stats.time_range_header += duration;
        stats.seeds += self.seeds.len();


        let (duration, _) = stage_time(|| {
            // std's stable sort: a read's seeds are few, and it sorts them in place without the
            // scratch buffer glidesort allocates per call. Same order (both stable).
            crate::align::common::sort_short_by_key(&mut self.seeds, |seed: &Seed| {
                (seed.rval, seed.rpos)//, seed.length) // seed.offset
            });
        });
        stats.time_seed_sorting += duration;

        &self.seeds
    }

    fn retrieve(&self) -> &[Seed] {
        &self.seeds
    }

    fn take_seeds(&mut self) -> Vec<Seed> {
        std::mem::take(&mut self.seeds)
    }

    fn put_seeds(&mut self, seeds: Vec<Seed>) {
        self.seeds = seeds;
    }

    fn max_best_flex(&self) -> usize {
        self.max_best_flex
    }

    fn set_max_best_flex(&mut self, v: usize) {
        self.max_best_flex = v;
    }

    fn last_ties(&self) -> (usize, usize) {
        self.last_ties
    }
}

#[cfg(test)]
mod flank_distance_tests {
    use super::flank_distances_bytes;

    #[test]
    fn simd_distances_match_scalar() {
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x as u32 };
        let mut out = Vec::new();
        for n in 0..40usize {
            let flex = next();
            // Cells near the flank too, so that small distances and the zero distance occur.
            let cells: Vec<u32> = (0..n).map(|i| if i % 3 == 0 { flex ^ (next() & next() & next()) } else { next() }).collect();
            let bytes: Vec<u8> = cells.iter().flat_map(|c| c.to_le_bytes()).collect();
            let min = flank_distances_bytes(&bytes, flex, &mut out);
            let want: Vec<u8> = cells.iter().map(|&c| flexmap::packed::flank_distance(c, flex) as u8).collect();
            assert_eq!(out, want, "n {n}");
            assert_eq!(min, want.iter().map(|&d| d as u32).min().unwrap_or(u32::MAX), "n {n}");
        }
    }
}
