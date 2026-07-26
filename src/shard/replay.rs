//! The rejoin side of the cut: `RangeRecord` -> `[Seed]` (ProjectShard.md §9).
//!
//! This is `StdSeedExtractor` with `VRange` swapped for `RangeRecord` and the flex filter already
//! applied upstream. It hands `StdPairedAnchorExtractor::generate` exactly what the unsharded path
//! hands it -- `&[Seed]` per mate, sorted by `(rval, rpos)` -- so everything below the cut runs
//! unmodified.
//!
//! The walk itself is not reimplemented: `walk_ranges` is the same code the unsharded path runs
//! (hurdle 11). Only the range *source* differs.

use flexmap::VD;

use crate::align::data_structures::common::Seed;
use crate::align::process::seed_extractor::{walk_ranges, RangeVerdict, SeedRange};
use crate::flexalign::time;

use super::record::{GroupRecord, RangeRecord};

impl<const K: usize, const C: usize, const F: usize> SeedRange<K, C, F> for RangeRecord {
    /// `max_best_flex` is ignored: the shard pass already applied the flex filter and dropped what
    /// failed, so a `RangeRecord` on the wire has by definition survived it and can never be
    /// `Discarded` here. Re-testing would be wrong as well as redundant -- `count` (the number of
    /// headers tying at `min_dist`) is not recoverable from the record, and `cells.len()` is not a
    /// substitute for it.
    ///
    /// The rejoin does not get to re-decide this. `Manifest::verify_compatible` is what stops a
    /// run from replaying flex=16 evidence while believing it was written at flex=32.
    fn emit_seeds(
        &self,
        _max_best_flex: usize,
        _flank_slack: u32,
        _dist_buf: &mut Vec<u8>,
        seeds: &mut Vec<Seed>,
    ) -> RangeVerdict {
        if self.headered {
            for cell in &self.cells {
                let (value, rpos) = VD::get(*cell);
                seeds.push(Seed::from_flexmer::<K, C, F>(
                    self.qpos as usize,
                    rpos,
                    value,
                    self.min_dist as u32,
                ));
            }
            RangeVerdict::Headered
        } else {
            for cell in &self.cells {
                let (value, rpos) = VD::get(*cell);
                seeds.push(Seed::from_coremer::<K, C, F>(self.qpos as usize, rpos, value));
            }
            RangeVerdict::Headerless
        }
    }
}

/// Rebuilds one mate's seeds from the evidence every shard wrote for it.
///
/// Shaped like `StdSeedExtractor` but not implementing `SeedExtractor`: that trait takes
/// `&[Range<F>]`, a live borrow of the map, which is precisely what the rejoin does not have.
#[derive(Clone)]
pub struct ReplaySeedExtractor<const K: usize, const C: usize, const F: usize> {
    seeds: Vec<Seed>,
    ranges: Vec<RangeRecord>,
    max_ranges: usize,
    /// Unused scratch to satisfy `walk_ranges` (the replay carries pre-resolved cells and never
    /// scans flanks); kept as a field so the call allocates nothing.
    flex_scratch: Vec<u8>,
}

impl<const K: usize, const C: usize, const F: usize> ReplaySeedExtractor<K, C, F> {
    pub fn new(max_ranges: usize) -> Self {
        Self { seeds: Vec::new(), ranges: Vec::new(), max_ranges, flex_scratch: Vec::new() }
    }

    /// Collects one mate's ranges out of the groups the N shard cursors yielded for this read.
    ///
    /// Groups arrive per shard, each holding the ranges whose core-mer that shard owned. Ordering
    /// across shards is arbitrary and meaningless -- `generate` re-establishes the global order.
    pub fn gather(&mut self, groups: &[GroupRecord], mate: u8) {
        self.ranges.clear();
        for g in groups.iter().filter(|g| g.mate == mate) {
            self.ranges.extend(g.ranges.iter().cloned());
        }
    }

    /// Restores the global rank order, walks the budget, and sorts the seeds.
    ///
    /// The sort key is `(n_positions, qpos)`, not `n_positions` alone. The unsharded path uses
    /// `sort_unstable_by_key(|r| positions.len())`, whose tie-break falls out of the pre-sort
    /// (k-mer) order -- reproducing that here would mean re-sorting by qpos first just to
    /// reconstruct an order nothing documents. Since sharded output is not required to be
    /// bit-identical (§4), the explicit tie-break is the better trade: one sort, deterministic
    /// across shard counts, which the unsharded path arguably is not across refactors (§3.3).
    pub fn generate(&mut self, time_sorting: &mut std::time::Duration) -> &[Seed] {
        self.seeds.clear();

        self.ranges.sort_unstable_by_key(|r| (r.n_positions, r.qpos));

        // usize::MAX: the flex filter lives in the shard pass. See `emit_seeds` above. flank_slack
        // is irrelevant here for the same reason (cells are pre-resolved on the wire).
        let _: (usize, usize, usize) = walk_ranges::<K, C, F, _>(
            &self.ranges,
            usize::MAX,
            0,
            self.max_ranges,
            &mut self.flex_scratch,
            &mut self.seeds,
        );

        let (duration, _) = time(|| {
            glidesort::sort_by_key(&mut self.seeds, |seed: &Seed| (seed.rval, seed.rpos));
        });
        *time_sorting += duration;

        &self.seeds
    }

    pub fn retrieve(&self) -> &[Seed] {
        &self.seeds
    }

    pub fn ranges(&self) -> &[RangeRecord] {
        &self.ranges
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(qpos: u32, n_positions: u32, headered: bool, min_dist: u8, cells: &[u64]) -> RangeRecord {
        RangeRecord { qpos, n_positions, headered, min_dist, cells: cells.to_vec() }
    }

    fn group(mate: u8, ranges: Vec<RangeRecord>) -> GroupRecord {
        GroupRecord { idx: 0, mate, ranges }
    }

    #[test]
    fn gather_takes_only_the_requested_mate() {
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        let groups = vec![
            group(0, vec![rec(1, 5, true, 0, &[1])]),
            group(1, vec![rec(2, 5, true, 0, &[2])]),
            group(0, vec![rec(3, 5, true, 0, &[3])]),
        ];
        r.gather(&groups, 0);
        assert_eq!(r.ranges().len(), 2);
        r.gather(&groups, 1);
        assert_eq!(r.ranges().len(), 1);
        assert_eq!(r.ranges()[0].qpos, 2);
    }

    #[test]
    fn gather_merges_every_shard_and_restores_rank_order() {
        // Shards see disjoint core-mers, so their groups arrive in no particular order. The
        // rejoin's job is to put the read's ranges back into one global rarest-first order.
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        let groups = vec![
            group(0, vec![rec(10, 90, true, 0, &[1])]),
            group(0, vec![rec(20, 3, true, 0, &[2])]),
            group(0, vec![rec(30, 40, true, 0, &[3])]),
        ];
        r.gather(&groups, 0);
        let mut d = std::time::Duration::ZERO;
        r.generate(&mut d);
        let order: Vec<u32> = r.ranges().iter().map(|x| x.n_positions).collect();
        assert_eq!(order, vec![3, 40, 90], "ranges must end up rarest-first");
    }

    #[test]
    fn ties_break_on_qpos_deterministically() {
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        let groups = vec![group(
            0,
            vec![rec(50, 7, true, 0, &[1]), rec(10, 7, true, 0, &[2]), rec(30, 7, true, 0, &[3])],
        )];
        r.gather(&groups, 0);
        let mut d = std::time::Duration::ZERO;
        r.generate(&mut d);
        let order: Vec<u32> = r.ranges().iter().map(|x| x.qpos).collect();
        assert_eq!(order, vec![10, 30, 50], "equal n_positions must order by qpos");
    }

    #[test]
    fn headerless_ranges_do_not_consume_the_budget() {
        // §3.2: only headered ranges count against --ranges. With max_ranges=1, the headered range
        // spends the budget, but headerless ranges sorted before it still emit.
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(1);
        let groups = vec![group(
            0,
            vec![
                rec(1, 1, false, 0, &[0x11]),
                rec(2, 2, false, 0, &[0x22]),
                rec(3, 9, true, 0, &[0x33]),
                rec(4, 10, true, 0, &[0x44]),
            ],
        )];
        r.gather(&groups, 0);
        let mut d = std::time::Duration::ZERO;
        let seeds = r.generate(&mut d);
        // Two headerless + the first headered; the second headered is past the cutoff.
        assert_eq!(seeds.len(), 3);
    }

    #[test]
    fn the_budget_cuts_off_headered_ranges() {
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(2);
        let ranges: Vec<RangeRecord> =
            (0..10).map(|i| rec(i, i + 1, true, 0, &[i as u64 + 1])).collect();
        r.gather(&[group(0, ranges)], 0);
        let mut d = std::time::Duration::ZERO;
        let seeds = r.generate(&mut d);
        assert_eq!(seeds.len(), 2, "only max_ranges headered ranges should be retrieved");
    }

    #[test]
    fn seeds_come_out_sorted_by_rval_rpos() {
        // The contract with StdPairedAnchorExtractor (§3).
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        let cells: Vec<u64> = (0..32u64).map(|i| (i * 2_654_435_761) & ((1 << 60) - 1)).collect();
        r.gather(&[group(0, vec![rec(1, 40, false, 0, &cells)])], 0);
        let mut d = std::time::Duration::ZERO;
        let seeds = r.generate(&mut d);
        assert_eq!(seeds.len(), 32);
        assert!(
            seeds.windows(2).all(|w| (w[0].rval, w[0].rpos) <= (w[1].rval, w[1].rpos)),
            "seeds must be sorted by (rval, rpos)"
        );
    }

    #[test]
    fn generate_clears_between_calls() {
        // §3.4: the unsharded recovery path appends to a non-cleared buffer and re-emits every
        // seed the first walk produced. The replay must not inherit that -- it has no recovery
        // path at all (§4), and a stale buffer would duplicate seeds across reads.
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        let groups = vec![group(0, vec![rec(1, 5, true, 0, &[7, 8])])];
        let mut d = std::time::Duration::ZERO;

        r.gather(&groups, 0);
        let first = r.generate(&mut d).len();
        r.gather(&groups, 0);
        let second = r.generate(&mut d).len();

        assert_eq!(first, second, "a second generate must not accumulate seeds");
    }

    #[test]
    fn a_read_absent_from_every_shard_yields_nothing() {
        let mut r = ReplaySeedExtractor::<31, 15, 16>::new(15);
        r.gather(&[], 0);
        let mut d = std::time::Duration::ZERO;
        assert!(r.generate(&mut d).is_empty());
    }
}
