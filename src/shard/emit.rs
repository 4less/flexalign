//! The shard side of the cut: `Range` -> `RangeRecord` (ProjectShard.md §4, §7).
//!
//! A shard pass runs the real `StdKmerExtractor` and `StdRangeExtractor` against a `ShardedDB`,
//! then stops. Instead of building seeds it writes down what it found, and the rejoin builds the
//! seeds later.
//!
//! Two things make that cheap. First, `Seed::from_flexmer::<K, C, F>(qpos, rpos, value, dist)` is a
//! pure function of the cell plus `qpos` and `dist`, and both of those are per-*range*, not
//! per-seed -- so all the per-seed metadata hoists into the range header and the payload is the raw
//! `VCell` verbatim. Note the signature below: the emitter is generic over `F` alone. It needs no
//! `K` or `C` because it never constructs a `Seed`.
//!
//! Second, the shard applies the flex filter itself and emits **nothing** for ranges that fail
//! (§4). That is where the design spends its exactness budget: staying exact would force every
//! range with `16 < count <= 128` to carry up to 128 cells (~1 KB against a ~21 B typical range),
//! on every read touching a repetitive k-mer, to be discarded unread except on the rare read that
//! triggers recovery.
//!
//! The min-distance resolution is *not* reimplemented here -- it comes from `Range::resolve_flex`,
//! the same one-pass helper the unsharded walk uses, so both modes select the same cells by
//! construction (hurdle 11). The wire carries one distance per range, so the shard always resolves
//! at `flank_slack = 0` (the widened window is an unsharded-only option).

use crate::align::process::range_extractor::Range;
use crate::align::process::seed_extractor::{RangeVerdict, MAX_KEPT_TIES};

use super::record::RangeRecord;

/// What a shard pass did with one read's ranges. Deliberately not the pipeline's `Stats`: those
/// counters accumulate per-extractor and a sharded run would count them N times over (hurdle 9).
/// The rejoin owns the logical stats; a shard pass reports only its own work.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EmitStats {
    /// Headered ranges written. These consume the rejoin's `--ranges` budget.
    pub headered: usize,
    /// Headerless ranges written. Always emitted: they are the rarest evidence and cost <=2 cells.
    pub headerless: usize,
    /// Ranges dropped by the flex filter. Never reach the wire (§4).
    pub discarded: usize,
    /// Ranges dropped by the per-shard headered cap (§5). Normally zero.
    pub capped: usize,
    pub cells: usize,
}

/// A range in the global `--ranges` rarest set must also be in its own shard's rarest set, so a
/// shard may cap its headered output without affecting what the rejoin would have kept.
///
/// At N=10 a read yields ~20 ranges, ~2 per shard, so this never binds -- it is not a win, it just
/// bounds the pathological case for free. Headerless ranges are exempt: they are the rarest
/// evidence and cost at most two cells.
pub const DEFAULT_MAX_HEADERED: usize = 15;

/// The unsharded seed walk's recovery (`StdSeedExtractor::generate`): a mate with fewer than
/// `min_ranges` headered ranges surviving the flex filter, and at least one range cut by it, is
/// walked again at `flex` -- `128`, under the `--mask-flank-mult` cap -- and those seeds are added.
#[derive(Debug, Clone, Copy)]
pub struct Recovery {
    pub min_ranges: usize,
    pub flex: usize,
}

/// Turns one mate's ranges into the records this shard owes the rejoin.
///
/// With `recovery`, a mate with fewer than `min_ranges` surviving headered ranges *in this shard*
/// also gets a recovery record for every range the flex filter cut: its cells at `recovery.flex`,
/// or none (a marker) if it fails that as well or falls past the per-shard cap there. Only the
/// rejoin can tell whether the mate triggers recovery -- its surviving ranges are spread over the
/// shards -- but a mate below `min_ranges` overall is below it in every shard, so every shard has
/// written what the rejoin then needs.
///
/// `ranges` must be ordered rarest-first (`StdRangeExtractor` sorts by `positions.len()`
/// ascending). The cap relies on that order; the rejoin re-establishes it globally after merging.
///
/// Appends to `out` -- a read with no surviving range in this shard appends nothing, which is the
/// common case and is why a group is only written when it has content.
pub fn emit_ranges<const F: usize>(
    ranges: &[Range<F>],
    max_best_flex: usize,
    recovery: Option<Recovery>,
    max_headered: usize,
    dist_buf: &mut Vec<u8>,
    out: &mut Vec<RangeRecord>,
) -> EmitStats {
    let mut stats = EmitStats::default();
    let mut idx = [0u32; MAX_KEPT_TIES];
    // Which ranges the flex filter cut, for the recovery records below.
    let mut cut: Vec<bool> = Vec::new();

    for (r, range) in ranges.iter().enumerate() {
        match range.vrange.has_header() {
            true => {
                // flank_slack = 0: the sharded record stores one distance per range.
                let (min_dist, count) = match range.resolve_flex(max_best_flex, 0, dist_buf, &mut idx) {
                    None => {
                        stats.discarded += 1;
                        if recovery.is_some() {
                            cut.resize(ranges.len(), false);
                            cut[r] = true;
                        }
                        continue;
                    }
                    Some(mc) => mc,
                };
                if stats.headered >= max_headered {
                    stats.capped += 1;
                    continue;
                }

                debug_assert!(min_dist <= F as u32, "dist {} exceeds {} flank bases", min_dist, F);
                let cells: Vec<u64> =
                    idx[..count].iter().map(|&i| range.vrange.cell(i as usize)).collect();

                stats.headered += 1;
                stats.cells += cells.len();
                out.push(RangeRecord {
                    qpos: range.qpos as u32,
                    n_positions: range.vrange.len() as u32,
                    headered: true,
                    min_dist: min_dist as u8,
                    recovery: false,
                    cells,
                });
            }
            false => {
                let cells: Vec<u64> = (0..range.vrange.len()).map(|i| range.vrange.cell(i)).collect();
                stats.headerless += 1;
                stats.cells += cells.len();
                out.push(RangeRecord {
                    qpos: range.qpos as u32,
                    n_positions: range.vrange.len() as u32,
                    headered: false,
                    // Meaningless without headers; the replay reads it only when `headered`.
                    min_dist: 0,
                    recovery: false,
                    cells,
                });
            }
        }
    }

    // Survivors past the cap are survivors too: a capped mate has max_headered >= min_ranges.
    let survived = stats.headered + stats.capped;
    if let Some(rec) = recovery {
        if stats.discarded > 0 && survived < rec.min_ranges {
            let mut big = vec![0u32; rec.flex.max(MAX_KEPT_TIES)];
            // Headered survivors at the relaxed threshold, in rank order: the per-shard cap holds
            // for the recovery walk as for the first (a range in its global rarest set is in its
            // shard's).
            let mut kept = 0usize;
            for (r, range) in ranges.iter().enumerate() {
                if !range.vrange.has_header() {
                    continue;
                }
                let resolved = range.resolve_flex(rec.flex, 0, dist_buf, &mut big);
                if !cut[r] {
                    kept += resolved.is_some() as usize;
                    continue;
                }
                let cells: Vec<u64> = match resolved {
                    Some((_, count)) if kept < max_headered => {
                        kept += 1;
                        big[..count].iter().map(|&i| range.vrange.cell(i as usize)).collect()
                    }
                    _ => Vec::new(),
                };
                let min_dist = resolved.map_or(0, |(d, _)| d);
                stats.cells += cells.len();
                out.push(RangeRecord {
                    qpos: range.qpos as u32,
                    n_positions: range.vrange.len() as u32,
                    headered: true,
                    min_dist: min_dist as u8,
                    recovery: true,
                    cells,
                });
            }
        }
    }

    stats
}

/// The verdict `emit_ranges` reached for a range, for tests and for reasoning about parity with
/// the unsharded walk.
pub fn verdict_of(record: &RangeRecord) -> RangeVerdict {
    if record.headered {
        RangeVerdict::Headered
    } else {
        RangeVerdict::Headerless
    }
}
