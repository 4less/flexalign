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
//! The min-distance resolution is *not* reimplemented here -- it comes from `Range::min_flank_dist`
//! and `Range::cells_at_dist`, the same helpers the unsharded walk uses, so both modes select the
//! same cells by construction (hurdle 11).

use crate::align::process::range_extractor::Range;
use crate::align::process::seed_extractor::RangeVerdict;

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

/// Turns one mate's ranges into the records this shard owes the rejoin.
///
/// `ranges` must be ordered rarest-first (`StdRangeExtractor` sorts by `positions.len()`
/// ascending). The cap relies on that order; the rejoin re-establishes it globally after merging.
///
/// Appends to `out` -- a read with no surviving range in this shard appends nothing, which is the
/// common case and is why a group is only written when it has content.
pub fn emit_ranges<const F: usize>(
    ranges: &[Range<F>],
    max_best_flex: usize,
    max_headered: usize,
    out: &mut Vec<RangeRecord>,
) -> EmitStats {
    let mut stats = EmitStats::default();

    for range in ranges {
        match range.vrange.header {
            Some(headers) => {
                debug_assert!(
                    !headers.is_empty(),
                    "a headered block has size > HEADER_THRESHOLD, so it always has a header"
                );

                let (min_dist, count) = range.min_flank_dist(headers);
                if count > max_best_flex {
                    stats.discarded += 1;
                    continue;
                }
                if stats.headered >= max_headered {
                    stats.capped += 1;
                    continue;
                }

                debug_assert!(min_dist <= F as u32, "dist {} exceeds {} flank bases", min_dist, F);
                let cells: Vec<u64> = range.cells_at_dist(headers, min_dist).collect();

                stats.headered += 1;
                stats.cells += cells.len();
                out.push(RangeRecord {
                    qpos: range.qpos as u32,
                    n_positions: range.vrange.positions.len() as u32,
                    headered: true,
                    min_dist: min_dist as u8,
                    cells,
                });
            }
            None => {
                let cells: Vec<u64> = range.vrange.positions.iter().map(|c| c.0).collect();
                stats.headerless += 1;
                stats.cells += cells.len();
                out.push(RangeRecord {
                    qpos: range.qpos as u32,
                    n_positions: range.vrange.positions.len() as u32,
                    headered: false,
                    // Meaningless without headers; the replay reads it only when `headered`.
                    min_dist: 0,
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
