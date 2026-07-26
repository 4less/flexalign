//! The shard pass driver (ProjectShard.md §13.7).
//!
//! One *pass* runs the pipeline front-half -- `StdKmerExtractor` then `StdRangeExtractor` -- against
//! one shard's key-range view of the index ([`ShardedDB`](super::db::ShardedDB)), and stops. Instead
//! of building seeds it feeds each mate's ranges through [`emit_ranges`] and records what survived,
//! keyed by the read's deterministic [`ReadId::ordinal`]. The rejoin (§13.8) rebuilds the seeds.
//!
//! This is the **in-RAM evidence path** the design says to build first (§10): a pass returns a
//! [`ShardEvidence`] held in memory, so the whole cut can be exercised and validated before the
//! on-disk chunk format and `GroupReader` cursors are wired up. The read id is the flat, file-derived
//! ordinal (see [`ReadId`]); the rejoin gathers a read's evidence by looking that ordinal up in each
//! shard, which is why a pass sorts its groups by `(ordinal, mate)` before returning.
//!
//! Per the design (hurdle 9), a pass reports only its own IO/emit counters, never the pipeline
//! `Stats`: those accumulate per-extractor and would be counted N times over across shards. The
//! `Stats` threaded into the extractors here is a scratch value, discarded.

use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::{read_fastq_paired_end_state_par_ids, Merge};
use bioreader::sequence::fastq_record::RefFastqRecord;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

use crate::align::common::{KmerExtractor, RangeExtractor};
use crate::align::process::kmer_extractor::StdKmerExtractor;
use crate::align::process::range_extractor::StdRangeExtractor;
use crate::align::stats::Stats;
use crate::database::common::FlexalignDatabase;

use super::emit::{emit_ranges, EmitStats};
use super::record::RangeRecord;

/// One (read, mate) that had at least one surviving range in this shard, keyed by the read's flat
/// file-order [`ReadId::ordinal`]. The rejoin gathers a read by looking this ordinal up per shard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceGroup {
    pub ordinal: u64,
    pub mate: u8,
    pub ranges: Vec<RangeRecord>,
}

/// Everything one shard pass found, in `(ordinal, mate)` order.
///
/// Groups are sparse: a (read, mate) with no surviving range in this shard contributes nothing,
/// which is the common case. Lookup is a binary search over the sorted groups.
#[derive(Debug, Clone, Default)]
pub struct ShardEvidence {
    /// Sorted by `(ordinal, mate)`.
    pub groups: Vec<EvidenceGroup>,
    pub emit: EmitStats,
    /// Read pairs seen by this pass (both mates count as one). Purely an IO stat.
    pub read_pairs: usize,
}

impl ShardEvidence {
    /// The ranges this shard holds for `(ordinal, mate)`, or `None` if it saw no hit there.
    pub fn ranges_for(&self, ordinal: u64, mate: u8) -> Option<&[RangeRecord]> {
        self.groups
            .binary_search_by_key(&(ordinal, mate), |g| (g.ordinal, g.mate))
            .ok()
            .map(|i| self.groups[i].ranges.as_slice())
    }
}

/// Per-worker accumulator. Groups arrive out of order (workers take batches as they win the lock),
/// so they are sorted once at the end of the pass rather than maintained in order.
#[derive(Clone, Default)]
struct PassState {
    groups: Vec<EvidenceGroup>,
    emit: EmitStats,
    read_pairs: usize,
}

impl Merge for PassState {
    fn merge_from(&mut self, other: &mut Self) {
        self.groups.append(&mut other.groups);
        self.emit.headered += other.emit.headered;
        self.emit.headerless += other.emit.headerless;
        self.emit.discarded += other.emit.discarded;
        self.emit.capped += other.emit.capped;
        self.emit.cells += other.emit.cells;
        self.read_pairs += other.read_pairs;
    }
}

/// Runs one shard pass over a paired FASTQ against `db` (a shard's key-range view of the index) and
/// returns its evidence in RAM.
///
/// `db` is typically a [`ShardedDB`](super::db::ShardedDB), but any `FlexalignDatabase` works -- a
/// view over the whole key space yields the evidence a single unsharded prefilter would, which is
/// what the tests compare against. `read1`/`read2` are the two mate streams (already decompressed if
/// gzipped); `buffer_size` **must** match across every pass and the rejoin so the read ids line up.
pub fn run_shard_pass<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    // `Clone` because the reader clones the per-worker closure, which owns a `StdRangeExtractor`;
    // that type's derived `Clone` bounds `D` even though it only holds a `&D`. Real databases
    // (`DB`, `ShardedDB<DB>`) are `Clone`, so this costs nothing in practice.
    D: FlexalignDatabase + Sync + Clone,
    R: std::io::Read + Send,
>(
    read1: R,
    read2: R,
    db: &D,
    buffer_size: usize,
    num_threads: u32,
    max_best_flex: usize,
    max_headered: usize,
) -> ShardEvidence {
    // Captured by the closure and cloned per worker, so each worker gets its own extractors and its
    // own scratch Stats. The extractors borrow `db`; `thread::scope` inside the reader keeps that
    // borrow valid without requiring `'static`.
    let kmer_extractor = StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default();
    let range_extractor = StdRangeExtractor::<K, C, F, D>::new(db);

    let mut kmer_extractor = kmer_extractor;
    let mut range_extractor = range_extractor;
    let mut scratch = Stats::default();
    let mut dist_buf: Vec<u8> = Vec::new();

    let worker = move |rec1: &RefFastqRecord,
                       rec2: &RefFastqRecord,
                       id: ReadId,
                       state: &mut PassState| {
        for (mate, rec) in [(0u8, rec1), (1u8, rec2)] {
            let kmers = kmer_extractor.generate(rec, &mut scratch);
            let ranges = range_extractor.generate(kmers, rec, &mut scratch);

            let mut records = Vec::new();
            let est = emit_ranges::<F>(ranges, max_best_flex, max_headered, &mut dist_buf, &mut records);

            state.emit.headered += est.headered;
            state.emit.headerless += est.headerless;
            state.emit.discarded += est.discarded;
            state.emit.capped += est.capped;
            state.emit.cells += est.cells;

            if !records.is_empty() {
                state.groups.push(EvidenceGroup { ordinal: id.ordinal, mate, ranges: records });
            }
        }
        state.read_pairs += 1;
    };

    let state = read_fastq_paired_end_state_par_ids(
        read1,
        read2,
        buffer_size,
        num_threads,
        PassState::default(),
        worker,
    )
    .expect("paired fastq read");

    let mut groups = state.groups;
    // A read-mate is emitted at most once, so `(ordinal, mate)` is a unique, total key.
    groups.sort_unstable_by_key(|g| (g.ordinal, g.mate));

    ShardEvidence { groups, emit: state.emit, read_pairs: state.read_pairs }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::shard::db::{ShardRange, ShardedDB};
    use crate::shard::test_util::{make_pair, MockDb, C, F, K, KEY_SPACE, L, S};

    fn pass<D: FlexalignDatabase + Sync + Clone>(
        d1: &[u8],
        d2: &[u8],
        db: &D,
        buffer_size: usize,
        threads: u32,
    ) -> ShardEvidence {
        run_shard_pass::<K, C, F, S, L, D, &[u8]>(d1, d2, db, buffer_size, threads, 32, 15)
    }

    /// Canonical form for comparison: every group's ranges sorted by qpos (equal `n_positions` makes
    /// the extractor's sort order arbitrary), groups already sorted by `(ordinal, mate)`.
    fn canonical(mut ev: ShardEvidence) -> Vec<(u64, u8, Vec<RangeRecord>)> {
        ev.groups
            .iter_mut()
            .map(|g| {
                g.ranges.sort_by_key(|r| r.qpos);
                (g.ordinal, g.mate, g.ranges.clone())
            })
            .collect()
    }

    #[test]
    fn evidence_is_deterministic_across_thread_counts() {
        let (d1, d2) = make_pair(3000);
        let full = ShardedDB::new(MockDb::new(), ShardRange::new(0, KEY_SPACE));

        let one = canonical(pass(&d1, &d2, &full, 512, 1));
        let many = canonical(pass(&d1, &d2, &full, 512, 16));

        assert!(!one.is_empty(), "mock should produce evidence");
        assert_eq!(one, many, "sharded evidence must not depend on thread count");
    }

    #[test]
    fn shards_partition_the_full_pass_evidence() {
        // Union of N key-sharded passes == one full-key-space pass. This is the core property: the
        // cut loses nothing and duplicates nothing (every core-mer lands in exactly one shard).
        let (d1, d2) = make_pair(3000);

        let full_db = ShardedDB::new(MockDb::new(), ShardRange::new(0, KEY_SPACE));
        let full = canonical(pass(&d1, &d2, &full_db, 4096, 4));

        let mut union: Vec<(u64, u8, Vec<RangeRecord>)> = Vec::new();
        for range in ShardRange::split_even(KEY_SPACE, 7) {
            let shard_db = ShardedDB::new(MockDb::new(), range);
            let ev = pass(&d1, &d2, &shard_db, 4096, 4);
            for g in ev.groups {
                union.push((g.ordinal, g.mate, g.ranges));
            }
        }
        // Merge the per-shard fragments of each (ordinal, mate) into one group, then canonicalize.
        union.sort_by_key(|(ord, mate, _)| (*ord, *mate));
        let mut merged: Vec<(u64, u8, Vec<RangeRecord>)> = Vec::new();
        for (ord, mate, mut ranges) in union {
            match merged.last_mut() {
                Some((o, m, rs)) if *o == ord && *m == mate => rs.append(&mut ranges),
                _ => merged.push((ord, mate, ranges)),
            }
        }
        for (_, _, rs) in &mut merged {
            rs.sort_by_key(|r| r.qpos);
        }

        assert_eq!(merged, full, "shard union differs from the full pass");
    }

    /// Not a test — a dumper so the shape of a shard pass's output is inspectable. Run with:
    /// `cargo +nightly test --lib shard::pass::tests::dump_example_evidence -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn dump_example_evidence() {
        let (d1, d2) = make_pair(8);
        let shards = ShardRange::split_even(KEY_SPACE, 2);
        for (i, range) in shards.iter().enumerate() {
            let db = ShardedDB::new(MockDb::new(), *range);
            let ev = pass(&d1, &d2, &db, 4096, 1);
            eprintln!(
                "\n== shard {} (keys {}..{}) ==  read_pairs={} emit={:?}",
                i, range.lo, range.hi, ev.read_pairs, ev.emit
            );
            for g in &ev.groups {
                eprintln!("  read#{} mate{}  ({} ranges)", g.ordinal, g.mate, g.ranges.len());
                for r in &g.ranges {
                    eprintln!(
                        "      qpos={:<3} n_positions={} headered={} min_dist={} cells={:?}",
                        r.qpos, r.n_positions, r.headered, r.min_dist, r.cells
                    );
                }
            }
        }
    }

    #[test]
    fn ordinals_match_file_order_and_are_within_range() {
        let n = 1000;
        let (d1, d2) = make_pair(n);
        let full = ShardedDB::new(MockDb::new(), ShardRange::new(0, KEY_SPACE));
        let ev = pass(&d1, &d2, &full, 700, 8);

        assert_eq!(ev.read_pairs, n, "every pair should be seen once per pass");
        for g in &ev.groups {
            assert!(g.ordinal < n as u64, "ordinal {} out of range", g.ordinal);
            assert!(g.mate == 0 || g.mate == 1);
            for r in &g.ranges {
                assert!(!r.headered, "mock emits only headerless ranges");
                assert_eq!(r.n_positions, 1);
            }
        }
    }
}
