//! Seed-level acceptance: sharded-rejoin seeds vs. the unsharded pipeline on a real index
//! (ProjectShard.md §13.9).
//!
//! With exactness deliberately relaxed (§4), this seed diff *is* the spec: the cut is accepted when
//! the only reads whose seeds differ are the ones the design predicts -- reads where the global
//! `--ranges` budget binds on a tie (§3.3) or the unsharded recovery path fires (§4). This module
//! runs both sides over the same reads and the same index and returns the diff for inspection; it
//! does not assert, because "identical" is not the contract on real data.
//!
//! Both sides share `walk_ranges` (hurdle 11); the difference under test is purely the cut --
//! splitting ranges across key-shards, writing them to the evidence format, and rebuilding them --
//! not two implementations of the seeding.

use std::collections::HashMap;

use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::{read_fastq_paired_end_state_par_ids, Merge};
use bioreader::sequence::fastq_record::RefFastqRecord;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

use crate::align::common::{KmerExtractor, RangeExtractor, SeedExtractor};
use crate::align::data_structures::common::Seed;
use crate::align::process::kmer_extractor::StdKmerExtractor;
use crate::align::process::range_extractor::StdRangeExtractor;
use crate::align::process::seed_extractor::StdSeedExtractor;
use crate::align::stats::Stats;
use crate::database::common::FlexalignDatabase;

use super::db::{ShardRange, ShardedDB};
use super::emit::EmitStats;
use super::pass::{run_shard_pass, ShardEvidence};
use super::rejoin::run_rejoin_seeds;

/// A total, order-independent identity for a seed (`Seed` has no `Eq`/`Ord`).
type SeedKey = (u64, u64, u32, u8, u8, u8);

fn seed_key(s: &Seed) -> SeedKey {
    (s.rpos, s.rval, s.qpos, s.mismatch, s.length, s.flag)
}

fn multiset(seeds: &[Seed]) -> Vec<SeedKey> {
    let mut v: Vec<SeedKey> = seeds.iter().map(seed_key).collect();
    v.sort_unstable();
    v
}

/// The outcome of comparing the two paths' seeds, read by read.
#[derive(Debug, Clone, Default)]
pub struct AcceptanceReport {
    /// Reads that produced at least one seed in either path.
    pub reads_with_seeds: usize,
    /// Reads whose fwd *and* rev seed multisets are identical between the paths.
    pub identical: usize,
    /// Reads that differ, with a few ordinals for inspection.
    pub differing: usize,
    pub differing_examples: Vec<u64>,
    pub unsharded_seed_total: usize,
    pub sharded_seed_total: usize,
    /// Emit counters summed over the shard passes.
    pub emit: EmitStats,
}

impl AcceptanceReport {
    pub fn agreement(&self) -> f64 {
        if self.reads_with_seeds == 0 {
            1.0
        } else {
            self.identical as f64 / self.reads_with_seeds as f64
        }
    }
}

/// Runs the unsharded front-half (kmer -> range -> `StdSeedExtractor`) over the whole index,
/// single-threaded, returning each read's seed multisets per mate keyed by ordinal.
fn unsharded_seeds<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    D: FlexalignDatabase + Sync + Clone,
>(
    read1: &[u8],
    read2: &[u8],
    db: &D,
    buffer_size: usize,
    max_best_flex: usize,
    max_range_size: usize,
    min_ranges: usize,
    max_ranges: usize,
) -> HashMap<u64, (Vec<SeedKey>, Vec<SeedKey>)> {
    #[derive(Default, Clone)]
    struct Collector {
        rows: Vec<(u64, Vec<SeedKey>, Vec<SeedKey>)>,
    }
    impl Merge for Collector {
        fn merge_from(&mut self, other: &mut Self) {
            self.rows.append(&mut other.rows);
        }
    }

    let mut kmer = StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default();
    let mut range = StdRangeExtractor::<K, C, F, D>::new(db);
    let mut seed = StdSeedExtractor::<K, C, F>::new(max_best_flex, max_range_size, min_ranges, max_ranges, 0, 0);
    let mut scratch = Stats::default();

    let collected = read_fastq_paired_end_state_par_ids(
        read1,
        read2,
        buffer_size,
        1, // single thread: attribute seeds to an ordinal without merging bookkeeping
        Collector::default(),
        move |rec1: &RefFastqRecord, rec2: &RefFastqRecord, id: ReadId, st: &mut Collector| {
            let k1 = kmer.generate(rec1, &mut scratch);
            let r1 = range.generate(k1, rec1, &mut scratch);
            let s1 = multiset(seed.generate(r1, &mut scratch));

            let k2 = kmer.generate(rec2, &mut scratch);
            let r2 = range.generate(k2, rec2, &mut scratch);
            let s2 = multiset(seed.generate(r2, &mut scratch));

            st.rows.push((id.ordinal, s1, s2));
        },
    )
    .unwrap();

    collected.rows.into_iter().map(|(o, f, r)| (o, (f, r))).collect()
}

type SeedMap = HashMap<u64, (Vec<SeedKey>, Vec<SeedKey>)>;

/// Diffs the two seed maps read by read, folding in the shard passes' emit counters.
fn diff(unsharded: &SeedMap, sharded: &SeedMap, emit: EmitStats) -> AcceptanceReport {
    let empty = (Vec::new(), Vec::new());
    let mut report = AcceptanceReport { emit, ..Default::default() };

    let mut ordinals: Vec<u64> = unsharded.keys().chain(sharded.keys()).copied().collect();
    ordinals.sort_unstable();
    ordinals.dedup();

    for ordinal in ordinals {
        let u = unsharded.get(&ordinal).unwrap_or(&empty);
        let s = sharded.get(&ordinal).unwrap_or(&empty);

        report.unsharded_seed_total += u.0.len() + u.1.len();
        report.sharded_seed_total += s.0.len() + s.1.len();

        if u.0.is_empty() && u.1.is_empty() && s.0.is_empty() && s.1.is_empty() {
            continue;
        }
        report.reads_with_seeds += 1;

        if u.0 == s.0 && u.1 == s.1 {
            report.identical += 1;
        } else {
            report.differing += 1;
            if report.differing_examples.len() < 20 {
                report.differing_examples.push(ordinal);
            }
        }
    }
    report
}

/// Runs one shard pass per provided shard DB and rejoins the evidence into a per-read seed map,
/// returning it alongside the summed emit counters. `shards` may be `ShardedDB` views or physical
/// `PhysicalShardDB`s -- anything implementing the database trait.
fn sharded_seed_map<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    SD: FlexalignDatabase + Sync + Clone,
>(
    read1: &[u8],
    read2: &[u8],
    shards: &[SD],
    buffer_size: usize,
    threads: u32,
    max_best_flex: usize,
    max_ranges: usize,
    max_headered: usize,
) -> (SeedMap, EmitStats) {
    let mut emit = EmitStats::default();
    let evidence: Vec<ShardEvidence> = shards
        .iter()
        .map(|shard| {
            let ev = run_shard_pass::<K, C, F, S, L, SD, &[u8]>(
                read1, read2, shard, buffer_size, threads, max_best_flex, max_headered,
            );
            emit.headered += ev.emit.headered;
            emit.headerless += ev.emit.headerless;
            emit.discarded += ev.emit.discarded;
            emit.capped += ev.emit.capped;
            emit.cells += ev.emit.cells;
            ev
        })
        .collect();

    let map = run_rejoin_seeds::<K, C, F>(&evidence, max_ranges)
        .into_iter()
        .map(|r| (r.ordinal, (multiset(&r.seeds_fwd), multiset(&r.seeds_rev))))
        .collect();

    (map, emit)
}

/// Compares sharded-rejoin seeds against unsharded seeds over the same reads and index, sharding the
/// index into `n_shards` in-RAM [`ShardedDB`] views that borrow the one index (`ShardedDB<&D>`, no
/// per-shard copy). Interpret the diff against §4 rather than expecting zero differences.
pub fn compare_seeds<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    D: FlexalignDatabase + Sync + Clone,
>(
    read1: &[u8],
    read2: &[u8],
    db: &D,
    n_shards: usize,
    buffer_size: usize,
    threads: u32,
    max_best_flex: usize,
    max_range_size: usize,
    min_ranges: usize,
    max_ranges: usize,
    max_headered: usize,
) -> AcceptanceReport {
    let key_space = 1u64 << (2 * C);
    let views: Vec<ShardedDB<&D>> = ShardRange::split_even(key_space, n_shards)
        .into_iter()
        .map(|range| ShardedDB::new(db, range))
        .collect();

    compare_seeds_with_shards::<K, C, F, S, L, ShardedDB<&D>, D>(
        read1, read2, &views, db, buffer_size, threads, max_best_flex, max_range_size, min_ranges,
        max_ranges, max_headered,
    )
}

/// Like [`compare_seeds`], but the caller supplies the shard databases -- e.g. physical
/// `PhysicalShardDB`s loaded from a slice. `oracle_db` is the full index the unsharded run uses.
pub fn compare_seeds_with_shards<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    SD: FlexalignDatabase + Sync + Clone,
    D: FlexalignDatabase + Sync + Clone,
>(
    read1: &[u8],
    read2: &[u8],
    shards: &[SD],
    oracle_db: &D,
    buffer_size: usize,
    threads: u32,
    max_best_flex: usize,
    max_range_size: usize,
    min_ranges: usize,
    max_ranges: usize,
    max_headered: usize,
) -> AcceptanceReport {
    let (sharded, emit) = sharded_seed_map::<K, C, F, S, L, SD>(
        read1, read2, shards, buffer_size, threads, max_best_flex, max_ranges, max_headered,
    );
    let unsharded = unsharded_seeds::<K, C, F, S, L, D>(
        read1, read2, oracle_db, buffer_size, max_best_flex, max_range_size, min_ranges, max_ranges,
    );
    diff(&unsharded, &sharded, emit)
}
