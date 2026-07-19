//! The rejoin driver — in-RAM seed reconstruction (ProjectShard.md §13.8, §10).
//!
//! Given the per-shard evidence a set of shard passes produced (see [`ShardEvidence`]), the rejoin
//! walks reads in ordinal order, gathers each read's ranges from every shard, and rebuilds its
//! seeds with [`ReplaySeedExtractor`] — the same `walk_ranges` the unsharded path runs. The result
//! is `&[Seed]` per mate, sorted by `(rval, rpos)`, exactly what the unsharded `StdSeedExtractor`
//! hands the anchor stage, so everything below the cut runs unchanged.
//!
//! This is the in-RAM path (§10): no chunk files and no `GroupReader` cursors — a read's evidence is
//! a per-shard lookup by ordinal. It stops at seeds; feeding them into the reference-comparison
//! stages (`AnchorExtractor` onward) reuses today's `ModularPE` and needs a real index + reference,
//! which is the acceptance/integration step (§13.9), not this one.

use std::time::Duration;

use crate::align::data_structures::common::Seed;

use super::pass::ShardEvidence;
use super::record::GroupRecord;
use super::replay::ReplaySeedExtractor;

/// One read's reconstructed seeds, per mate. Only reads with at least one seed are produced.
#[derive(Debug, Clone)]
pub struct RejoinedRead {
    pub ordinal: u64,
    pub seeds_fwd: Vec<Seed>,
    pub seeds_rev: Vec<Seed>,
}

/// Gathers a read's evidence from every shard into `groups` (one `GroupRecord` per shard/mate that
/// saw it). `idx` carries the ordinal; the in-RAM replay keys on `mate`, not the chunk-relative idx,
/// so any stable value works here.
fn gather_read(shards: &[ShardEvidence], ordinal: u64, groups: &mut Vec<GroupRecord>) {
    groups.clear();
    for shard in shards {
        for mate in [0u8, 1] {
            if let Some(ranges) = shard.ranges_for(ordinal, mate) {
                groups.push(GroupRecord { idx: ordinal as u32, mate, ranges: ranges.to_vec() });
            }
        }
    }
}

/// Rebuilds seeds for every read any shard has evidence for, in ordinal order.
///
/// `max_ranges` is the `--ranges` budget, applied globally by the replay after merging shards (only
/// headered ranges consume it). Single-threaded by design: this is the debuggable in-RAM path;
/// parallelising the rejoin per bioreader batch is a scale-out concern (§13.10).
pub fn run_rejoin_seeds<const K: usize, const C: usize, const F: usize>(
    shards: &[ShardEvidence],
    max_ranges: usize,
) -> Vec<RejoinedRead> {
    let mut ordinals: Vec<u64> =
        shards.iter().flat_map(|s| s.groups.iter().map(|g| g.ordinal)).collect();
    ordinals.sort_unstable();
    ordinals.dedup();

    let mut replay_fwd = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut replay_rev = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut groups = Vec::new();
    let mut sort_time = Duration::ZERO;

    let mut out = Vec::with_capacity(ordinals.len());
    for ordinal in ordinals {
        gather_read(shards, ordinal, &mut groups);

        replay_fwd.gather(&groups, 0);
        let seeds_fwd = replay_fwd.generate(&mut sort_time).to_vec();

        replay_rev.gather(&groups, 1);
        let seeds_rev = replay_rev.generate(&mut sort_time).to_vec();

        if !seeds_fwd.is_empty() || !seeds_rev.is_empty() {
            out.push(RejoinedRead { ordinal, seeds_fwd, seeds_rev });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::{HashMap, HashSet};

    use bioreader::fastq_byte_reader::ReadId;
    use bioreader::parallel::fastq::{read_fastq_paired_end_state_par_ids, Merge};
    use bioreader::sequence::fastq_record::RefFastqRecord;
    use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

    use crate::align::common::{KmerExtractor, RangeExtractor, SeedExtractor};
    use crate::align::process::kmer_extractor::StdKmerExtractor;
    use crate::align::process::range_extractor::StdRangeExtractor;
    use crate::align::process::seed_extractor::StdSeedExtractor;
    use crate::align::stats::Stats;
    use crate::shard::db::{ShardRange, ShardedDB};
    use crate::shard::pass::run_shard_pass;
    use crate::shard::test_util::{make_pair, MockDb, C, F, K, KEY_SPACE, L, S};

    const MAX_RANGES: usize = 1000; // high enough never to bind on this mock (all headerless)

    // A total, order-independent identity for a seed (Seed has no Eq/Ord).
    type SeedKey = (u64, u64, u32, u8, u8, u8);

    fn seed_key(s: &Seed) -> SeedKey {
        (s.rpos, s.rval, s.qpos, s.mismatch, s.length, s.flag)
    }

    fn multiset(seeds: &[Seed]) -> Vec<SeedKey> {
        let mut v: Vec<SeedKey> = seeds.iter().map(seed_key).collect();
        v.sort_unstable();
        v
    }

    /// Runs the *unsharded* front-half (kmer -> range -> `StdSeedExtractor`) over the whole index,
    /// single-threaded, returning each read's seeds per mate keyed by ordinal. This is the oracle
    /// the rejoin must reproduce.
    fn unsharded_seeds(d1: &[u8], d2: &[u8]) -> HashMap<u64, (Vec<SeedKey>, Vec<SeedKey>)> {
        #[derive(Default, Clone)]
        struct Collector {
            rows: Vec<(u64, Vec<SeedKey>, Vec<SeedKey>)>,
        }
        impl Merge for Collector {
            fn merge_from(&mut self, other: &mut Self) {
                self.rows.append(&mut other.rows);
            }
        }

        let db = MockDb::new();
        let mut kmer = StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default();
        let mut range = StdRangeExtractor::<K, C, F, MockDb>::new(&db);
        let mut seed = StdSeedExtractor::<K, C, F>::new(32, 1 << 20, 0, MAX_RANGES);
        let mut scratch = Stats::default();

        let collected = read_fastq_paired_end_state_par_ids(
            d1,
            d2,
            4096,
            1, // single thread: simplest way to attribute seeds to an ordinal
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

    fn shard_evidence(d1: &[u8], d2: &[u8], n_shards: usize, threads: u32) -> Vec<ShardEvidence> {
        ShardRange::split_even(KEY_SPACE, n_shards)
            .into_iter()
            .map(|range| {
                let db = ShardedDB::new(MockDb::new(), range);
                run_shard_pass::<K, C, F, S, L, ShardedDB<MockDb>, &[u8]>(
                    d1, d2, &db, 4096, threads, 32, 15,
                )
            })
            .collect()
    }

    #[test]
    fn rejoin_seeds_match_the_unsharded_pipeline() {
        // The seed-level acceptance (§13.9): the cut + rejoin reconstructs exactly the seeds the
        // unsharded pipeline would build, for every read and both mates.
        let (d1, d2) = make_pair(3000);

        let oracle = unsharded_seeds(&d1, &d2);
        let evidence = shard_evidence(&d1, &d2, 6, 4);
        let rejoined = run_rejoin_seeds::<K, C, F>(&evidence, MAX_RANGES);

        assert!(!rejoined.is_empty(), "the rejoin produced no seeds");

        for rr in &rejoined {
            let (want_fwd, want_rev) =
                oracle.get(&rr.ordinal).expect("rejoined read absent from unsharded run");
            assert_eq!(&multiset(&rr.seeds_fwd), want_fwd, "fwd seeds differ at ordinal {}", rr.ordinal);
            assert_eq!(&multiset(&rr.seeds_rev), want_rev, "rev seeds differ at ordinal {}", rr.ordinal);
        }

        // And nothing the unsharded run found is missing from the rejoin.
        let rejoined_ordinals: HashSet<u64> = rejoined.iter().map(|r| r.ordinal).collect();
        for (ordinal, (fwd, rev)) in &oracle {
            if !fwd.is_empty() || !rev.is_empty() {
                assert!(
                    rejoined_ordinals.contains(ordinal),
                    "ordinal {} had seeds unsharded but none after rejoin",
                    ordinal
                );
            }
        }
    }

    #[test]
    fn rejoin_is_invariant_to_shard_count() {
        // The rejoin's output must not depend on how the key space was partitioned.
        let (d1, d2) = make_pair(2000);

        let a = run_rejoin_seeds::<K, C, F>(&shard_evidence(&d1, &d2, 1, 4), MAX_RANGES);
        let b = run_rejoin_seeds::<K, C, F>(&shard_evidence(&d1, &d2, 13, 4), MAX_RANGES);

        let as_map = |v: &[RejoinedRead]| -> HashMap<u64, (Vec<SeedKey>, Vec<SeedKey>)> {
            v.iter().map(|r| (r.ordinal, (multiset(&r.seeds_fwd), multiset(&r.seeds_rev)))).collect()
        };
        assert_eq!(as_map(&a), as_map(&b), "rejoin output changed with shard count");
    }
}
