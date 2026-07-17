//! A key-range view over an existing index (ProjectShard.md §2).
//!
//! `FMKeys` is direct-addressed over all `2^(2C)` core-mers in key order, so a shard is just a key
//! range `[lo, hi)` -- a contiguous slice of both the key and value arrays. This module provides
//! the *logical* half of that: `ShardedDB` answers `get_vrange` only for keys it owns and `None`
//! for everything else, which is all a shard pass needs to behave as though it held a sliced index.
//!
//! Slicing the file physically (so the pass only pays for its own keys) is a separate step; the
//! logical view is what makes the pass testable and is what the in-RAM path (§10) uses.
//!
//! Key-sharding is **bit-exact per key**: a key lives entirely in one shard, so `positions.len()`,
//! `HEADER_THRESHOLD` and the flank headers all see exactly what they see unsharded. That is the
//! property `shards_partition_the_key_space` pins down, and it is the reason every semantic
//! question in ProjectShard is about *cross-key* decisions only.

use flexmap::values::VRange;

use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::options::Options;

/// A half-open key range `[lo, hi)`. Boundaries are multiples of `KEY_ALIGN` so a shard never
/// splits a control block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ShardRange {
    pub lo: u64,
    pub hi: u64,
}

/// `FMKeys` packs 16 keys per control block (4 head cells + 16 body cells per 20-cell block), so
/// shard boundaries must land on a multiple of 16 to keep each block wholly inside one shard.
pub const KEY_ALIGN: u64 = 16;

impl ShardRange {
    pub fn new(lo: u64, hi: u64) -> Self {
        debug_assert!(lo <= hi, "shard range {}..{} is inverted", lo, hi);
        debug_assert!(lo % KEY_ALIGN == 0, "shard lo {} is not block-aligned", lo);
        Self { lo, hi }
    }

    pub fn contains(&self, cmer: u64) -> bool {
        cmer >= self.lo && cmer < self.hi
    }

    pub fn len(&self) -> u64 {
        self.hi - self.lo
    }

    pub fn is_empty(&self) -> bool {
        self.lo == self.hi
    }

    /// Splits `key_space` into `n` block-aligned ranges covering it exactly once.
    ///
    /// This splits by **key count**, which is only a proxy for the thing that matters -- bytes of
    /// value payload resident per pass. Keys are not uniformly populated, so a hot range makes one
    /// shard fatter than its peers. §2 gives the real answer: the control-head array is already a
    /// prefix sum over value offsets, so byte-balanced boundaries are a binary search rather than a
    /// scan. That needs the built index in hand; this split needs nothing, and it is what the
    /// in-RAM path and the tests run on.
    pub fn split_even(key_space: u64, n: usize) -> Vec<ShardRange> {
        assert!(n > 0, "need at least one shard");
        let blocks = key_space.div_ceil(KEY_ALIGN);
        (0..n as u64)
            .map(|i| {
                let lo = (blocks * i / n as u64) * KEY_ALIGN;
                let hi = ((blocks * (i + 1) / n as u64) * KEY_ALIGN).min(key_space);
                ShardRange::new(lo, hi)
            })
            .collect()
    }
}

/// An existing database restricted to one shard's keys.
///
/// Drops into `StdRangeExtractor` unchanged -- it is generic over `D: FlexalignDatabase` -- so a
/// shard pass reuses the real range lookup rather than a parallel implementation of it.
///
/// `build`/`save`/`load` are not supported: this is a **view** over an already-built index, not an
/// index. They panic rather than delegate, because delegating would silently operate on the whole
/// underlying database while the caller believed it held a shard.
#[derive(Debug, Clone)]
pub struct ShardedDB<D> {
    inner: D,
    range: ShardRange,
}

impl<D: FlexalignDatabase> ShardedDB<D> {
    pub fn new(inner: D, range: ShardRange) -> Self {
        Self { inner, range }
    }

    pub fn range(&self) -> ShardRange {
        self.range
    }

    pub fn inner(&self) -> &D {
        &self.inner
    }

    pub fn into_inner(self) -> D {
        self.inner
    }
}

impl<D: FlexalignDatabase> FlexalignDatabase for ShardedDB<D> {
    /// The whole point: a key outside this shard reads as absent, exactly as it would if the
    /// shard's slice of the index were the only thing loaded.
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange<'_>> {
        if !self.range.contains(canonical_kmer) {
            return None;
        }
        self.inner.get_vrange(canonical_kmer)
    }

    // References are not sharded (ProjectShard.md scope) -- every shard sees all of them, and the
    // rejoin's reference-based stages need them intact.
    fn get_rid(&self, reference: &str) -> Option<&usize> {
        self.inner.get_rid(reference)
    }

    fn get_rname(&self, id: usize) -> Option<&str> {
        self.inner.get_rname(id)
    }

    fn get_reference(&self, id: usize) -> Option<&[u8]> {
        self.inner.get_reference(id)
    }

    fn build(_options: &Options) -> Self {
        unimplemented!("ShardedDB is a view over a built index; build the inner database and wrap it")
    }

    fn save(&self, _paths: &DBPaths, _version: u32) -> Result<(), std::io::Error> {
        unimplemented!("ShardedDB is a view; save the inner database instead")
    }

    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!("ShardedDB is a view; load the inner database and wrap it with a ShardRange")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_covers_the_key_space_exactly_once() {
        for n in [1usize, 2, 3, 7, 10, 16, 64] {
            let space = 1u64 << 20;
            let shards = ShardRange::split_even(space, n);
            assert_eq!(shards.len(), n);
            assert_eq!(shards[0].lo, 0, "n={}", n);
            assert_eq!(shards.last().unwrap().hi, space, "n={}", n);
            for w in shards.windows(2) {
                assert_eq!(w[0].hi, w[1].lo, "gap or overlap at n={}", n);
            }
            assert_eq!(shards.iter().map(|s| s.len()).sum::<u64>(), space, "n={}", n);
        }
    }

    #[test]
    fn split_boundaries_are_block_aligned() {
        // A boundary inside a 16-key control block would put one block in two shards.
        for n in [1usize, 3, 7, 10, 64] {
            for s in ShardRange::split_even(1u64 << 20, n) {
                assert_eq!(s.lo % KEY_ALIGN, 0, "lo {} unaligned at n={}", s.lo, n);
            }
        }
    }

    #[test]
    fn shards_partition_the_key_space() {
        // The load-bearing property: every key is owned by exactly one shard, so a sharded run
        // sees each key's block exactly as an unsharded run does.
        let space = 4096u64;
        let shards = ShardRange::split_even(space, 7);
        for k in 0..space {
            let owners = shards.iter().filter(|s| s.contains(k)).count();
            assert_eq!(owners, 1, "key {} owned by {} shards", k, owners);
        }
    }

    #[test]
    fn split_handles_more_shards_than_blocks() {
        // 64 keys = 4 blocks, but 8 shards: some must come out empty rather than overlapping.
        let shards = ShardRange::split_even(64, 8);
        assert_eq!(shards.iter().map(|s| s.len()).sum::<u64>(), 64);
        assert!(shards.iter().any(|s| s.is_empty()));
        for k in 0..64u64 {
            assert_eq!(shards.iter().filter(|s| s.contains(k)).count(), 1, "key {}", k);
        }
    }

    #[test]
    fn range_is_half_open() {
        let r = ShardRange::new(16, 32);
        assert!(!r.contains(15));
        assert!(r.contains(16));
        assert!(r.contains(31));
        assert!(!r.contains(32), "hi must be exclusive");
    }
}
