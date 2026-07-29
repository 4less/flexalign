//! Physical index slicing (ProjectShard.md §2, §13.6).
//!
//! Cuts a built index into N self-contained shard blobs so a shard pass loads only its slice.
//! Both halves of the index are sharded:
//!
//! * **Keys** — direct-addressed, uniform per control block, so a shard is a contiguous block range
//!   `[lo_block, hi_block]`. Because [`ShardRange`] boundaries are multiples of `CELLS_PER_BODY`,
//!   a key `k` in `[lo, hi)` maps to the same block/slot as `k - lo` does in the slice, so lookups
//!   reuse flexmap's own `vrange` verbatim on `k - lo`.
//! * **Values** — a contiguous key range references a **contiguous** value range `[v_lo, v_hi)`
//!   (control-head offsets are monotone), so each shard carries only its own values. The sliced
//!   keys' control heads are **rebased** by `-v_lo` so they index the sliced values from 0. This is
//!   where the memory win comes from: values scale with the data, and a shard now holds only its
//!   share of them.
//!
//! Each shard is written as one [`FlexmapBlob`] (`.shardN.blob`) and **memory-mapped** on load, so a
//! pass pages in only what it touches and the mapping is released the moment the shard is dropped.

use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use flexmap::keys::FMKeys;
use flexmap::values::{FMValues, VRange};
use serde::{Deserialize, Serialize};

use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::options::Options;

use super::db::ShardRange;

/// A control-block head is a `u64` stored as four `u16` `KCell`s (fixed by the flexmap layout;
/// `FMKeys::CELLS_PER_HEAD` is private, hence mirrored here).
const CELLS_PER_HEAD: usize = 4;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ShardEntry {
    pub lo: u64,
    pub hi: u64,
    /// Self-contained shard blob: rebased key slice + this shard's value slice.
    pub blob_file: String,
}

/// Records how an index was sliced: the scheme parameters, the shared reference files, and each
/// shard's key range and blob.
///
/// The shard COUNT is part of every filename -- `<index>.s<n>.shards.json` and
/// `<index>.s<n>.shard<i>.blob` -- so slicings at different counts sit side by side next to the
/// index, alongside the unsharded blob, instead of overwriting one another. With the count omitted,
/// re-slicing the same reference at a different N silently clobbered the previous slicing (a
/// 4-shard slice overwrote shard0/shard1 of a 2-shard one and orphaned shard2/shard3), so switching
/// back meant a full re-slice.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SliceManifest {
    pub c: usize,
    pub f: usize,
    pub cells_per_body: u64,
    pub header_threshold: usize,
    pub id2ref_file: String,
    pub ref2id_file: String,
    pub shards: Vec<ShardEntry>,
}

impl SliceManifest {
    /// Filename stem shared by one slicing's artifacts: `<index>.s<n>`.
    pub fn stem_for(paths: &DBPaths, n_shards: usize) -> String {
        format!("{}.s{}", paths.index_path.to_string_lossy(), n_shards)
    }

    pub fn path_for(paths: &DBPaths, n_shards: usize) -> String {
        format!("{}.shards.json", Self::stem_for(paths, n_shards))
    }

    pub fn load(path: &str) -> std::io::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        Ok(serde_json::from_str(&text).expect("valid shard manifest"))
    }
}

/// Byte-balanced shard boundaries: block-aligned key ranges chosen so each shard carries roughly the
/// same number of value bytes.
///
/// The control-head array is a prefix sum of value offsets, so a block range `[lo, hi)` references
/// `head(hi) - head(lo)` values. We binary-search that monotone array for the first block whose
/// cumulative offset reaches each `i/n` fraction of the total. Boundaries land on whole control
/// blocks (16 keys), so no value block is ever split. This replaces an even key-space split, which
/// piles most of the payload into one shard whenever the populated keys are skewed.
fn split_by_value_bytes<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(
    full: &Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
    key_space: u64,
    block_cells: usize,
    shift: usize,
    n: usize,
) -> Vec<ShardRange> {
    assert!(n > 0, "need at least one shard");
    let total_blocks = key_space >> shift; // number of control blocks
    let head = |b: u64| full.keys.get_control_header_value(b as usize * block_cells);
    let total_values = head(total_blocks); // sentinel head == total value payload

    let mut bounds = Vec::with_capacity(n + 1);
    bounds.push(0u64);
    let mut prev = 0u64;
    for i in 1..n as u64 {
        let target = (total_values as u128 * i as u128 / n as u128) as u64;
        // Smallest block `b` in (prev, total_blocks] with head(b) >= target.
        let (mut lo, mut hi) = (prev, total_blocks);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if head(mid) < target {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let b = lo.max(prev + 1).min(total_blocks); // keep boundaries strictly increasing
        bounds.push(b);
        prev = b;
    }
    bounds.push(total_blocks);

    bounds
        .windows(2)
        .map(|w| ShardRange::new((w[0] << shift).min(key_space), (w[1] << shift).min(key_space)))
        .collect()
}

/// Slices the built index of `paths` into `n_shards` self-contained shard blobs (keys **and**
/// values), writes a [`SliceManifest`], and returns it. Loads the full index once (the one-time
/// offline slicing cost); each shard blob is then independently memory-mappable.
pub fn slice_index<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(
    paths: &DBPaths,
    n_shards: usize,
) -> std::io::Result<SliceManifest> {
    let full = Flexmap::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::load(
        &paths.index_keys_file(),
        &paths.index_values_file(),
    );

    let block_cells = CELLS_PER_HEAD + CELLS_PER_BODY as usize;
    let shift = CELLS_PER_BODY.ilog2() as usize;
    let key_space = 1u64 << (2 * C);
    // The shard count is in the stem, so slicings at different counts cannot overwrite each other.
    let stem = SliceManifest::stem_for(paths, n_shards);

    let mut shards = Vec::new();
    for (i, range) in split_by_value_bytes(&full, key_space, block_cells, shift, n_shards)
        .into_iter()
        .enumerate()
    {
        if range.is_empty() {
            continue;
        }
        let lo_block = (range.lo >> shift) as usize;
        let hi_block = ((range.hi - 1) >> shift) as usize;
        let start = lo_block * block_cells;
        // Include the next block's head: `vrange` of the last key reads it for its value-end.
        let end = ((hi_block + 1) * block_cells + CELLS_PER_HEAD).min(full.keys.data.len());

        // Value range this shard references: first block's head offset .. next block's head offset.
        // Control-head offsets are monotone, so this is contiguous.
        let v_lo = full.keys.get_control_header_value(lo_block * block_cells) as usize;
        let v_hi = full.keys.get_control_header_value((hi_block + 1) * block_cells) as usize;

        // Rebased key slice: subtract v_lo from every control head so it indexes the sliced values.
        let mut keys = FMKeys::<C, CELLS_PER_BODY> { data: full.keys.data[start..end].to_vec() };
        let mut head = 0;
        while head + CELLS_PER_HEAD <= keys.data.len() {
            let abs = keys.get_control_header_value(head);
            keys.set_control_header_value(head, abs - v_lo as u64);
            head += block_cells;
        }

        let values =
            FMValues::<F, HEADER_THRESHOLD> { data: full.values.data[v_lo..v_hi].to_vec() };
        let shard_map = Flexmap::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD> { keys, values };

        let blob_file = format!("{}.shard{}.blob", stem, i);
        shard_map.save_blob(&blob_file);
        shards.push(ShardEntry { lo: range.lo, hi: range.hi, blob_file });
    }

    let manifest = SliceManifest {
        c: C,
        f: F,
        cells_per_body: CELLS_PER_BODY,
        header_threshold: HEADER_THRESHOLD,
        id2ref_file: paths.id2reference_path.to_string_lossy().into_owned(),
        ref2id_file: paths.reference2id_path.to_string_lossy().into_owned(),
        shards,
    };

    std::fs::write(
        SliceManifest::path_for(paths, n_shards),
        serde_json::to_string_pretty(&manifest).expect("serialize manifest"),
    )?;

    Ok(manifest)
}

/// One physical shard: a memory-mapped blob holding this shard's rebased keys and its own value
/// slice. Answers `get_vrange` only for its `[lo, hi)` key range. The mapping is released when the
/// shard is dropped, so a pass that finishes with a shard frees it immediately.
///
/// Reference lookups are stubbed: a shard pass stops at the range lookup (the rejoin's
/// reference-comparison stages use the full database). `build`/`save`/`load` panic -- a shard is
/// produced by [`slice_index`] and read via [`PhysicalShardDB::load`].
#[derive(Clone)]
pub struct PhysicalShardDB<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> {
    blob: FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
    lo: u64,
    hi: u64,
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    PhysicalShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    /// Memory-maps shard `shard_idx`'s blob (keys + values). Near-instant; nothing is read until
    /// queried.
    pub fn load(manifest: &SliceManifest, shard_idx: usize) -> Self {
        let entry = &manifest.shards[shard_idx];
        Self {
            blob: FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::mmap_from_file(&entry.blob_file),
            lo: entry.lo,
            hi: entry.hi,
        }
    }

    pub fn range(&self) -> ShardRange {
        ShardRange::new(self.lo, self.hi)
    }
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexalignDatabase for PhysicalShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange<'_>> {
        if canonical_kmer < self.lo || canonical_kmer >= self.hi {
            return None;
        }
        // `lo` is block-aligned, so `k - lo` selects the same block/slot in the slice, and the
        // rebased heads make the resulting offsets index this shard's own value slice.
        self.blob.get_vrange(canonical_kmer - self.lo)
    }

    fn get_rid(&self, _reference: &str) -> Option<&usize> {
        None
    }
    fn get_rname(&self, _id: usize) -> Option<&str> {
        None
    }
    fn get_reference(&self, _id: usize) -> Option<&[u8]> {
        None
    }
    fn build(_options: &Options) -> Self {
        unimplemented!("a physical shard is produced by slice_index, not built")
    }
    fn save(&self, _paths: &DBPaths, _version: u32) -> Result<(), std::io::Error> {
        unimplemented!("slice_index writes shard blobs")
    }
    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!("use PhysicalShardDB::load(manifest, shard_idx)")
    }
}
