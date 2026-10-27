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
//!
//! A **compact** index (`--index-format compact`, the default) is sliced by flexmap's own
//! [`PackedFlexmap::slice`]: keys `[lo, hi)` become keys `0..hi - lo` of a self-contained compact
//! index holding only their values (`.shardN.compact`). Its key table holds the used keys alone, so
//! boundaries need no block alignment; they are still multiples of `KEY_ALIGN` so one boundary rule
//! serves both formats. A compact shard is read into memory on load and freed when dropped.
//!
//! A **hashed** index (`--index-format hashed`) is sliced by flexmap's
//! [`PackedFlexmap::slice_hashed`]: keys `[lo, hi)` become keys `0..hi - lo` of a self-contained
//! hashed index (`.shardN.hashed`), whose key table homes a key by multiply-high and so takes a key
//! space of any size. It is read into memory on load, as a compact shard is.

use std::sync::Arc;

use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use flexmap::keys::FMKeys;
use flexmap::packed::{CompactKeys, HashedKeys, PackedFlexmap, SlotLocator};
use flexmap::range::AnyRange;
use flexmap::values::FMValues;
use serde::{Deserialize, Serialize};

use crate::database::common::{DBPaths, FlexalignDatabase, IndexFormat};
use crate::database::flexmap::FlexmapStore;
use crate::database::hashed_lines::HashedLines;
use crate::options::Options;

use super::db::{ShardRange, KEY_ALIGN};

/// A control-block head is a `u64` stored as four `u16` `KCell`s (fixed by the flexmap layout;
/// `FMKeys::CELLS_PER_HEAD` is private, hence mirrored here).
const CELLS_PER_HEAD: usize = 4;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ShardEntry {
    pub lo: u64,
    pub hi: u64,
    /// Self-contained shard index (rebased key slice + this shard's values), in the manifest's
    /// format. Called `blob_file` in manifests written before packed shards.
    #[serde(alias = "blob_file")]
    pub file: String,
}

fn blob_format() -> IndexFormat {
    IndexFormat::Blob
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
    /// The shards' index format. Manifests written before packed shards hold blobs.
    #[serde(default = "blob_format")]
    pub format: IndexFormat,
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

    /// Shard counts already sliced for this reference, ascending.
    ///
    /// The slicings are discoverable from the reference path alone -- `<index>.s<n>.shards.json` --
    /// so a caller can tell whether a sharded index exists without being told. That is what lets
    /// the main binary pick the sharded path automatically instead of requiring a separate
    /// executable, and what lets several shard counts coexist for the same reference.
    pub fn available(paths: &DBPaths) -> Vec<usize> {
        let stem = format!("{}.s", paths.index_path.to_string_lossy());
        let dir = match paths.index_path.parent() {
            Some(d) => d,
            None => return Vec::new(),
        };
        let mut found: Vec<usize> = std::fs::read_dir(dir)
            .into_iter()
            .flatten()
            .flatten()
            .filter_map(|e| {
                let p = e.path();
                let s = p.to_string_lossy().into_owned();
                let rest = s.strip_prefix(&stem)?.strip_suffix(".shards.json")?;
                rest.parse::<usize>().ok()
            })
            .collect();
        found.sort_unstable();
        found.dedup();
        found
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

/// Byte-balanced shard boundaries for a packed index: as [`split_by_value_bytes`], with the
/// cumulative cell count walked from the key table (one pass over the used keys) instead of read
/// from control heads. Boundaries are rounded down to `KEY_ALIGN` and kept strictly increasing.
fn split_packed_by_cells<K: SlotLocator, const C: usize, const F: usize, const T: usize>(
    full: &PackedFlexmap<K, C, F, T>,
    key_space: u64,
    n: usize,
) -> Vec<ShardRange> {
    assert!(n > 0, "need at least one shard");
    let total = full.values.total_slots();
    let mut bounds = vec![0u64];
    let mut cells_before = 0u64;
    full.for_each_key_in(0, key_space, &mut |key, _, slots| {
        let i = bounds.len() as u64;
        if i < n as u64 && cells_before >= (total as u128 * i as u128 / n as u128) as u64 {
            let b = key / KEY_ALIGN * KEY_ALIGN;
            if b > *bounds.last().unwrap() {
                bounds.push(b);
            }
        }
        cells_before += slots;
    });
    bounds.push(key_space);
    bounds.windows(2).map(|w| ShardRange::new(w[0], w[1])).collect()
}

/// Cuts a packed index into `n_shards` shard files named after `stem`: compact ones, or hashed
/// ones when `hashed` (the format of the index sliced).
fn slice_packed<K: SlotLocator, const C: usize, const F: usize, const T: usize>(
    full: &PackedFlexmap<K, C, F, T>,
    stem: &str,
    n_shards: usize,
    hashed: bool,
) -> std::io::Result<Vec<ShardEntry>> {
    let key_space = 1u64 << (2 * C);
    let mut shards = Vec::new();
    for (i, range) in split_packed_by_cells(full, key_space, n_shards).into_iter().enumerate() {
        if range.is_empty() {
            continue;
        }
        let file = if hashed {
            let file = format!("{}.shard{}.hashed", stem, i);
            full.slice_hashed(range.lo, range.hi)
                .ok_or_else(|| std::io::Error::other(format!("shard {i}: its hashed key table does not fit 64-bit entries")))?
                .save(&file)?;
            file
        } else {
            let file = format!("{}.shard{}.compact", stem, i);
            full.slice(range.lo, range.hi).save(&file)?;
            file
        };
        shards.push(ShardEntry { lo: range.lo, hi: range.hi, file });
    }
    Ok(shards)
}

/// Slices the built index of `paths` into `n_shards` self-contained shard indexes (keys **and**
/// values) of `format`, writes a [`SliceManifest`], and returns it. Loads the full
/// index once (the one-time offline slicing cost); each shard is then loaded on its own.
pub fn slice_index<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(
    paths: &DBPaths,
    n_shards: usize,
    format: IndexFormat,
) -> std::io::Result<SliceManifest> {
    let stem = SliceManifest::stem_for(paths, n_shards);
    let shards = match format {
        IndexFormat::Blob => slice_blob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(paths, n_shards)?,
        IndexFormat::Compact | IndexFormat::Hashed => {
            // The packed index, converted from the blob first if that is all there is.
            let (full, note) = FlexmapStore::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::load(paths, format);
            eprintln!("Slicing the {format:?} index ({note})");
            match &full {
                FlexmapStore::Compact(p) => slice_packed(&**p, &stem, n_shards, false)?,
                FlexmapStore::Hashed(p, _) => slice_packed(&**p, &stem, n_shards, true)?,
                _ => unreachable!("FlexmapStore::load of a packed format returns a packed index"),
            }
        }
    };

    let manifest = SliceManifest {
        format,
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

/// The blob slicing: rebased control-block ranges of the split key/value files.
fn slice_blob<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(
    paths: &DBPaths,
    n_shards: usize,
) -> std::io::Result<Vec<ShardEntry>> {
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

        let file = format!("{}.shard{}.blob", stem, i);
        shard_map.save_blob(&file);
        shards.push(ShardEntry { lo: range.lo, hi: range.hi, file });
    }

    Ok(shards)
}

/// A shard's index: a memory-mapped blob, or a compact or hashed index read into memory (a hashed
/// one with flexalign's reader of its key table, as the unsharded index).
#[derive(Clone)]
enum ShardStore<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize> {
    Blob(FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>),
    Compact(Arc<PackedFlexmap<CompactKeys, C, F, HEADER_THRESHOLD>>),
    Hashed(Arc<PackedFlexmap<HashedKeys, C, F, HEADER_THRESHOLD>>, HashedLines),
}

/// One physical shard: a memory-mapped blob or a compact index holding this shard's rebased keys and
/// its own value slice. Answers `get_vrange` only for its `[lo, hi)` key range. The shard is
/// released when dropped, so a pass that finishes with a shard frees it immediately.
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
    store: ShardStore<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
    lo: u64,
    hi: u64,
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    PhysicalShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    /// Opens shard `shard_idx`: memory-maps a blob (near-instant; nothing is read until queried) or
    /// reads a compact shard into memory on `threads` threads.
    pub fn load(manifest: &SliceManifest, shard_idx: usize, threads: usize) -> Self {
        let entry = &manifest.shards[shard_idx];
        let store = match manifest.format {
            IndexFormat::Blob => ShardStore::Blob(FlexmapBlob::mmap_from_file(&entry.file)),
            IndexFormat::Compact => ShardStore::Compact(Arc::new(
                PackedFlexmap::load(&entry.file, threads).unwrap_or_else(|e| panic!("read {}: {e}", entry.file)),
            )),
            IndexFormat::Hashed => {
                let p: PackedFlexmap<HashedKeys, C, F, HEADER_THRESHOLD> =
                    PackedFlexmap::load(&entry.file, threads).unwrap_or_else(|e| panic!("read {}: {e}", entry.file));
                let lines = HashedLines::new(&p.keys);
                ShardStore::Hashed(Arc::new(p), lines)
            }
        };
        Self { store, lo: entry.lo, hi: entry.hi }
    }

    pub fn range(&self) -> ShardRange {
        ShardRange::new(self.lo, self.hi)
    }
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexalignDatabase for PhysicalShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    fn get_vrange(&self, canonical_kmer: u64) -> Option<flexmap::range::AnyRange<'_>> {
        if canonical_kmer < self.lo || canonical_kmer >= self.hi {
            return None;
        }
        // `lo` is block-aligned, so `k - lo` selects the same block/slot in a blob slice, and the
        // rebased heads make the resulting offsets index this shard's own value slice. A packed
        // slice holds keys `0..hi - lo` outright.
        let key = canonical_kmer - self.lo;
        match &self.store {
            ShardStore::Blob(b) => b.get_vrange(key).map(AnyRange::from),
            ShardStore::Compact(p) => p.get(key).map(AnyRange::from),
            ShardStore::Hashed(p, h) => h.locate(key).map(|(start, slots)| AnyRange::from(p.values.range(start, slots))),
        }
    }

    /// A compact shard fetches a read's in-range keys ahead, as the unsharded index does. Keys
    /// outside the shard are answered `None` without touching it (a packed table must not be asked
    /// past its keys).
    fn for_each_vrange<'s>(&'s self, keys: &[u64], f: &mut dyn FnMut(usize, Option<AnyRange<'s>>)) {
        fn prefetched<'s, K: SlotLocator, const C: usize, const F: usize, const T: usize>(
            p: &'s PackedFlexmap<K, C, F, T>,
            keys: &[u64],
            lo: u64,
            hi: u64,
            f: &mut dyn FnMut(usize, Option<AnyRange<'s>>),
        ) {
            let mut inside = Vec::with_capacity(keys.len());
            let mut at = Vec::with_capacity(keys.len());
            for (i, &k) in keys.iter().enumerate() {
                if k >= lo && k < hi {
                    inside.push(k - lo);
                    at.push(i);
                }
            }
            // In key order: the keys outside the shard before each inside one, as `None`.
            let mut next = 0;
            p.for_each_prefetched(&inside, |j, r| {
                for i in next..at[j] {
                    f(i, None);
                }
                f(at[j], r.map(AnyRange::from));
                next = at[j] + 1;
            });
            for i in next..keys.len() {
                f(i, None);
            }
        }
        match &self.store {
            ShardStore::Compact(p) => prefetched(&**p, keys, self.lo, self.hi, f),
            ShardStore::Hashed(p, h) => {
                // As `prefetched`, through the hashed reader: the in-range keys, rebased, then the
                // keys outside the shard as `None`, all in key order.
                let mut inside = Vec::with_capacity(keys.len());
                let mut at = Vec::with_capacity(keys.len());
                for (i, &k) in keys.iter().enumerate() {
                    if k >= self.lo && k < self.hi {
                        inside.push(k - self.lo);
                        at.push(i);
                    }
                }
                let mut next = 0;
                h.for_each(&**p, &inside, |j, r| {
                    for i in next..at[j] {
                        f(i, None);
                    }
                    f(at[j], r.map(AnyRange::from));
                    next = at[j] + 1;
                });
                for i in next..keys.len() {
                    f(i, None);
                }
            }
            ShardStore::Blob(_) => {
                for (i, &k) in keys.iter().enumerate() {
                    f(i, self.get_vrange(k));
                }
            }
        }
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
        unimplemented!("slice_index writes shard files")
    }
    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!("use PhysicalShardDB::load(manifest, shard_idx)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifests_from_before_packed_shards_read_as_blob() {
        let old = r#"{"c":15,"f":16,"cells_per_body":16,"header_threshold":2,
            "id2ref_file":"r.flex.id2ref","ref2id_file":"r.flex.ref2id",
            "shards":[{"lo":0,"hi":16,"blob_file":"r.flex.index.s1.shard0.blob"}]}"#;
        let m: SliceManifest = serde_json::from_str(old).expect("old manifest");
        assert_eq!(m.format, IndexFormat::Blob);
        assert_eq!(m.shards[0].file, "r.flex.index.s1.shard0.blob");

        let new = serde_json::to_string(&SliceManifest { format: IndexFormat::Compact, ..m }).unwrap();
        assert!(new.contains(r#""format":"compact""#) && new.contains(r#""file":"#), "{new}");
        assert_eq!(serde_json::from_str::<SliceManifest>(&new).unwrap().format, IndexFormat::Compact);
    }
}
