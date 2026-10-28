//! Online shards (ProjectShard.md §2): a shard is a key range of the index **as it lies in its
//! bundle**, read alone. Nothing is cut, rewritten or copied.
//!
//! * **Bounds** -- `<bundle>.s<n>.bounds.json` holds the `n + 1` keys that cut the key space into
//!   `n` ranges of about equal value bytes. They come from the key table alone (flexmap's
//!   `split_in_file`, about a gigabyte on the full protal DB), are computed once per shard count and
//!   reused while the bundle is unchanged (its length and modification time are stored with them).
//! * **A shard** -- flexmap's `read_window` reads, for keys `[lo, hi)`, the part of the key table
//!   their lookups touch and the values of their cells, each in one parallel read into huge pages
//!   (never a mapping: the bundle may be on a network filesystem). The parts keep the whole index's
//!   addressing, so a lookup is the unsharded lookup, through the same code ([`HashedLines`] for a
//!   hashed table) -- asked only for keys of the range.
//!
//! A packed index (compact or hashed) is sharded this way. A blob index has no bundle and is not
//! sharded.

use std::sync::Arc;
use std::time::UNIX_EPOCH;

use flexmap::packed::{CompactKeys, HashedKeys, PackedFlexmap};
use flexmap::range::AnyRange;
use serde::{Deserialize, Serialize};

use crate::database::bundle::{self, Bundle};
use crate::database::common::{DBPaths, FlexalignDatabase, IndexFormat};
use crate::database::hashed_lines::HashedLines;
use crate::database::repeat::{RepeatIndex, RepeatLookup};
use crate::options::Options;

use super::db::ShardRange;

/// The keys cutting a bundle's index into shards, and the bundle they were computed for.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ShardBounds {
    pub format: IndexFormat,
    pub bundle_len: u64,
    pub bundle_modified_ns: u128,
    /// `n + 1` non-decreasing keys from 0 to the key space.
    pub bounds: Vec<u64>,
}

fn bundle_stamp(bundle: &str) -> std::io::Result<(u64, u128)> {
    let meta = std::fs::metadata(bundle)?;
    let modified = meta.modified()?.duration_since(UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    Ok((meta.len(), modified))
}

impl ShardBounds {
    pub fn path_for(paths: &DBPaths, format: IndexFormat, n_shards: usize) -> String {
        format!("{}.s{}.bounds.json", paths.bundle_file(format), n_shards)
    }

    /// Shard counts with bounds next to the bundle of `format`, ascending.
    pub fn available(paths: &DBPaths, format: IndexFormat) -> Vec<usize> {
        let bundle = paths.bundle_file(format);
        let stem = format!("{bundle}.s");
        let dir = match std::path::Path::new(&bundle).parent() {
            Some(d) => d.to_path_buf(),
            None => return Vec::new(),
        };
        let mut found: Vec<usize> = std::fs::read_dir(dir)
            .into_iter()
            .flatten()
            .flatten()
            .filter_map(|e| {
                let s = e.path().to_string_lossy().into_owned();
                s.strip_prefix(&stem)?.strip_suffix(".bounds.json")?.parse::<usize>().ok()
            })
            .collect();
        found.sort_unstable();
        found.dedup();
        found
    }

    /// The bounds for `n_shards` of the bundle of `format`: read from their file when it is there
    /// and the bundle has not changed since; otherwise computed from the key table and saved.
    /// Returns them and whether they were computed.
    pub fn load_or_compute<const C: usize, const F: usize, const T: usize>(
        paths: &DBPaths,
        format: IndexFormat,
        n_shards: usize,
        threads: usize,
    ) -> std::io::Result<(Self, bool)> {
        let bundle = paths.bundle_file(format);
        let (bundle_len, bundle_modified_ns) = bundle_stamp(&bundle)?;
        let path = Self::path_for(paths, format, n_shards);
        if let Ok(text) = std::fs::read_to_string(&path) {
            if let Ok(b) = serde_json::from_str::<Self>(&text) {
                if b.format == format && b.bundle_len == bundle_len && b.bundle_modified_ns == bundle_modified_ns && b.bounds.len() == n_shards + 1 {
                    return Ok((b, false));
                }
            }
        }
        let (offset, _) = Bundle::section_in_file(&bundle, bundle::INDEX)?;
        let file = std::fs::File::open(&bundle)?;
        let bounds = match format {
            IndexFormat::Hashed => PackedFlexmap::<HashedKeys, C, F, T>::split_in_file(&file, offset, n_shards, threads)?,
            IndexFormat::Compact => PackedFlexmap::<CompactKeys, C, F, T>::split_in_file(&file, offset, n_shards, threads)?,
            IndexFormat::Blob => return Err(std::io::Error::other("a blob index is not sharded")),
        };
        let b = Self { format, bundle_len, bundle_modified_ns, bounds };
        std::fs::write(&path, serde_json::to_string_pretty(&b).expect("serialize shard bounds"))?;
        Ok((b, true))
    }

    /// The non-empty key ranges, in key order.
    pub fn ranges(&self) -> Vec<ShardRange> {
        self.bounds.windows(2).filter(|w| w[0] < w[1]).map(|w| ShardRange::new(w[0], w[1])).collect()
    }
}

/// A shard's index: a window of a compact or hashed index (a hashed one with flexalign's reader
/// of its key table, as the unsharded index).
#[derive(Clone)]
enum ShardStore<const C: usize, const F: usize, const HEADER_THRESHOLD: usize> {
    Compact(Arc<PackedFlexmap<CompactKeys, C, F, HEADER_THRESHOLD>>),
    Hashed(Arc<PackedFlexmap<HashedKeys, C, F, HEADER_THRESHOLD>>, HashedLines),
}

/// One shard: the keys `[lo, hi)` of the bundle's index, read alone ([`ShardBounds`]). Answers
/// `get_vrange` for its keys, as the whole index does, and `None` for every other key without
/// touching memory it does not hold. Freed when dropped, so a pass that finishes with a shard
/// releases it before the next loads.
///
/// Reference lookups are stubbed: a shard pass stops at the range lookup (the join's
/// reference-comparison stages use the reference from the bundle). `build`/`save`/`load` panic.
#[derive(Clone)]
pub struct WindowShardDB<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> {
    store: ShardStore<C, F, HEADER_THRESHOLD>,
    lo: u64,
    hi: u64,
    /// `--repeat-index`: this shard's keys of the repeat index (`RepeatIndex::load_window`).
    repeat: Option<Arc<RepeatIndex>>,
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    WindowShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    /// Reads keys `range` of the index in the bundle of `format`, on `threads` threads.
    pub fn load(paths: &DBPaths, format: IndexFormat, range: ShardRange, threads: usize) -> std::io::Result<Self> {
        let bundle = paths.bundle_file(format);
        let (offset, _) = Bundle::section_in_file(&bundle, bundle::INDEX)?;
        let file = std::fs::File::open(&bundle)?;
        let store = match format {
            IndexFormat::Hashed => {
                let p = PackedFlexmap::<HashedKeys, C, F, HEADER_THRESHOLD>::read_window(&file, offset, range.lo, range.hi, threads)?;
                let lines = HashedLines::new(&p.keys);
                ShardStore::Hashed(Arc::new(p), lines)
            }
            IndexFormat::Compact => ShardStore::Compact(Arc::new(
                PackedFlexmap::<CompactKeys, C, F, HEADER_THRESHOLD>::read_window(&file, offset, range.lo, range.hi, threads)?,
            )),
            IndexFormat::Blob => return Err(std::io::Error::other("a blob index is not sharded")),
        };
        Ok(Self { store, lo: range.lo, hi: range.hi, repeat: None })
    }

    /// Answers the dropped cores of this shard's keys from `repeat` (their window of the repeat
    /// index), as the unsharded index does from the whole of it.
    pub fn with_repeat(mut self, repeat: RepeatIndex) -> Self {
        self.repeat = Some(Arc::new(repeat));
        self
    }

    pub fn range(&self) -> ShardRange {
        ShardRange::new(self.lo, self.hi)
    }

    /// Bytes of the index this shard holds in memory.
    pub fn memory_bytes(&self) -> usize {
        let (k, v) = match &self.store {
            ShardStore::Compact(p) => p.memory_bytes(),
            ShardStore::Hashed(p, _) => p.memory_bytes(),
        };
        k + v + self.repeat.as_ref().map_or(0, |r| r.memory_bytes())
    }
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexalignDatabase for WindowShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    fn get_vrange(&self, canonical_kmer: u64) -> Option<AnyRange<'_>> {
        if canonical_kmer < self.lo || canonical_kmer >= self.hi {
            return None;
        }
        // The window keeps the whole index's addressing: the key is looked up as it is.
        match &self.store {
            ShardStore::Compact(p) => p.get(canonical_kmer).map(AnyRange::from),
            ShardStore::Hashed(p, h) => h.locate(canonical_kmer).map(|(start, slots)| AnyRange::from(p.values.range(start, slots))),
        }
    }

    /// A read's keys of this shard are fetched ahead, as the unsharded index does; the others are
    /// answered `None` without touching the window (it holds nothing for them).
    fn for_each_vrange<'s>(&'s self, keys: &[u64], f: &mut dyn FnMut(usize, Option<AnyRange<'s>>)) {
        let mut inside = Vec::with_capacity(keys.len());
        let mut at = Vec::with_capacity(keys.len());
        for (i, &k) in keys.iter().enumerate() {
            if k >= self.lo && k < self.hi {
                inside.push(k);
                at.push(i);
            }
        }
        // In key order: the keys outside the shard before each inside one, as `None`.
        let mut next = 0;
        let mut emit = |j: usize, r: Option<AnyRange<'s>>| {
            for i in next..at[j] {
                f(i, None);
            }
            f(at[j], r);
            next = at[j] + 1;
        };
        match &self.store {
            ShardStore::Compact(p) => p.for_each_prefetched(&inside, |j, r| emit(j, r.map(AnyRange::from))),
            ShardStore::Hashed(p, h) => h.for_each(&**p, &inside, |j, r| emit(j, r.map(AnyRange::from))),
        }
        for i in next..keys.len() {
            f(i, None);
        }
    }

    fn repeat_lookup_batch<'s>(&'s self, queries: &[(u64, u32)], max: usize, f: &mut dyn FnMut(usize, RepeatLookup<'s>)) {
        // The window holds this shard's dropped cores only: any other key is not a repeat here.
        if let Some(r) = &self.repeat {
            r.lookup_batch(queries, max, f);
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
        unimplemented!("a shard is a window of a built index, not built")
    }
    fn save(&self, _paths: &DBPaths, _version: u32) -> Result<(), std::io::Error> {
        unimplemented!("a shard is a window of the bundle; there is nothing to save")
    }
    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!("use WindowShardDB::load(paths, format, range, threads)")
    }
}
