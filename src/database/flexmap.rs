use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}};

use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use flexmap::packed::{CompactKeys, HashedKeys, PackedFlexmap, SlotLocator};
use flexmap::range::AnyRange;
use std::sync::Arc;
use savefile::{load, save};

use crate::flexalign::time;

use super::bundle::{self, Bundle, NameTable, Names};
use super::refseq2::{self, RefSeq2};
use super::hashed_lines::HashedLines;
use super::common::{DBPaths, load_references, FlexalignDatabase, IndexFormat, index_format};


/// The full-mode index. `Owned` is a freshly built map still on the heap (build-then-align in one
/// run); `Mmap` is the memory-mapped blob loaded from disk -- open is near-instant and only the
/// keys/values actually queried are paged in, so a cold multi-GB index no longer costs a full read.
/// `Compact` and `Hashed` are flexmap's packed index (runs use `Hashed`), read
/// into huge pages on all threads: the same answers in about half the memory, and a lookup pipeline
/// that fetches a read's keys and values ahead. All answer `get_vrange`.
#[derive(Clone)]
pub enum FlexmapStore<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize> {
    Owned(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>),
    Mmap(FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>),
    Compact(Arc<PackedFlexmap<CompactKeys, C, F, HEADER_THRESHOLD>>),
    /// With a faster reader of its key table (`HashedLines`), which shares the table's memory.
    Hashed(Arc<PackedFlexmap<HashedKeys, C, F, HEADER_THRESHOLD>>, HashedLines),
    /// No index: a reference-only load (`DB::load_reference_only`) of a packed index, which has no
    /// blob to map. Every key reads as absent.
    Absent,
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexmapStore<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    #[inline]
    fn get_vrange(&self, canonical_kmer: u64) -> Option<AnyRange<'_>> {
        match self {
            FlexmapStore::Owned(f) => f.get_vrange(canonical_kmer).map(AnyRange::from),
            FlexmapStore::Mmap(b) => b.get_vrange(canonical_kmer).map(AnyRange::from),
            FlexmapStore::Compact(p) => p.get(canonical_kmer).map(AnyRange::from),
            FlexmapStore::Hashed(p, h) => h.locate(canonical_kmer).map(|(start, slots)| AnyRange::from(p.values.range(start, slots))),
            FlexmapStore::Absent => None,
        }
    }

    /// `f(i, get_vrange(keys[i]))` in order, with the index's batch prefetch where it has one.
    #[inline]
    fn for_each_vrange<'s>(&'s self, keys: &[u64], f: &mut dyn FnMut(usize, Option<AnyRange<'s>>)) {
        match self {
            FlexmapStore::Compact(p) => p.for_each_prefetched(keys, |i, r| f(i, r.map(AnyRange::from))),
            FlexmapStore::Hashed(p, h) => h.for_each(p, keys, |i, r| f(i, r.map(AnyRange::from))),
            FlexmapStore::Mmap(b) => {
                let mut i = 0;
                b.for_each_vrange_prefetched(keys, |_, r| {
                    f(i, r.map(AnyRange::from));
                    i += 1;
                });
            }
            FlexmapStore::Owned(m) => {
                for (i, &k) in keys.iter().enumerate() {
                    f(i, m.get_vrange(k).map(AnyRange::from));
                }
            }
            FlexmapStore::Absent => {
                for i in 0..keys.len() {
                    f(i, None);
                }
            }
        }
    }

    fn hashed(index: PackedFlexmap<HashedKeys, C, F, HEADER_THRESHOLD>) -> Self {
        let lines = HashedLines::new(&index.keys);
        FlexmapStore::Hashed(Arc::new(index), lines)
    }

    /// Fetches the key-table lines of `keys`, ahead of a lookup of them.
    #[inline]
    fn prefetch_keys(&self, keys: &[u64]) {
        match self {
            FlexmapStore::Compact(p) => keys.iter().for_each(|&k| p.keys.prefetch_key(k)),
            FlexmapStore::Hashed(_, h) => keys.iter().for_each(|&k| h.prefetch_key(k)),
            _ => {}
        }
    }

    /// The packed index of `format` from this one (a built map or a blob); `self` if it already is one.
    fn into_packed(self, format: IndexFormat) -> Self {
        let compact = match (&self, format) {
            (_, IndexFormat::Blob) | (FlexmapStore::Compact(_), IndexFormat::Compact) | (FlexmapStore::Hashed(..), IndexFormat::Hashed) => return self,
            (FlexmapStore::Owned(m), _) => Arc::new(PackedFlexmap::from_flexmap(m)),
            (FlexmapStore::Mmap(b), _) => Arc::new(PackedFlexmap::from_blob(b)),
            (FlexmapStore::Compact(p), _) => Arc::clone(p),
            (FlexmapStore::Hashed(..), _) => panic!("a hashed index cannot be converted to a compact one; rebuild"),
            (FlexmapStore::Absent, _) => panic!("no index to convert"),
        };
        match format {
            IndexFormat::Compact => FlexmapStore::Compact(compact),
            _ => Self::hashed(PackedFlexmap::from_compact(&compact).expect("the hashed key table's fields fit 64 bits")),
        }
    }

    /// Writes a packed index to `file` (nothing for the blob forms).
    fn save_packed(&self, file: &str) -> std::io::Result<()> {
        match self {
            FlexmapStore::Compact(p) => p.save(file),
            FlexmapStore::Hashed(p, _) => p.save(file),
            _ => Ok(()),
        }
    }

    /// Bytes of the index in memory (key table + values), for the load report.
    fn memory_bytes(&self) -> usize {
        match self {
            FlexmapStore::Owned(m) => m.keys.data.len() * 2 + m.values.data.len() * 8,
            FlexmapStore::Mmap(b) => b.file_size(),
            FlexmapStore::Compact(p) => p.memory_bytes().0 + p.memory_bytes().1,
            FlexmapStore::Hashed(p, _) => p.memory_bytes().0 + p.memory_bytes().1,
            FlexmapStore::Absent => 0,
        }
    }

    /// The index of `format` for `paths`: the packed file when it is there; otherwise the blob,
    /// converted (and the packed file written, so the next run reads it). Returns the index and
    /// what was done, for the load report.
    pub fn load(paths: &DBPaths, format: IndexFormat) -> (Self, String) {
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        let blob_file = paths.index_blob_file();
        if format == IndexFormat::Blob {
            return (FlexmapStore::Mmap(FlexmapBlob::mmap_from_file(&blob_file)), format!("mmap of {blob_file}"));
        }
        let file = paths.index_packed_file(format);
        if Self::packed_compatible(paths, format) {
            let store = match format {
                IndexFormat::Hashed if index_mmap() => Self::hashed(PackedFlexmap::mmap(&file).expect("mappable hashed index")),
                IndexFormat::Hashed => Self::hashed(PackedFlexmap::load(&file, threads).expect("readable hashed index")),
                _ if index_mmap() => FlexmapStore::Compact(Arc::new(PackedFlexmap::mmap(&file).expect("mappable compact index"))),
                _ => FlexmapStore::Compact(Arc::new(PackedFlexmap::load(&file, threads).expect("readable compact index"))),
            };
            return (store, format!("read {file} on {threads} threads"));
        }
        let store = FlexmapStore::Mmap(FlexmapBlob::mmap_from_file(&blob_file)).into_packed(format);
        let note = match store.save_packed(&file) {
            Ok(()) => format!("converted {blob_file}, wrote {file}"),
            Err(e) => {
                log::warn!("Could not write the packed index {file} ({e}); converting again next run.");
                format!("converted {blob_file} (not saved)")
            }
        };
        (store, note)
    }

    /// Whether `paths` holds a packed index of `format` that this build reads.
    pub fn packed_compatible(paths: &DBPaths, format: IndexFormat) -> bool {
        let file = paths.index_packed_file(format);
        match format {
            IndexFormat::Hashed => PackedFlexmap::<HashedKeys, C, F, HEADER_THRESHOLD>::is_compatible(&file),
            _ => PackedFlexmap::<CompactKeys, C, F, HEADER_THRESHOLD>::is_compatible(&file),
        }
    }
}

/// Concatenated reference sequences + an offset table, memory-mapped. `get` returns a zero-copy
/// slice for a reference id, so a DB load can `mmap` this blob instead of re-parsing the multi-GB
/// reference FASTA into owned records every run (the dominant DB-load cost on a large marker DB).
const REFSEQ_MAGIC: [u8; 8] = *b"FREFSEQ1";

pub struct RefSeqBlob {
    // Shared so the DB wrapper stays cheaply cloneable per worker; the raw pointers remain valid
    // views into it.
    backing: RefBacking,
    offsets: *const u64, // n+1 little-endian offsets into the data section; offsets[0] == 0
    data: *const u8,
    n: usize,
}

/// What a [`RefSeqBlob`]'s bytes live in: its own mapped file, or the reference section
/// (offset, length) of a bundle read in or mapped.
#[derive(Clone)]
enum RefBacking {
    File(std::sync::Arc<memmap2::Mmap>),
    Bundle(flexmap::memory::SharedBytes, usize, usize),
}

unsafe impl Send for RefSeqBlob {}
unsafe impl Sync for RefSeqBlob {}
impl Clone for RefSeqBlob {
    fn clone(&self) -> Self {
        Self { backing: self.backing.clone(), offsets: self.offsets, data: self.data, n: self.n }
    }
}
impl RefSeqBlob {
    #[inline]
    fn get(&self, id: usize) -> &[u8] {
        unsafe {
            let start = (*self.offsets.add(id)) as usize;
            let end = (*self.offsets.add(id + 1)) as usize;
            std::slice::from_raw_parts(self.data.add(start), end - start)
        }
    }
    fn len(&self) -> usize { self.n }

    /// Drop this mapping's resident pages (`MADV_DONTNEED`). Read-only and file-backed, so nothing
    /// is lost: a later access re-faults from the page cache. Used by `--ref-budget` to bound the
    /// resident set to one batch's worth of reference instead of everything the run ever touched.
    /// A bundle read into memory has no file behind it, so its pages stay.
    ///
    /// There is deliberately no `prefault` counterpart here: the up-front reference read goes
    /// through `common::warm_file`, which reads the file with several threads (~2 GB/s) rather than
    /// touching mapped pages one at a time on a single thread (~1 GB/s). Same bytes, half the time.
    fn release(&self) {
        // `DontNeed` is `unchecked` in memmap2 because on a writable private mapping it discards
        // un-written changes. These mappings are read-only and file-backed (`Mmap`, never
        // `MmapMut`), so there is nothing to lose: the pages are clean, and a later `get` re-faults
        // them from the page cache.
        match &self.backing {
            RefBacking::File(m) => unsafe {
                let _ = m.unchecked_advise(memmap2::UncheckedAdvice::DontNeed);
            },
            RefBacking::Bundle(bytes, offset, len) => {
                if let Some(m) = bytes.mapping() {
                    // Whole pages inside the section only: the index shares the mapping.
                    let start = (offset + 4095) / 4096 * 4096;
                    let end = (offset + len) / 4096 * 4096;
                    if end > start {
                        unsafe {
                            let _ = m.unchecked_advise_range(memmap2::UncheckedAdvice::DontNeed, start, end - start);
                        }
                    }
                }
            }
        }
    }

    /// Writes the reference section: `MAGIC(8) | n(8) | offsets((n+1)*8) | data`. Little-endian;
    /// `offsets[0] == 0`. Returns the bytes written.
    fn write_to<'s, W: std::io::Write>(f: &mut W, n: usize, seq: &dyn Fn(usize) -> &'s [u8]) -> std::io::Result<u64> {
        f.write_all(&REFSEQ_MAGIC)?;
        f.write_all(&(n as u64).to_le_bytes())?;
        let mut off: u64 = 0;
        f.write_all(&off.to_le_bytes())?;
        for i in 0..n {
            off += seq(i).len() as u64;
            f.write_all(&off.to_le_bytes())?;
        }
        for i in 0..n {
            f.write_all(seq(i))?;
        }
        Ok(16 + 8 * (n as u64 + 1) + off)
    }

    /// The reference section as its own file (the layout of [`write_to`](Self::write_to)).
    fn write(filename: &str, records: &[OwnedFastaRecord]) -> std::io::Result<()> {
        let mut f = BufWriter::new(File::create(filename)?);
        Self::write_to(&mut f, records.len(), &|i| records[i].seq())?;
        std::io::Write::flush(&mut f)
    }

    /// Pointers into `bytes`, which hold a reference section; `None` if they do not.
    fn parse(bytes: &[u8]) -> Option<(*const u64, *const u8, usize)> {
        if bytes.len() < 16 || bytes[0..8] != REFSEQ_MAGIC {
            return None;
        }
        let n = u64::from_le_bytes(bytes[8..16].try_into().ok()?) as usize;
        let offsets_start = 16usize;
        let data_start = offsets_start + (n + 1) * 8;
        if bytes.len() < data_start {
            return None;
        }
        let base = bytes.as_ptr();
        // offsets_start (16) is 8-byte aligned and the section starts on a page, so the u64 reads
        // are aligned.
        Some((unsafe { base.add(offsets_start) } as *const u64, unsafe { base.add(data_start) }, n))
    }

    fn mmap_from_file(filename: &str) -> Option<Self> {
        let file = File::open(filename).ok()?;
        let mmap = unsafe { memmap2::Mmap::map(&file).ok()? };
        let (offsets, data, n) = Self::parse(&mmap)?;
        Some(Self { backing: RefBacking::File(std::sync::Arc::new(mmap)), offsets, data, n })
    }

    fn from_bundle(bundle: &Bundle) -> Option<Self> {
        let (offset, len) = bundle.section(bundle::REFSEQ)?;
        let (offsets, data, n) = Self::parse(&bundle.bytes.bytes()[offset..offset + len])?;
        Some(Self { backing: RefBacking::Bundle(bundle.bytes.clone(), offset, len), offsets, data, n })
    }
}

/// Where the reference sequences live: freshly parsed owned records (build, or a load before the
/// blob exists) or a memory-mapped [`RefSeqBlob`]. Both answer `get`/`len`.
#[derive(Clone)]
pub enum RefStore {
    Owned(Vec<OwnedFastaRecord>),
    Mmap(RefSeqBlob),
    /// 2 bits a base, decoded per read into this thread's scratch (`refseq2`).
    Packed(RefSeq2),
}
impl RefStore {
    /// Sequence `id` for alignment: a long packed sequence (a chromosome) is decoded a window at
    /// a time (`RefView`), a short one whole, as `get` does.
    #[inline]
    fn view(&self, id: usize) -> crate::database::refview::RefView<'_> {
        use crate::database::refview::RefView;
        match self {
            RefStore::Packed(p) if p.seq_len(id) > refseq2::MAX_PACKED_LEN => {
                RefView::Packed { store: p, id, len: p.seq_len(id) }
            }
            _ => RefView::Bytes(self.get(id)),
        }
    }

    /// Sequence `id` whole (a packed one decoded whole: for writing it out, not per read).
    #[inline]
    fn get(&self, id: usize) -> &[u8] {
        match self {
            RefStore::Owned(v) => v[id].seq(),
            RefStore::Mmap(b) => b.get(id),
            RefStore::Packed(p) => p.get(id),
        }
    }

    /// The length of sequence `id`; a packed one is not decoded.
    #[inline]
    fn seq_len(&self, id: usize) -> usize {
        match self {
            RefStore::Packed(p) => p.seq_len(id),
            _ => self.get(id).len(),
        }
    }

    /// The bundle's reference: the 2-bit section when it has one, else the bytes.
    fn from_bundle(bundle: &Bundle) -> Option<Self> {
        match bundle.section(bundle::REFSEQ2) {
            Some((offset, len)) => RefSeq2::from_shared(&bundle.bytes, offset, len).map(RefStore::Packed),
            None => RefSeqBlob::from_bundle(bundle).map(RefStore::Mmap),
        }
    }
    fn release(&self) {
        if let RefStore::Mmap(b) = self {
            b.release();
        }
    }
    fn len(&self) -> usize {
        match self {
            RefStore::Owned(v) => v.len(),
            RefStore::Mmap(b) => b.len(),
            RefStore::Packed(p) => p.len(),
        }
    }
}


#[repr(C)]
#[derive(Clone)]
pub struct DB<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> {
    flexmap: FlexmapStore<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
    rid_to_rname: Names,
    rname_to_rid: HashMap<String, usize>,
    references: RefStore,
}

impl<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> {
    /// Whether an existing on-disk index blob can be loaded by this build with these const params.
    /// `false` (e.g. after a flexmap key-format change bumped the blob version, or a C/F/params
    /// change) means the index must be rebuilt rather than loaded silently against a mismatched
    /// format -- callers OR this into their build-vs-load decision. Cheap: reads only the header.
    pub fn index_compatible(paths: &DBPaths) -> bool {
        let blob = paths.valid_paths()
            && FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::is_compatible(&paths.index_blob_file());
        let names = paths.names_valid();
        match index_format() {
            IndexFormat::Blob => blob,
            // The bundle of this format; or a packed index of it, or a blob to convert to one, with
            // the separate name files.
            format => Self::bundle_compatible(paths, format)
                || (names && (blob || Self::packed_compatible(paths, format))),
        }
    }

    fn packed_compatible(paths: &DBPaths, format: IndexFormat) -> bool {
        FlexmapStore::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::packed_compatible(paths, format)
    }

    /// The index of `index_format()` (see [`FlexmapStore::load`]).
    fn load_index(paths: &DBPaths) -> (FlexmapStore<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, String) {
        FlexmapStore::load(paths, index_format())
    }

    /// Load ONLY the reference (sequences + names), not the k-mer index.
    ///
    /// The sharded pipeline's rejoin/align half replays the anchors the shard passes already found
    /// (see `shard::disk::run_disk_rejoin_align`: it `gather`s evidence and calls `align_from_seeds`,
    /// never a range extractor), so `get_vrange` is never invoked there. The index is therefore not
    /// needed at all in the join -- only the reference bases to extend/align against and the names to
    /// write. Reading the 30 GB index here, as `load` does, is pure waste.
    ///
    /// The index blob is memory-mapped so the `flexmap` field is a valid value, but it is never
    /// warmed and never queried, so **zero index bytes are read** -- only the reference is. Calling
    /// `get_vrange` on the result would fault the index in on demand; the join never does.
    pub fn load_reference_only(paths: &DBPaths) -> Self {
        // A packed index's bundle holds the reference and the names after the index: read those
        // sections alone.
        let bundle_file = paths.bundle_file(index_format());
        if std::path::Path::new(&bundle_file).exists() {
            let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
            let first = if Bundle::section_in_file(&bundle_file, bundle::REFSEQ2).is_ok() { bundle::REFSEQ2 } else { bundle::REFSEQ };
            let (d, bundle) = time(|| {
                Bundle::open_from(&bundle_file, first, threads).unwrap_or_else(|e| panic!("bundle {bundle_file}: {e}"))
            });
            let references = RefStore::from_bundle(&bundle).expect("bundle has the reference");
            let names = NameTable::from_bundle(&bundle).expect("bundle has the names");
            eprintln!(
                "Reference-only load (join phase): {:.2} GB of {bundle_file} (reference and names) in {:.2}s -- index NOT read",
                bundle.bytes.len() as f64 / 1e9,
                d.as_secs_f64()
            );
            let rname_to_rid: HashMap<String, usize> = if crate::GOLDSTD_EVAL {
                (1..names.len()).map(|i| (names.get(i).unwrap().to_string(), i)).collect()
            } else {
                HashMap::new()
            };
            return Self { flexmap: FlexmapStore::Absent, rid_to_rname: Names::Table(names), rname_to_rid, references };
        }
        let io_before = super::common::disk_read_bytes();

        // mmap only: opens the file into address space, reads no data pages. A packed index has no
        // blob, and the join never queries one anyway.
        let blob = paths.index_blob_file();
        let flexmap = if std::path::Path::new(&blob).exists() {
            FlexmapStore::Mmap(FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::mmap_from_file(&blob))
        } else {
            FlexmapStore::Absent
        };

        // Names are part of the reference, needed to label output records.
        let (d_names, rid_to_rname) = time(|| -> Vec<String> {
            load(&mut buffered(&paths.id2reference_path, "id2ref"), crate::GLOBAL_VERSION)
                .expect("Valid id2ref")
        });

        let mapped = RefSeqBlob::mmap_from_file(&paths.refseq_file());
        let references = match mapped {
            Some(blob) => {
                // The reference is read UP FRONT, as everywhere else: eager is the default and
                // `--lazy-ref` is the explicit opt-out. The join needs the reference and nothing
                // else, so loading it here is loading exactly what this phase uses.
                let (_d, ref_bytes) = time(|| {
                    if super::common::eager_ref() { super::common::warm_file(&paths.refseq_file()) } else { 0 }
                });
                let gb = ref_bytes as f64 / (1u64 << 30) as f64;
                eprintln!(
                    "Reference-only load (join phase): names {:.2}s ({} entries), reference {} -- \
                     index NOT read",
                    d_names.as_secs_f64(), rid_to_rname.len(),
                    if super::common::eager_ref() { format!("{gb:.1} GB, sequential") }
                    else { "--lazy-ref: paged in during alignment".to_string() },
                );
                RefStore::Mmap(blob)
            }
            None => panic!(
                "reference-sequence blob {} is missing; build the index once (which writes it) \
                 before a sharded run", paths.refseq_file()
            ),
        };

        let read_gb = (super::common::disk_read_bytes().saturating_sub(io_before)) as f64 / (1u64 << 30) as f64;
        eprintln!("  -> {read_gb:.1} GB read from disk (reference only; the {} index is mapped but untouched)",
                  paths.index_blob_file());

        Self { flexmap, rid_to_rname: Names::Owned(rid_to_rname), rname_to_rid: HashMap::new(), references }
    }
}

impl<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>  FlexalignDatabase for DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> {
    fn get_rid(&self, reference: &str) -> Option<&usize> {
        self.rname_to_rid.get(reference)
    }

    fn get_rname(&self, id: usize) -> Option<&str> {
        self.rid_to_rname.get(id)
    }

    fn get_vrange(&self, canonical_kmer: u64) -> Option<AnyRange<'_>> {
        self.flexmap.get_vrange(canonical_kmer)
    }

    fn for_each_vrange<'s>(&'s self, keys: &[u64], f: &mut dyn FnMut(usize, Option<AnyRange<'s>>)) {
        self.flexmap.for_each_vrange(keys, f)
    }

    fn prefetch_keys(&self, keys: &[u64]) {
        self.flexmap.prefetch_keys(keys)
    }

    fn release_references(&self) {
        self.references.release();
    }

    fn build(options: &crate::options::Options) -> Self {
        let db_paths = DBPaths::new(&options.reference);

        let result = flexmap::build::default_build::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(
            &options.reference, options.args.max_range_size
        );

        let (flexmap, rname_to_rid, rid_to_rname) = match result {
            Ok(result) => {
                (result.0, result.1, result.2)
            },
            Err(_) => todo!(),
        };

        let references_file = &mut File::open(&db_paths.reference_path).expect("Working references file");
        let references = load_references(references_file, &rname_to_rid, &rid_to_rname);

        let references = match references {
            Ok(references) => references,
            Err(why) => panic!("Could not load references {}", why),
        };
        // A packed format is converted here, so the run aligns against what it will save and load.
        let (d_pack, flexmap) = time(|| FlexmapStore::Owned(flexmap).into_packed(index_format()));
        if index_format() != IndexFormat::Blob {
            eprintln!(
                "Packed the index ({:?} keys) in {:.1}s: {:.2} GB",
                index_format(), d_pack.as_secs_f64(), flexmap.memory_bytes() as f64 / 1e9
            );
        }
        Self {
            flexmap,
            rid_to_rname: Names::Owned(rid_to_rname),
            rname_to_rid,
            references: RefStore::Owned(references),
        }
    }

    fn load(paths: &super::common::DBPaths, version: u32) -> Self {
        let _ = version;
        if index_format() != IndexFormat::Blob && Self::bundle_compatible(paths, index_format()) {
            return Self::load_bundle(paths, index_format());
        }
        let io_before = crate::database::common::disk_read_bytes();

        let rid2rname_file = &mut buffered(&paths.id2reference_path, "id2ref");
        let rname2rid_file = &mut buffered(&paths.reference2id_path, "ref2id");
        let references_file = &mut File::open(&paths.reference_path).expect("Working references file");

        // The blob is memory-mapped (open is near-instant, and it is read in full below); a packed
        // index is read into memory on all threads, or converted from the blob on first use.
        let (d_index, (flexmap, index_note)) = time(|| Self::load_index(paths));
        let index_note = format!("{index_note}, {:.2} GB", flexmap.memory_bytes() as f64 / 1e9);


        // let config = bincode::config::standard();
        // let flexmap = decode_from_reader(map_reader, config).expect("Valid reference database");


        // The reference-name table is DESERIALISED into the heap -- 14.5M owned Strings on the
        // protal marker DB -- and is the dominant cost of a load. It was previously untimed, which
        // made the whole load look like it cost only the two mmap calls above (microseconds).
        let (d_names, rid_to_rname) = time(|| -> Vec<String> {
            load(rid2rname_file, version).expect("Valid reference database")
        });

        // Reference sequences: `mmap` the pre-built blob when it exists (near-instant); otherwise
        // parse the FASTA once (the slow path) and write the blob so the next load is fast. The
        // rname->rid map is only needed to attribute FASTA records to ids and for gold-standard
        // evaluation, so on the fast (mmap) path we skip loading that 14.5M-entry map entirely
        // unless GOLDSTD_EVAL is compiled in.
        let (d_refmap, mapped) = time(|| RefSeqBlob::mmap_from_file(&paths.refseq_file()));
        // The INDEX is read in full, here, sequentially. Lazy paging is wrong for it: k-mer lookups
        // hit the whole extent in essentially random order, so deferring the transfer does not avoid
        // it -- it just converts one streaming read into millions of 4K faults, and hides the cost
        // inside alignment where it looks like compute.
        let (d_prefault, idx_bytes) = time(|| match &flexmap {
            FlexmapStore::Mmap(_) => {
                super::common::warm_file(&paths.index_blob_file())
            }
            _ => 0,
        });

        // The REFERENCE is read up front too, by default: "load the database" should mean the whole
        // database, and the same argument applies -- the bytes are read either way, and doing it in
        // one sequential pass lets the kernel read ahead and keeps the cost inside the load figure
        // instead of smeared through alignment. --lazy-ref opts out for runs bounded by RAM rather
        // than time, where most of the reference is genuinely never touched.
        let (d_refwarm, ref_bytes) = time(|| match (&mapped, super::common::eager_ref()) {
            (Some(_), true) => super::common::warm_file(&paths.refseq_file()),
            _ => 0,
        });

        let mut d_refparse = std::time::Duration::ZERO;
        // Assigned in every branch below, so no initialiser is needed.
        let d_rev;
        let mut ref_mode = "mmap of blob";
        let (references, rname_to_rid) = match mapped {
            Some(blob) => {
                let (d, rname_to_rid) = time(|| {
                    if crate::GOLDSTD_EVAL {
                        load(rname2rid_file, version).expect("Valid reference database")
                    } else {
                        HashMap::new()
                    }
                });
                d_rev = d;
                (RefStore::Mmap(blob), rname_to_rid)
            }
            None => {
                ref_mode = "parsed from FASTA";
                let (d, rname_to_rid) = time(|| -> HashMap<String, usize> {
                    load(rname2rid_file, version).expect("Valid reference database")
                });
                d_rev = d;
                let (duration, parsed) =
                    time(|| load_references(references_file, &rname_to_rid, &rid_to_rname));
                d_refparse = duration;
                let parsed = match parsed {
                    Ok(r) => r,
                    Err(why) => panic!("Could not load references {}", why),
                };
                // Persist the blob so subsequent loads mmap it instead of re-parsing the FASTA.
                if let Err(e) = RefSeqBlob::write(&paths.refseq_file(), &parsed) {
                    log::warn!("Could not write reference-sequence blob ({e}); will re-parse next run.");
                }
                (RefStore::Owned(parsed), rname_to_rid)
            }
        };

        // One unambiguous breakdown, so nobody has to guess which step a load's seconds went to.
        // Only the DESERIALISED parts occupy heap here; the mmapped blobs are address space until
        // alignment touches them, and those page faults are charged to read processing, not to this.
        let secs = |d: std::time::Duration| d.as_secs_f64();
        let total = secs(d_index) + secs(d_names) + secs(d_refmap) + secs(d_refparse) + secs(d_rev)
            + secs(d_prefault) + secs(d_refwarm);
        let gb = |b: u64| b as f64 / (1u64 << 30) as f64;
        let read_gb = gb(crate::database::common::disk_read_bytes().saturating_sub(io_before));
        eprintln!(
            "Database load {total:.2}s total\n               index                 {:>8.3}s  ({})\n               reference names       {:>8.3}s  (deserialised, {} entries)\n               reference sequences   {:>8.3}s  ({})\n               name->id map          {:>8.3}s  ({})\n  \
             reading index in      {:>8.3}s  ({:.1} GB, sequential)\n  \
             reading reference in  {:>8.3}s  ({})\n  \
             -> {:.1} GB actually read from disk during load",
            secs(d_index), index_note,
            secs(d_names), rid_to_rname.len(),
            secs(d_refmap) + secs(d_refparse), ref_mode,
            secs(d_rev),
            if rname_to_rid.is_empty() { "skipped" } else { "deserialised" },
            secs(d_prefault), gb(idx_bytes),
            secs(d_refwarm),
            if super::common::eager_ref() {
                format!("{:.1} GB, sequential", gb(ref_bytes))
            } else {
                "--lazy-ref: paged in during alignment instead".to_string()
            },
            read_gb,
        );

        Self {
            flexmap,
            rid_to_rname: Names::Owned(rid_to_rname),
            rname_to_rid,
            references,
        }
    }
    
    fn save(&self, paths: &DBPaths, version: u32) -> Result<(), std::io::Error> {
        if matches!(self.flexmap, FlexmapStore::Compact(_) | FlexmapStore::Hashed(..)) {
            return self.save_bundle(paths);
        }
        self.save_separate(paths, version)
    }

    fn get_reference(&self, id: usize) -> Option<crate::database::refview::RefView<'_>> {
        if id < self.references.len() {
            Some(self.references.view(id))
        } else {
            None
        }
    }

    fn get_reference_len(&self, id: usize) -> Option<usize> {
        (id < self.references.len()).then(|| self.references.seq_len(id))
    }
}

impl<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> {
    /// The index, reference and names as separate files, as before bundles (and for the blob
    /// format, and for shard slicing, which reads the name files and the reference blob).
    pub fn save_separate(&self, paths: &DBPaths, version: u32) -> Result<(), std::io::Error> {
        let _ = version;
        match &self.flexmap {
            // A freshly built map: write the single-file blob (mmap-loaded next run) and, for the
            // shard slicer which still reads the split key array, the two-file form as well.
            FlexmapStore::Owned(f) => {
                f.save_blob(&paths.index_blob_file());
                f.save(&paths.index_keys_file(), &paths.index_values_file());
            }
            // Already a blob (e.g. re-saving a loaded index): copy its bytes back out.
            FlexmapStore::Mmap(b) => {
                std::fs::write(paths.index_blob_file(), b.backing_bytes())?;
            }
            // A packed index is one file, read in full at load; no blob is written for it.
            packed @ (FlexmapStore::Compact(_) | FlexmapStore::Hashed(..)) => {
                packed.save_packed(&paths.index_packed_file(index_format()))?;
            }
            FlexmapStore::Absent => {}
        }

        // The memory-mappable reference-sequence blob, so the next load skips re-parsing the FASTA.
        // A blob loaded from that file is already on disk; owned records and a bundle's section
        // are written out.
        if !matches!(&self.references, RefStore::Mmap(RefSeqBlob { backing: RefBacking::File(_), .. })) {
            let mut f = BufWriter::new(File::create(paths.refseq_file())?);
            RefSeqBlob::write_to(&mut f, self.references.len(), &|i| {
                refseq2::reset_scratch();
                self.references.get(i)
            })?;
            std::io::Write::flush(&mut f)?;
        }

        let names: Vec<String> = match &self.rid_to_rname {
            Names::Owned(v) => v.clone(),
            Names::Table(t) => (0..t.len()).map(|i| t.get(i).unwrap().to_string()).collect(),
        };
        let mut file = BufWriter::new(File::create(&paths.id2reference_path)
            .unwrap_or_else(|why| panic!("couldn't open {}: {}", paths.id2reference_path.display(), why)));
        let _ = save(&mut file, version, &names);
        let _ = std::io::Write::flush(&mut file);

        // A bundle load keeps no name -> id map; it is the inverse of the names.
        let rname_to_rid: HashMap<String, usize> = if self.rname_to_rid.is_empty() {
            names.iter().enumerate().skip(1).map(|(i, n)| (n.clone(), i)).collect()
        } else {
            self.rname_to_rid.clone()
        };
        let mut file = BufWriter::new(File::create(&paths.reference2id_path)
            .unwrap_or_else(|why| panic!("couldn't open {}: {}", paths.reference2id_path.display(), why)));
        let _ = save(&mut file, version, &rname_to_rid);
        let _ = std::io::Write::flush(&mut file);

        Ok(())
    }

    /// Writes the packed index, the reference and the names as one bundle
    /// (`paths.bundle_file(format)`).
    fn save_bundle(&self, paths: &DBPaths) -> Result<(), std::io::Error> {
        let file = paths.bundle_file(index_format());
        // The reference is stored at 2 bits a base (`refseq2::packable`): a marker database's
        // sequences are handed out whole, a genome's a window at a time (`RefView`). Decoding a
        // packed reference to write it again goes through this thread's scratch, reset per
        // sequence.
        let packed = refseq2::packable(self.references.len(), &|i| {
            refseq2::reset_scratch();
            self.references.get(i)
        });
        refseq2::reset_scratch();
        bundle::write(
            &file,
            |w| match &self.flexmap {
                FlexmapStore::Compact(p) => p.write_to(w),
                FlexmapStore::Hashed(p, _) => p.write_to(w),
                _ => unreachable!("a bundle holds a packed index"),
            },
            |w| {
                // Each sequence is used and dropped before the next is asked for, so a packed
                // store's scratch can be reset per sequence.
                let seq = |i: usize| {
                    refseq2::reset_scratch();
                    self.references.get(i)
                };
                if packed {
                    refseq2::write_to(w, self.references.len(), &seq)
                } else {
                    RefSeqBlob::write_to(w, self.references.len(), &seq)
                }
            },
            &|i| self.rid_to_rname.get(i).unwrap_or(""),
            self.rid_to_rname.len(),
            if packed { bundle::REFSEQ2 } else { bundle::REFSEQ },
        )
    }

    /// Whether `paths` has a bundle of `format` whose index this build reads.
    fn bundle_compatible(paths: &DBPaths, format: IndexFormat) -> bool {
        match bundle::Bundle::peek_index(&paths.bundle_file(format)) {
            Some(head) => match format {
                IndexFormat::Hashed => PackedFlexmap::<HashedKeys, C, F, HEADER_THRESHOLD>::is_compatible_header(&head),
                _ => PackedFlexmap::<CompactKeys, C, F, HEADER_THRESHOLD>::is_compatible_header(&head),
            },
            None => false,
        }
    }

    /// The whole database from its bundle: one parallel read into huge pages (or one mapping with
    /// `--index-mmap`), of which the index, the reference and the names are views.
    fn load_bundle(paths: &DBPaths, format: IndexFormat) -> Self {
        let file = paths.bundle_file(format);
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        let io_before = crate::database::common::disk_read_bytes();
        let (d_read, bundle) = time(|| Bundle::open(&file, index_mmap(), threads).unwrap_or_else(|e| panic!("bundle {file}: {e}")));
        let mmap = bundle.mapped;
        let (d_views, (flexmap, references, names)) = time(|| {
            let (offset, _) = bundle.section(bundle::INDEX).expect("bundle has an index");
            let flexmap = match format {
                IndexFormat::Hashed => FlexmapStore::hashed(PackedFlexmap::from_shared(&bundle.bytes, offset)),
                _ => FlexmapStore::Compact(Arc::new(PackedFlexmap::from_shared(&bundle.bytes, offset))),
            };
            let references = RefStore::from_bundle(&bundle).expect("bundle has the reference");
            let names = NameTable::from_bundle(&bundle).expect("bundle has the names");
            (flexmap, references, names)
        });
        // The name -> id map only serves the gold-standard evaluation.
        let rname_to_rid: HashMap<String, usize> = if crate::GOLDSTD_EVAL {
            (1..names.len()).map(|i| (names.get(i).unwrap().to_string(), i)).collect()
        } else {
            HashMap::new()
        };
        let gb = bundle.bytes.len() as f64 / 1e9;
        let read_gb = crate::database::common::disk_read_bytes().saturating_sub(io_before) as f64 / 1e9;
        eprintln!(
            "Database load {:.2}s total: bundle {file} ({gb:.2} GB)\n  \
             {:<10} {:>8.3}s  ({})\n  \
             {:<10} {:>8.3}s  (index {:.2} GB, {} reference sequences, {} names: views, nothing copied)\n  \
             -> {read_gb:.2} GB read from disk during load",
            d_read.as_secs_f64() + d_views.as_secs_f64(),
            if mmap { "mapped" } else { "read" },
            d_read.as_secs_f64(),
            if mmap {
                format!("mmap, populated on {threads} threads")
            } else if index_mmap() {
                format!("{threads} threads, huge pages; network filesystem, so read, not mapped")
            } else {
                format!("{threads} threads, huge pages")
            },
            "views",
            d_views.as_secs_f64(),
            flexmap.memory_bytes() as f64 / 1e9,
            references.len(),
            names.len(),
        );
        Self { flexmap, rid_to_rname: Names::Table(names), rname_to_rid, references }
    }
}


/// `--index-mmap`, set once at start-up (the loader has no options in reach).
pub static INDEX_MMAP: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn index_mmap() -> bool {
    INDEX_MMAP.load(std::sync::atomic::Ordering::Relaxed)
}

/// A name table opened for savefile, buffered: savefile reads each entry with its own small reads,
/// so a bare `File` costs two syscalls a name -- 5.4 M for 2.7 M names, seconds on Lustre.
fn buffered(path: &std::path::Path, what: &str) -> BufReader<File> {
    BufReader::with_capacity(1 << 20, File::open(path).unwrap_or_else(|e| panic!("{what} file {}: {e}", path.display())))
}
