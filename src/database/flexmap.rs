use std::{collections::HashMap, fs::File, io::BufWriter};

use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use flexmap::packed::{CompactKeys, HashedKeys, PackedFlexmap, SlotLocator};
use flexmap::range::AnyRange;
use std::sync::Arc;
use savefile::{load, save};

use crate::flexalign::time;

use super::hashed_lines::HashedLines;
use super::common::{DBPaths, load_references, FlexalignDatabase, IndexFormat, index_format};


/// The full-mode index. `Owned` is a freshly built map still on the heap (build-then-align in one
/// run); `Mmap` is the memory-mapped blob loaded from disk -- open is near-instant and only the
/// keys/values actually queried are paged in, so a cold multi-GB index no longer costs a full read.
/// `Compact` and `Hashed` are flexmap's packed index (`--index-format`, compact by default), read
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
    // `Arc` so the DB wrapper stays cheaply cloneable per worker while sharing one mapping; the raw
    // pointers remain valid views into it.
    backing: std::sync::Arc<memmap2::Mmap>,
    offsets: *const u64, // n+1 little-endian offsets into the data section; offsets[0] == 0
    data: *const u8,
    n: usize,
}
unsafe impl Send for RefSeqBlob {}
unsafe impl Sync for RefSeqBlob {}
impl Clone for RefSeqBlob {
    fn clone(&self) -> Self {
        Self { backing: std::sync::Arc::clone(&self.backing), offsets: self.offsets, data: self.data, n: self.n }
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
    ///
    /// There is deliberately no `prefault` counterpart here: the up-front reference read goes
    /// through `common::warm_file`, which reads the file with several threads (~2 GB/s) rather than
    /// touching mapped pages one at a time on a single thread (~1 GB/s). Same bytes, half the time.
    fn release(&self) {
        // `DontNeed` is `unchecked` in memmap2 because on a writable private mapping it discards
        // un-written changes. This mapping is read-only and file-backed (`Mmap`, never `MmapMut`),
        // so there is nothing to lose: the pages are clean, and a later `get` re-faults them from
        // the page cache.
        unsafe {
            let _ = self.backing.unchecked_advise(memmap2::UncheckedAdvice::DontNeed);
        }
    }

    /// Layout: `MAGIC(8) | n(8) | offsets((n+1)*8) | data`. Little-endian; `offsets[0] == 0`.
    fn write(filename: &str, records: &[OwnedFastaRecord]) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = BufWriter::new(File::create(filename)?);
        f.write_all(&REFSEQ_MAGIC)?;
        f.write_all(&(records.len() as u64).to_le_bytes())?;
        let mut off: u64 = 0;
        f.write_all(&off.to_le_bytes())?;
        for r in records {
            off += r.seq().len() as u64;
            f.write_all(&off.to_le_bytes())?;
        }
        for r in records {
            f.write_all(r.seq())?;
        }
        Ok(())
    }

    fn mmap_from_file(filename: &str) -> Option<Self> {
        let file = File::open(filename).ok()?;
        let mmap = unsafe { memmap2::Mmap::map(&file).ok()? };
        if mmap.len() < 16 || mmap[0..8] != REFSEQ_MAGIC {
            return None;
        }
        let n = u64::from_le_bytes(mmap[8..16].try_into().ok()?) as usize;
        let offsets_start = 16usize;
        let data_start = offsets_start + (n + 1) * 8;
        if mmap.len() < data_start {
            return None;
        }
        let base = mmap.as_ptr();
        // offsets_start (16) is 8-byte aligned and the mmap base is page-aligned, so the u64 reads
        // are aligned.
        let offsets = unsafe { base.add(offsets_start) } as *const u64;
        let data = unsafe { base.add(data_start) };
        Some(Self { backing: std::sync::Arc::new(mmap), offsets, data, n })
    }
}

/// Where the reference sequences live: freshly parsed owned records (build, or a load before the
/// blob exists) or a memory-mapped [`RefSeqBlob`]. Both answer `get`/`len`.
#[derive(Clone)]
pub enum RefStore {
    Owned(Vec<OwnedFastaRecord>),
    Mmap(RefSeqBlob),
}
impl RefStore {
    #[inline]
    fn get(&self, id: usize) -> &[u8] {
        match self {
            RefStore::Owned(v) => v[id].seq(),
            RefStore::Mmap(b) => b.get(id),
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
    rid_to_rname: Vec<String>,
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
            // A packed index of this format, or a blob to convert to one at load.
            IndexFormat::Compact => names && (blob || Self::packed_compatible(paths, IndexFormat::Compact)),
            IndexFormat::Hashed => names && (blob || Self::packed_compatible(paths, IndexFormat::Hashed)),
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
            load(&mut File::open(&paths.id2reference_path).expect("Working id2ref file"), crate::GLOBAL_VERSION)
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

        Self { flexmap, rid_to_rname, rname_to_rid: HashMap::new(), references }
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
        if (id as usize) < self.rid_to_rname.len() {
            return Some(&self.rid_to_rname[id as usize])
        } else {
            None
        }
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
            rid_to_rname,
            rname_to_rid,
            references: RefStore::Owned(references),
        }
    }

    fn load(paths: &super::common::DBPaths, version: u32) -> Self {
        let _ = version;
        let io_before = crate::database::common::disk_read_bytes();

        let rid2rname_file = &mut File::open(&paths.id2reference_path).expect("Working id2ref file");
        let rname2rid_file = &mut File::open(&paths.reference2id_path).expect("Working ref2id file");
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
            rid_to_rname,
            rname_to_rid,
            references,
        }
    }
    
    fn save(&self, paths: &DBPaths, version: u32) -> Result<(), std::io::Error> {
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

        // Write the memory-mappable reference-sequence blob so the next load skips re-parsing the
        // FASTA. Only possible when the sequences are still owned in RAM (a fresh build); a loaded
        // index is already backed by the blob.
        if let RefStore::Owned(records) = &self.references {
            RefSeqBlob::write(&paths.refseq_file(), records)?;
        }

        let mut file = match File::create(&paths.id2reference_path) {
            Err(why) => panic!(
                "couldn't open {}: {}",
                paths.id2reference_path.display(),
                why
            ),
            Ok(file) => file,
        };
        let _ = save(&mut file, version, &self.rid_to_rname);
    
    
        let mut file = match File::create(&paths.reference2id_path) {
            Err(why) => panic!(
                "couldn't open {}: {}",
                paths.id2reference_path.display(),
                why
            ),
            Ok(file) => file,
        };
        let _ = save(&mut file, version, &self.rname_to_rid);

        // let mut ser = PureCopySerializer::<16, 8, 16, 1024, _>::new();
        // let storage = ser.serialize(&self.flexmap);

        Ok(())
    }
    
    fn get_reference(&self, id: usize) -> Option<&[u8]> {
        if id < self.references.len() {
            Some(self.references.get(id))
        } else {
            None
        }
    }
}


/// `--index-mmap`, set once at start-up (the loader has no options in reach).
pub static INDEX_MMAP: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn index_mmap() -> bool {
    INDEX_MMAP.load(std::sync::atomic::Ordering::Relaxed)
}
