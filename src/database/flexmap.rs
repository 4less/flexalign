use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}};

use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use savefile::{load, save};

use crate::flexalign::time;

use super::common::{DBPaths, load_references, FlexalignDatabase};


/// The full-mode index. `Owned` is a freshly built map still on the heap (build-then-align in one
/// run); `Mmap` is the memory-mapped blob loaded from disk -- open is near-instant and only the
/// keys/values actually queried are paged in, so a cold multi-GB index no longer costs a full read.
/// Both answer `get_vrange`.
#[derive(Clone)]
pub enum FlexmapStore<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize> {
    Owned(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>),
    Mmap(FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>),
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexmapStore<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    #[inline]
    fn get_vrange(&self, canonical_kmer: u64) -> Option<flexmap::values::VRange<'_>> {
        match self {
            FlexmapStore::Owned(f) => f.get_vrange(canonical_kmer),
            FlexmapStore::Mmap(b) => b.get_vrange(canonical_kmer),
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
    /// is lost: a later access re-faults from the page cache. Used by `--lazy-ref` to bound the
    /// resident set to one batch's worth of reference instead of everything the run ever touched.
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
        FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::is_compatible(&paths.index_blob_file())
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

    fn get_vrange(&self, canonical_kmer: u64) -> Option<flexmap::values::VRange<'_>> {
        self.flexmap.get_vrange(canonical_kmer)
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
        Self {
            flexmap: FlexmapStore::Owned(flexmap),
            rid_to_rname,
            rname_to_rid,
            references: RefStore::Owned(references),
        }
    }

    fn load(paths: &super::common::DBPaths, version: u32) -> Self {
        let _ = version;

        let rid2rname_file = &mut File::open(&paths.id2reference_path).expect("Working id2ref file");
        let rname2rid_file = &mut File::open(&paths.reference2id_path).expect("Working ref2id file");
        let references_file = &mut File::open(&paths.reference_path).expect("Working references file");

        // Memory-map the single-file blob: open is near-instant and only the queried keys/values
        // are paged in, instead of reading the whole multi-GB index into the heap.
        let (duration, flexmap) = time(|| {
            FlexmapStore::Mmap(
                FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::mmap_from_file(&paths.index_blob_file()),
            )
        });
        eprintln!("Memory-mapping index took {:?}", duration);

        // let config = bincode::config::standard();
        // let flexmap = decode_from_reader(map_reader, config).expect("Valid reference database");


        let rid_to_rname: Vec<String> = load(rid2rname_file, version).expect("Valid reference database");

        // Reference sequences: `mmap` the pre-built blob when it exists (near-instant); otherwise
        // parse the FASTA once (the slow path) and write the blob so the next load is fast. The
        // rname->rid map is only needed to attribute FASTA records to ids and for gold-standard
        // evaluation, so on the fast (mmap) path we skip loading that 14.5M-entry map entirely
        // unless GOLDSTD_EVAL is compiled in.
        let (duration, mapped) = time(|| RefSeqBlob::mmap_from_file(&paths.refseq_file()));
        let (references, rname_to_rid) = match mapped {
            Some(blob) => {
                eprintln!("Memory-mapping references took {:?}", duration);
                let rname_to_rid = if crate::GOLDSTD_EVAL {
                    load(rname2rid_file, version).expect("Valid reference database")
                } else {
                    HashMap::new()
                };
                (RefStore::Mmap(blob), rname_to_rid)
            }
            None => {
                let rname_to_rid: HashMap<String, usize> =
                    load(rname2rid_file, version).expect("Valid reference database");
                let (duration, parsed) =
                    time(|| load_references(references_file, &rname_to_rid, &rid_to_rname));
                eprintln!("Loading references (parsing FASTA) took {:?}", duration);
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
