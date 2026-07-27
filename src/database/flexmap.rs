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
    references: Vec<OwnedFastaRecord>,
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
            references,
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
        let rname_to_rid: HashMap<String, usize> = load(rname2rid_file, version).expect("Valid reference database");

        let (duration, references) = time(|| {
            load_references(references_file, &rname_to_rid, &rid_to_rname)
        });
        eprintln!("Loading references took {:?}", duration);

        let references = match references {
            Ok(references) => references,
            Err(why) => panic!("Could not load references {}", why),
        };

        Self {
            flexmap,
            rid_to_rname,
            rname_to_rid,
            references: references,
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
        Some(self.references[id].seq())
    }
}
