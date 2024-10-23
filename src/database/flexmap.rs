use std::{collections::HashMap, fs::File, io::BufReader};

use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::flexmap::{Flexmap, FlexmapHash, VRangeGetter};
use savefile::{load, save};
use ser_raw::{storage, CompleteSerializer, PureCopySerializer, Serialize, SerializeWith, Serializer};

use crate::flexalign::time;

use super::common::{DBPaths, load_references, FlexalignDatabase};


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
    flexmap: Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
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

    fn get_vrange(&self, canonical_kmer: u64) -> Option<flexmap::values::VRange> {
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
            flexmap,
            rid_to_rname,
            rname_to_rid,
            references,
        }
    }

    fn load(paths: &super::common::DBPaths, version: u32) -> Self {
        let map_file = &mut File::open(&paths.index_path).expect("Working flexmap file");
        let mut map_reader = BufReader::new(map_file);

        let rid2rname_file = &mut File::open(&paths.id2reference_path).expect("Working id2ref file");
        let rname2rid_file = &mut File::open(&paths.reference2id_path).expect("Working ref2id file");
        let references_file = &mut File::open(&paths.reference_path).expect("Working references file");

        let flexmap = load(&mut map_reader, version).expect("Valid reference database");

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
        let mut file = match File::create(&paths.index_path) {
            Err(why) => panic!("couldn't open {}: {}", paths.index_path.display(), why),
            Ok(file) => file,
        };
        let _ = save(&mut file, version, &self.flexmap);
    
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




#[derive(Clone)]
pub struct DBHash<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize,
> {
    flexmap: FlexmapHash<C, F, HEADER_THRESHOLD>,
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
    const HEADER_THRESHOLD: usize,
>  FlexalignDatabase for DBHash<K, C, F, S, L, HEADER_THRESHOLD> {
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

    fn get_vrange(&self, canonical_kmer: u64) -> Option<flexmap::values::VRange> {
        self.flexmap.get_vrange(canonical_kmer)
    }

    fn build(options: &crate::options::Options) -> Self {
        let db_paths = DBPaths::new(&options.reference);

        let result = flexmap::build::hash_build::<K, C, F, S, L, HEADER_THRESHOLD>(
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
            flexmap,
            rid_to_rname,
            rname_to_rid,
            references,
        }
    }

    fn load(paths: &super::common::DBPaths, version: u32) -> Self {
        let map_file = &mut File::open(&paths.index_path).expect("Working flexmap file");
        let mut map_reader = BufReader::new(map_file);

        let rid2rname_file = &mut File::open(&paths.id2reference_path).expect("Working id2ref file");
        let rname2rid_file = &mut File::open(&paths.reference2id_path).expect("Working ref2id file");
        let references_file = &mut File::open(&paths.reference_path).expect("Working references file");

        let flexmap = load(&mut map_reader, version).expect("Valid reference database");

        // let config = bincode::config::standard();
        // let flexmap = decode_from_reader(map_reader, config).expect("Valid reference database");


        let rid_to_rname: Vec<String> = load(rid2rname_file, version).expect("Valid reference database");
        let rname_to_rid: HashMap<String, usize> = load(rname2rid_file, version).expect("Valid reference database");

        let references = load_references(references_file, &rname_to_rid, &rid_to_rname);

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
        let mut file = match File::create(&paths.index_path) {
            Err(why) => panic!("couldn't open {}: {}", paths.index_path.display(), why),
            Ok(file) => file,
        };
        let _ = save(&mut file, version, &self.flexmap);
    
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

        Ok(())
    }
    
    fn get_reference(&self, id: usize) -> Option<&[u8]> {
        Some(self.references[id].seq())
    }
}