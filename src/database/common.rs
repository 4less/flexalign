use std::{collections::HashMap, io::Read, path::{Path, PathBuf}, sync::{Arc, Mutex}};
use bioreader::{fasta_byte_reader::FastaByteReader, fasta_reader::FastaReader, sequence::fasta_record::OwnedFastaRecord};
use flexmap::values::VRange;
use crate::options::Options;

const INDEX_EXTENSION: &str = ".flex.index";
const ID2REF_MAP_EXTENSION: &str = ".flex.id2ref";
const REF2ID_MAP_EXTENSION: &str = ".flex.ref2id";


pub struct DBPaths {
    pub reference_path: PathBuf,
    pub index_path: PathBuf,
    pub reference2id_path: PathBuf,
    pub id2reference_path: PathBuf,
}

impl DBPaths {
    pub fn new(reference_path: impl AsRef<Path>) -> Self {
        let reference_path = match reference_path.as_ref().exists() {
            true => reference_path,
            false => panic!("Reference file {} does not exist", reference_path.as_ref().display()),
        };

        let index_path = PathBuf::from(reference_path.as_ref().display().to_string() + INDEX_EXTENSION);
        let id2reference_path = PathBuf::from(reference_path.as_ref().display().to_string() + ID2REF_MAP_EXTENSION);
        let reference2id_path = PathBuf::from(reference_path.as_ref().display().to_string() + REF2ID_MAP_EXTENSION);
        
        DBPaths {
            reference_path: reference_path.as_ref().to_path_buf(),
            index_path,
            reference2id_path,
            id2reference_path,
        }
    }

    pub fn valid_paths(&self) -> bool {
        Path::exists(&self.reference_path) &
        Path::exists(&self.index_path) &
        Path::exists(&self.reference2id_path) &
        Path::exists(&self.id2reference_path)
    }
}

pub trait FlexalignDatabase {
    fn get_rid(&self, reference: &str) -> Option<&usize>;
    fn get_rname(&self, id: usize) -> Option<&str>;
    fn get_reference(&self, id: usize) -> Option<&[u8]>;
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange>;
    fn build(options: &Options) -> Self;
    fn save(&self, paths: &DBPaths, version: u32) -> Result<(), std::io::Error>;
    fn load(paths: &DBPaths, version: u32) -> Self;
}


pub fn load_references<R>(references_file: R, reference2id: &HashMap<String, usize>, id2reference: &Vec<String>) -> Result<Vec<OwnedFastaRecord>, std::io::Error> where R: Read {
    let buffer_size: usize = usize::pow(2, 24);
    let data = Mutex::new(FastaByteReader::new(references_file, buffer_size)?);
    let mut byte_reader = Arc::new(data);
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);

    let mut record = OwnedFastaRecord::new();

    let mut data = Vec::default();
    data.resize(id2reference.len(), OwnedFastaRecord::new());

    while let Some(()) = fasta_reader
        .load_batch_par(&mut byte_reader)
        .expect("Batch is invalid")
    {
        while let Some(_) = fasta_reader.next(&mut record) {
            if !record.valid_extended() {
                panic!("Record is not valid {:?}", record.to_string())
            }

            // let header = String::from_utf8_lossy(&record.head()[1..]).into_owned();
            let header = String::from_utf8_lossy(&record.head()[1..]).split(' ').next().unwrap().to_string();
            let reference_id = reference2id[&header];
            data[reference_id] = record.clone();
        }
    }
    Ok(data)
}