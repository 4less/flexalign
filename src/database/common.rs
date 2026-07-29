use std::{collections::HashMap, io::Read, path::{Path, PathBuf}, sync::{Arc, Mutex}};
use bioreader::{fasta_byte_reader::FastaByteReader, fasta_reader::FastaReader, sequence::fasta_record::OwnedFastaRecord};
use flexmap::values::VRange;
use crate::options::Options;

/// Budget in kB for resident reference pages; 0 = unbounded (the default).
///
/// A process-wide switch rather than a parameter because `FlexalignDatabase::load` is a trait
/// method taking only paths, and the mapping happens several layers below the option parsing.
/// Written once at startup, only read afterwards.
pub static REF_BUDGET_KB: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

pub fn ref_budget_kb() -> u64 {
    REF_BUDGET_KB.load(std::sync::atomic::Ordering::Relaxed)
}

/// Resolve `--ref-budget`: `auto`, a size like `8G`/`512M`, or absent (unbounded).
///
/// `auto` asks what this process is *actually* allowed to use, which is not the same as the
/// machine's RAM: under a cgroup (container, slurm cgroup, `systemd-run -p MemoryMax=`) the cap is
/// far below total memory, and that is exactly the case where being unbounded gets you OOM-killed.
/// So prefer the cgroup limit, fall back to MemAvailable, and take a fraction to leave room for the
/// non-reference parts of the run.
pub fn resolve_ref_budget(spec: Option<&str>) -> u64 {
    fn parse_size(s: &str) -> Option<u64> {
        let s = s.trim();
        let (num, mult) = match s.chars().last()?.to_ascii_uppercase() {
            'K' => (&s[..s.len() - 1], 1u64),
            'M' => (&s[..s.len() - 1], 1024),
            'G' => (&s[..s.len() - 1], 1024 * 1024),
            'T' => (&s[..s.len() - 1], 1024 * 1024 * 1024),
            _ => (s, 1024 * 1024), // bare number = GB, the unit anyone types for a reference budget
        };
        num.trim().parse::<f64>().ok().map(|v| (v * mult as f64) as u64)
    }

    fn cgroup_limit_kb() -> Option<u64> {
        let cg = std::fs::read_to_string("/proc/self/cgroup").ok()?;
        let path = cg.lines().find_map(|l| l.split(':').nth(2))?.to_string();
        let v = std::fs::read_to_string(format!("/sys/fs/cgroup{path}/memory.max")).ok()?;
        v.trim().parse::<u64>().ok().map(|b| b / 1024) // "max" fails to parse -> unbounded
    }

    fn mem_available_kb() -> Option<u64> {
        std::fs::read_to_string("/proc/meminfo").ok()?.lines()
            .find(|l| l.starts_with("MemAvailable:"))?
            .split_whitespace().nth(1)?.parse().ok()
    }

    match spec {
        None => 0,
        Some(s) if s.eq_ignore_ascii_case("auto") => {
            let avail = cgroup_limit_kb().or_else(mem_available_kb).unwrap_or(0);
            // 60%: the rejoin also holds anonymous memory (reads, evidence, output buffers), and a
            // budget that consumed everything available would just move the OOM somewhere else.
            (avail as f64 * 0.60) as u64
        }
        Some(s) => parse_size(s).unwrap_or(0),
    }
}

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

    /// The `Flexmap` index is stored as two files (keys + values); flexmap's
    /// `Flexmap::save`/`load` take them separately. Both are derived from
    /// `index_path` by suffix so the `.flex.index` naming stays the base.
    pub fn index_keys_file(&self) -> String {
        format!("{}.keys", self.index_path.to_string_lossy())
    }

    pub fn index_values_file(&self) -> String {
        format!("{}.values", self.index_path.to_string_lossy())
    }

    /// Single-file, memory-mappable index blob (keys + values). This is what full-mode loads.
    pub fn index_blob_file(&self) -> String {
        format!("{}.blob", self.index_path.to_string_lossy())
    }

    /// Memory-mappable reference-sequence blob (all sequences concatenated + an offset table), so a
    /// load can `mmap` it instead of re-parsing the multi-GB reference FASTA every run. Generated
    /// lazily on the first load that doesn't find it.
    pub fn refseq_file(&self) -> String {
        format!("{}.refseq", self.index_path.to_string_lossy())
    }

    pub fn valid_paths(&self) -> bool {
        Path::exists(&self.reference_path) &
        Path::exists(Path::new(&self.index_blob_file())) &
        Path::exists(Path::new(&self.index_keys_file())) &
        Path::exists(Path::new(&self.index_values_file())) &
        Path::exists(&self.reference2id_path) &
        Path::exists(&self.id2reference_path)
    }
}

pub trait FlexalignDatabase {
    fn get_rid(&self, reference: &str) -> Option<&usize>;
    fn get_rname(&self, id: usize) -> Option<&str>;
    fn get_reference(&self, id: usize) -> Option<&[u8]>;

    /// Drop resident pages of the reference-sequence mapping, if it is memory-mapped.
    ///
    /// Alignment touches a few hundred scattered bases per candidate, so over millions of
    /// candidates essentially the whole reference becomes resident and sets peak RSS -- above even
    /// a full shard. Calling this at a batch boundary bounds residency to one batch's worth
    /// instead. Default is a no-op (owned records cannot be dropped).
    fn release_references(&self) {}
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange<'_>>;
    fn build(options: &Options) -> Self;
    fn save(&self, paths: &DBPaths, version: u32) -> Result<(), std::io::Error>;
    fn load(paths: &DBPaths, version: u32) -> Self;
}

/// A shared reference to a database is itself a read-only database. This lets key-range views such
/// as `ShardedDB<&DB>` borrow one built index across many shards instead of cloning it -- essential
/// for sharding a multi-GB index without N copies of it in RAM. The mutating/persisting methods are
/// unreachable through a shared reference and panic; a shard pass only ever reads.
impl<D: FlexalignDatabase> FlexalignDatabase for &D {
    fn get_rid(&self, reference: &str) -> Option<&usize> {
        (**self).get_rid(reference)
    }
    fn get_rname(&self, id: usize) -> Option<&str> {
        (**self).get_rname(id)
    }
    fn get_reference(&self, id: usize) -> Option<&[u8]> {
        (**self).get_reference(id)
    }
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange<'_>> {
        (**self).get_vrange(canonical_kmer)
    }
    fn build(_options: &Options) -> Self {
        unimplemented!("&DB is a borrowed read-only view; build the owned database instead")
    }
    fn save(&self, _paths: &DBPaths, _version: u32) -> Result<(), std::io::Error> {
        unimplemented!("&DB is a borrowed read-only view")
    }
    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!("&DB is a borrowed read-only view; load the owned database instead")
    }
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