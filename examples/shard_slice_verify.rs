//! End-to-end check of the physical index slice (ProjectShard.md §5, §13.6).
//!
//! Slices a built index into N physical shard files, loads each as a `PhysicalShardDB` (keys slice +
//! shared values), and confirms the disk-loaded shards reproduce the unsharded pipeline's seeds --
//! i.e. the rebased `get_vrange` is correct and the RAM-bounded shard workflow works.
//!
//! Run:
//!   cargo +nightly run --release --example shard_slice_verify -- \
//!       <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards>

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

use flate2::read::GzDecoder;

use bioreader::utils::is_gzip;
use flexalign::database::common::{DBPaths, FlexalignDatabase};
use flexalign::database::flexmap::DB;
use flexalign::shard::acceptance::compare_seeds_with_shards;
use flexalign::shard::slice::{slice_index, PhysicalShardDB, SliceManifest};

const K: usize = 31;
const C: usize = 15;
const F: usize = 16;
const S: usize = 4;
const L: usize = 12;
const CELLS_PER_BODY: u64 = 16;
const HEADER_THRESHOLD: usize = 2;

const BUFFER_SIZE: usize = 1 << 24;
const THREADS: u32 = 4;
const MAX_BEST_FLEX: usize = 16;
const MAX_RANGE_SIZE: usize = 256;
const MAX_RANGES: usize = 15;
const MAX_HEADERED: usize = 15;

type Shard = PhysicalShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>;
type FullDb = DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>;

fn read_all(path: &str) -> Vec<u8> {
    let gz = is_gzip(path).unwrap_or(false);
    let file = File::open(path).unwrap_or_else(|e| panic!("open {}: {}", path, e));
    let mut buf = Vec::new();
    if gz {
        GzDecoder::new(file).read_to_end(&mut buf).expect("gunzip");
    } else {
        BufReader::new(file).read_to_end(&mut buf).expect("read");
    }
    buf
}

fn file_mib(path: &str) -> f64 {
    std::fs::metadata(path).map(|m| m.len() as f64 / (1 << 20) as f64).unwrap_or(0.0)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: {} <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards>", args[0]);
        std::process::exit(2);
    }
    let reference = &args[1];
    let n_shards: usize = args[4].parse().expect("n_shards must be an integer");
    let paths = DBPaths::new(PathBuf::from(reference));

    eprintln!("Slicing index into {n_shards} physical shards ...");
    let manifest: SliceManifest =
        slice_index::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(&paths, n_shards).expect("slice");

    let full_keys_mib = file_mib(&paths.index_keys_file());
    eprintln!("  full keys : {:.1} MiB", full_keys_mib);
    for (i, e) in manifest.shards.iter().enumerate() {
        eprintln!("  shard {i:>2} : blob {:>7.1} MiB  keyspace [{}, {})", file_mib(&e.blob_file), e.lo, e.hi);
    }
    eprintln!(
        "  -> each pass now loads ~1/{} of the {:.1} MiB keys array instead of all of it",
        n_shards, full_keys_mib
    );

    eprintln!("\nLoading full index (oracle) and physical shards ...");
    let full: FullDb = DB::load(&paths, 1);
    let shards: Vec<Shard> = (0..manifest.shards.len()).map(|i| Shard::load(&manifest, i)).collect();

    let reads1 = read_all(&args[2]);
    let reads2 = read_all(&args[3]);

    eprintln!("Comparing physical-shard rejoin seeds vs unsharded pipeline ...");
    let report = compare_seeds_with_shards::<K, C, F, S, L, Shard, FullDb>(
        &reads1, &reads2, &shards, &full, BUFFER_SIZE, THREADS, MAX_BEST_FLEX, MAX_RANGE_SIZE,
        0, // recovery off: any divergence is the budget-boundary/cap, not §4 recovery
        MAX_RANGES, MAX_HEADERED,
    );

    eprintln!("\n=== physical-shard acceptance ===");
    eprintln!("  reads with seeds : {}", report.reads_with_seeds);
    eprintln!("  identical        : {} ({:.4} agreement)", report.identical, report.agreement());
    eprintln!("  differing        : {}", report.differing);
    if !report.differing_examples.is_empty() {
        eprintln!("  differing ordinals: {:?}", report.differing_examples);
    }
    eprintln!("  seeds  unsharded={}  sharded={}", report.unsharded_seed_total, report.sharded_seed_total);
    eprintln!(
        "  shard emit: headered={} headerless={} discarded={} capped={} cells={}",
        report.emit.headered, report.emit.headerless, report.emit.discarded, report.emit.capped, report.emit.cells
    );
}
