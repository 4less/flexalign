//! Seed-level acceptance of the shard cut against the unsharded pipeline, on a real index
//! (ProjectShard.md §13.9).
//!
//! Loads a built index, then for the same reads compares the seeds produced by
//!   (shard passes -> evidence -> rejoin)   vs.   (unsharded kmer -> range -> StdSeedExtractor).
//!
//! Run:
//!   cargo +nightly run --release --example shard_acceptance -- \
//!       <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards>
//!
//! The reference's `.flex.index.keys/.values` must already exist (build once via a normal run).

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

use flate2::read::GzDecoder;

use bioreader::utils::is_gzip;
use flexalign::database::common::{DBPaths, FlexalignDatabase};
use flexalign::database::flexmap::DB;
use flexalign::shard::acceptance::{compare_seeds, AcceptanceReport};

// The default profile's compile-time scheme (see flexalign::run): K=31 C=15 F=16 S=4 L=12.
const K: usize = 31;
const C: usize = 15;
const F: usize = 16;
const S: usize = 4;
const L: usize = 12;
const CELLS_PER_BODY: u64 = 16;
const HEADER_THRESHOLD: usize = 2;

const BUFFER_SIZE: usize = 1 << 24; // must match between the shard passes and the unsharded pass
const THREADS: u32 = 4;
const MAX_BEST_FLEX: usize = 16;
const MAX_RANGE_SIZE: usize = 256;
const MAX_RANGES: usize = 15;
const MAX_HEADERED: usize = 15;

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

fn print_report(label: &str, r: &AcceptanceReport) {
    eprintln!("\n=== {label} ===");
    eprintln!("  reads with seeds : {}", r.reads_with_seeds);
    eprintln!(
        "  identical        : {} ({:.4} agreement)",
        r.identical,
        r.agreement()
    );
    eprintln!("  differing        : {}", r.differing);
    if !r.differing_examples.is_empty() {
        eprintln!("  differing ordinals (first {}): {:?}", r.differing_examples.len(), r.differing_examples);
    }
    eprintln!(
        "  seeds  unsharded={}  sharded={}",
        r.unsharded_seed_total, r.sharded_seed_total
    );
    eprintln!(
        "  shard emit: headered={} headerless={} discarded={} capped={} cells={}",
        r.emit.headered, r.emit.headerless, r.emit.discarded, r.emit.capped, r.emit.cells
    );
}

fn run(reads1: &[u8], reads2: &[u8], db: &DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>, n_shards: usize, min_ranges: usize) -> AcceptanceReport {
    compare_seeds::<K, C, F, S, L, DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>>(
        reads1,
        reads2,
        db,
        n_shards,
        BUFFER_SIZE,
        THREADS,
        MAX_BEST_FLEX,
        MAX_RANGE_SIZE,
        min_ranges,
        MAX_RANGES,
        MAX_HEADERED,
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!(
            "usage: {} <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards>",
            args[0]
        );
        std::process::exit(2);
    }
    let reference = &args[1];
    let reads1_path = &args[2];
    let reads2_path = &args[3];
    let n_shards: usize = args[4].parse().expect("n_shards must be an integer");

    let paths = DBPaths::new(PathBuf::from(reference));
    eprintln!("Loading index for {reference} ...");
    let db: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = DB::load(&paths, 1);

    eprintln!("Reading {reads1_path} / {reads2_path} ...");
    let reads1 = read_all(reads1_path);
    let reads2 = read_all(reads2_path);

    eprintln!("\nComparing sharded ({n_shards} shards) vs unsharded seeds ...");

    // min_ranges = 0 disables the unsharded recovery path, so any divergence is purely the
    // budget-boundary tie-break (§3.3) rather than the recovery differences (§4).
    let exact = run(&reads1, &reads2, &db, n_shards, 0);
    print_report("min_ranges=0 (recovery off)", &exact);

    // The pipeline default. Extra divergence here is the §4-predicted recovery-path set.
    let with_recovery = run(&reads1, &reads2, &db, n_shards, 4);
    print_report("min_ranges=4 (pipeline default)", &with_recovery);
}
