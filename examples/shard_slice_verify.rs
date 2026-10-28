//! End-to-end check of online shards (ProjectShard.md §5).
//!
//! Computes the bounds of N shards of a built bundle, reads each as a `WindowShardDB` (a window of
//! the bundle's index), and confirms the windows reproduce the unsharded pipeline's seeds -- i.e.
//! the window lookups answer as the whole index and the RAM-bounded shard workflow works.
//!
//! Run:
//!   cargo +nightly run --release --example shard_slice_verify -- \
//!       <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards> [hashed]

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

use flate2::read::GzDecoder;

use bioreader::utils::is_gzip;
use flexalign::database::common::{set_index_format, DBPaths, FlexalignDatabase, IndexFormat};
use flexalign::database::flexmap::DB;
use flexalign::shard::acceptance::compare_seeds_with_shards;
use flexalign::shard::slice::{ShardBounds, WindowShardDB};

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

type Shard = WindowShardDB<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>;
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: {} <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards> [hashed]", args[0]);
        std::process::exit(2);
    }
    let reference = &args[1];
    let n_shards: usize = args[4].parse().expect("n_shards must be an integer");
    let format = match args.get(5).map(String::as_str) {
        None | Some("hashed") => IndexFormat::Hashed,
        Some(other) => panic!("unknown index format {other}"),
    };
    // The oracle below loads the whole index in the format the shards are windows of.
    set_index_format(format);
    let paths = DBPaths::new(PathBuf::from(reference));

    let (bounds, _) = ShardBounds::load_or_compute::<C, F, HEADER_THRESHOLD>(&paths, format, n_shards, THREADS as usize)
        .expect("shard bounds");

    eprintln!("Loading full index (oracle) and window shards ...");
    let full: FullDb = DB::load(&paths, 1);
    let shards: Vec<Shard> = bounds
        .ranges()
        .into_iter()
        .map(|r| Shard::load(&paths, format, r, THREADS as usize).expect("window shard"))
        .collect();
    for (i, sh) in shards.iter().enumerate() {
        let r = sh.range();
        eprintln!("  shard {i:>2} : {:>7.1} MiB  keyspace [{}, {})", sh.memory_bytes() as f64 / (1 << 20) as f64, r.lo, r.hi);
    }

    let reads1 = read_all(&args[2]);
    let reads2 = read_all(&args[3]);

    eprintln!("Comparing window-shard rejoin seeds vs unsharded pipeline ...");
    let report = compare_seeds_with_shards::<K, C, F, S, L, Shard, FullDb>(
        &reads1, &reads2, &shards, &full, BUFFER_SIZE, THREADS, MAX_BEST_FLEX, MAX_RANGE_SIZE,
        0, // recovery off: any divergence is the budget-boundary/cap, not §4 recovery
        MAX_RANGES, MAX_HEADERED,
    );

    eprintln!("\n=== window-shard acceptance ===");
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
