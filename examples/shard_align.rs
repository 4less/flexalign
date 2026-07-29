//! End-to-end sharded alignment producing PAF, RAM-bounded (ProjectShard.md §5, §8, §13.10).
//!
//! Slices the index into N physical shards (one-time; reused if present), then for each shard
//! streams over *all* reads and writes its evidence to a per-shard chunk file on disk -- dropping
//! the shard before the next loads. The rejoin re-streams the reads and, per batch, reads just that
//! batch's evidence back from the shard files, rebuilds seeds, and aligns into PAF. Evidence is
//! never accumulated in RAM, so peak memory is ~one shard (mmapped) + one batch's evidence,
//! independent of the read count.
//!
//! Multi-sample batches put the SAMPLES INSIDE the shard loop, so each shard is loaded once for the
//! whole batch rather than once per sample -- the load is the dominant cost on a big reference, and
//! N separate invocations pay it N times over. The full index for the align half is likewise loaded
//! once and reused for every sample's rejoin.
//!
//! Run:
//!   cargo +nightly run --release --example shard_align -- \
//!       <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards> <out.paf>
//!
//!   # batch: one line per sample, `reads_1<TAB>reads_2<TAB>out` (blank/# lines ignored)
//!   cargo +nightly run --release --example shard_align -- \
//!       <reference.fna> <n_shards> --batch <samples.tsv>

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use clap::Parser;
use flate2::read::GzDecoder;

use bioreader::utils::is_gzip;
use flexalign::database::common::{DBPaths, FlexalignDatabase};
use flexalign::database::flexmap::DB;
use flexalign::io::output_buffer::{OutputBuffer, OutputTarget};
use flexalign::options::{Args, Options};
use flexalign::shard::disk::{count_batches, run_disk_rejoin_align, run_shard_pass_to_disk};
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

/// One sample of a batch: its read files, its output, and its reads once loaded.
struct Sample {
    r1: String,
    r2: String,
    out: String,
    reads1: Vec<u8>,
    reads2: Vec<u8>,
    n_batches: u64,
}

/// `reads_1<TAB>reads_2<TAB>out` per line; blank lines and `#` comments ignored.
fn parse_batch(path: &str) -> Vec<(String, String, String)> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("open batch manifest {}: {}", path, e));
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            panic!("batch manifest {}:{}: expected reads_1<TAB>reads_2<TAB>out", path, n + 1);
        }
        out.push((f[0].to_string(), f[1].to_string(), f[2].to_string()));
    }
    if out.is_empty() {
        panic!("batch manifest {} has no samples", path);
    }
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let batch_at = a.iter().position(|x| x == "--batch");
    let usage = || -> ! {
        eprintln!("usage: {} <reference.fna> <reads_1.fq[.gz]> <reads_2.fq[.gz]> <n_shards> <out.paf>", a[0]);
        eprintln!("       {} <reference.fna> <n_shards> --batch <samples.tsv>", a[0]);
        std::process::exit(2)
    };

    // Two accepted forms; the single-sample positional one is unchanged so existing callers keep
    // working, and the batch form is what makes N samples cost one shard load instead of N.
    let (reference, n_shards, pairs): (&String, usize, Vec<(String, String, String)>) =
        match batch_at {
            Some(i) => {
                if a.len() < 4 || i + 1 >= a.len() {
                    usage();
                }
                (&a[1], a[2].parse().expect("n_shards"), parse_batch(&a[i + 1]))
            }
            None => {
                if a.len() < 6 {
                    usage();
                }
                (&a[1], a[4].parse().expect("n_shards"),
                 vec![(a[2].clone(), a[3].clone(), a[5].clone())])
            }
        };

    // Emit SAM when the output path ends in `.sam` (PAF otherwise). The align half already supports
    // both via `make_output`; SAM is the RAM-bounded design paired with the important CIGAR output.
    // One batch cannot mix the two: the format is a property of the run, not of a single output.
    let sam = pairs[0].2.ends_with(".sam");
    if pairs.iter().any(|p| p.2.ends_with(".sam") != sam) {
        panic!("batch mixes .sam and .paf outputs; use one format per run");
    }
    let mut argv: Vec<&str> = vec!["flexalign", "-r", reference, "-1", &pairs[0].0, "-2", &pairs[0].1];
    if sam {
        argv.push("--sam");
    }
    // The align half is what maps the 15 GB refseq blob, so --lazy-ref matters most here.
    if let Some(i) = a.iter().position(|x| x == "--ref-budget") {
        argv.push("--ref-budget");
        argv.push(&a[i + 1]);
    } else if a.iter().any(|x| x == "--lazy-ref") {
        argv.push("--lazy-ref");
    }
    let args = Args::parse_from(argv);
    let options = Options::from_args(args);
    {
        let spec = options.args.ref_budget.clone()
            .or_else(|| options.args.lazy_ref.then(|| "auto".to_string()));
        let kb = flexalign::database::common::resolve_ref_budget(spec.as_deref());
        flexalign::database::common::REF_BUDGET_KB.store(kb, std::sync::atomic::Ordering::Relaxed);
        if kb > 0 {
            eprintln!("--ref-budget: capping resident reference at {:.1} GB", kb as f64 / 1048576.0);
        }
    }
    let paths = DBPaths::new(PathBuf::from(reference));

    // One-time slice (reused if present) -- the only step that loads the full index.
    let manifest = match SliceManifest::load(&SliceManifest::path_for(&paths, n_shards)) {
        Ok(m) if m.shards.len() == n_shards => {
            eprintln!("Reusing existing slice ({n_shards} shards)");
            m
        }
        _ => {
            eprintln!("Slicing index into {n_shards} physical shards ...");
            slice_index::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(&paths, n_shards).expect("slice")
        }
    };

    // All samples' reads are resident at once: the shard loop is the OUTER one, so every sample has
    // to be streamable while a given shard is loaded. That is the whole point -- a shard costs one
    // load for the batch instead of one per sample -- and reads are small next to a shard.
    let samples: Vec<Sample> = pairs
        .into_iter()
        .map(|(r1, r2, out)| {
            let reads1 = read_all(&r1);
            let reads2 = read_all(&r2);
            let n_batches = count_batches(reads1.as_slice(), reads2.as_slice(), BUFFER_SIZE, THREADS);
            Sample { r1, r2, out, reads1, reads2, n_batches }
        })
        .collect();
    eprintln!(
        "{} sample(s), batches: {}",
        samples.len(),
        samples.iter().map(|s| s.n_batches.to_string()).collect::<Vec<_>>().join(", ")
    );

    // Evidence is per (shard, sample) -- keyed on the sample's own output path so concurrent runs
    // and re-runs cannot collide.
    let evidence: Vec<Vec<String>> = (0..manifest.shards.len())
        .map(|i| samples.iter().map(|s| format!("{}.evidence.s{i}.bin", s.out)).collect())
        .collect();

    // Shard passes -> disk. One shard resident at a time, now amortised over every sample.
    let t_pass = Instant::now();
    for i in 0..manifest.shards.len() {
        let shard: Shard = Shard::load(&manifest, i);
        for (j, s) in samples.iter().enumerate() {
            let emit = run_shard_pass_to_disk::<K, C, F, S, L, Shard, &[u8]>(
                &s.reads1, &s.reads2, &shard, s.n_batches, &evidence[i][j], BUFFER_SIZE, THREADS,
                options.args.max_best_flex, MAX_HEADERED,
            )
            .expect("shard pass");
            eprintln!(
                "  shard {i} / sample {j} ({}): emit headered={} headerless={}",
                s.r1, emit.headered, emit.headerless
            );
        }
        // `shard` drops here (mmap unmaps), freeing it before the next shard loads.
    }
    let pass_secs = t_pass.elapsed().as_secs_f64();

    // Rejoin + align: streams evidence back per batch. The full index is loaded ONCE and reused for
    // every sample -- the second place a batch would otherwise repeat a multi-GB load.
    eprintln!("Loading full index for the align half ...");
    let full: FullDb = DB::load(&paths, 1);

    let t_align = Instant::now();
    let mut total_alignments = 0usize;
    for (j, s) in samples.iter().enumerate() {
        let sink = Arc::new(Mutex::new(OutputTarget::File(
            File::create(&s.out).unwrap_or_else(|e| panic!("create {}: {}", s.out, e)),
        )));
        let out_buffer = OutputBuffer::new(Arc::clone(&sink), BUFFER_SIZE);
        let ev: Vec<String> = (0..manifest.shards.len()).map(|i| evidence[i][j].clone()).collect();
        let stats = run_disk_rejoin_align::<K, C, F, S, L, FullDb, &[u8]>(
            &s.reads1, &s.reads2, &full, &options, out_buffer, &ev, options.args.ranges,
            BUFFER_SIZE, THREADS,
        )
        .expect("rejoin");
        total_alignments += stats.alignments_successful as usize;
        eprintln!("  sample {j}: alignments={} -> {}", stats.alignments_successful, s.out);
        for p in &ev {
            let _ = std::fs::remove_file(p);
        }
    }
    let align_secs = t_align.elapsed().as_secs_f64();

    eprintln!(
        "\nsharded done: {} sample(s), {} shard(s); passes {:.2}s + rejoin/align {:.2}s ; alignments={}",
        samples.len(),
        manifest.shards.len(),
        pass_secs,
        align_secs,
        total_alignments
    );
}
