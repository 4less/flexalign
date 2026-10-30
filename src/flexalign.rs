use std::time::{Duration, Instant};
use std::process::exit;

use crate::align::process_fastq::process_fastq_wrapper_modular;
use crate::database::flexmap::DB;
use crate::database::common::{DBPaths, FlexalignDatabase, IndexFormat};
use crate::options::{Args, Options};
use crate::GLOBAL_VERSION;


/// Whether the per-stage timers in the read pipeline are live. They wrap every stage of every
/// read, and two `Instant::now()` per stage was ~2% of a run, so they are off unless asked for:
/// `--time-log`, `--debug`, or `FLEXALIGN_TIMING=1`. Off, `time()` reports zero.
pub static STAGE_TIMING: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub fn set_stage_timing(on: bool) {
    STAGE_TIMING.store(on, std::sync::atomic::Ordering::Relaxed);
}

/// Always measures: for the few one-off timings (index and reference load).
#[inline]
pub fn time<F, T>(f: F) -> (Duration, T)
    where F: FnOnce() -> T {
    let start: Instant = Instant::now();
    let result = f();
    (start.elapsed(), result)
}

/// Measures only when the stage timers are on: for the per-read pipeline stages.
#[inline]
pub fn stage_time<F, T>(f: F) -> (Duration, T)
    where F: FnOnce() -> T {
    if !STAGE_TIMING.load(std::sync::atomic::Ordering::Relaxed) {
        return (Duration::ZERO, f());
    }
    let start: Instant = Instant::now();
    let result = f();
    (start.elapsed(), result)
}

// pub fn run_worker(args: Args) {
    
// }

/// The (reads_1, reads_2, output) triples a sharded run should process.
///
/// Output path follows the same rules as the unsharded path: `--output` names a file for a single
/// input pair and a directory for several, and the format is taken from the extension (`.sam` by
/// default, `.paf` with --paf) so a caller does not have to restate it.
fn sharded_pairs(options: &Options) -> Vec<(String, String, String)> {
    let ext = if options.args.emit_sam() { "sam" } else { "paf" };
    let fwd: Vec<String> = options.args.fwd.iter().filter(|f| !f.is_empty()).cloned().collect();
    let rev: Vec<String> = options.args.rev.iter().filter(|f| !f.is_empty()).cloned().collect();
    if fwd.is_empty() || fwd.len() != rev.len() {
        return Vec::new();
    }
    let multi = fwd.len() > 1;
    fwd.iter()
        .zip(rev.iter())
        .map(|(r1, r2)| {
            let out = match (&options.args.output, multi) {
                // One pair, explicit output: use it verbatim.
                (Some(o), false) => o.clone(),
                // Several pairs: --output is a directory, prefix inferred per input.
                (Some(dir), true) => {
                    let _ = std::fs::create_dir_all(dir);
                    let stem = std::path::Path::new(r1)
                        .file_name()
                        .map(|s| s.to_string_lossy().split('.').next().unwrap_or("out").to_string())
                        .unwrap_or_else(|| "out".to_string());
                    format!("{}/{}.{}", dir.trim_end_matches('/'), stem, ext)
                }
                // No --output: next to the reads.
                (None, _) => format!("{}.{}", r1.trim_end_matches(".gz"), ext),
            };
            (r1.clone(), r2.clone(), out)
        })
        .collect()
}

pub fn run(args: Args) {
    let options = Options::from_args(args);
    crate::align::common::MAPQ_PER_MISMATCH.store(options.args.mapq_per_mismatch, std::sync::atomic::Ordering::Relaxed);
    crate::align::common::MAPQ_UNIQUE_SETTING.store(options.args.mapq_unique, std::sync::atomic::Ordering::Relaxed);
    crate::align::modular_workflow::DISCORDANT_MAPQ_CAP.store(options.args.discordant_mapq_cap, std::sync::atomic::Ordering::Relaxed);
    crate::align::modular_workflow::DISCORDANT_ONLY_MAPPED.store(options.args.discordant_only_mapped, std::sync::atomic::Ordering::Relaxed);
    crate::align::common::MAPQ_UNIQUE_SCALE.store((options.args.mapq_unique_scale as f32).to_bits(), std::sync::atomic::Ordering::Relaxed);
    if let Some(path) = options.args.dump_candidates.as_ref() {
        use std::io::Write as _;
        let file = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(path).expect("--dump-candidates file"));
        // A .gz path is written compressed (fast level: the dump is large).
        let mut w: Box<dyn std::io::Write + Send> = if path.ends_with(".gz") {
            Box::new(flate2::write::GzEncoder::new(file, flate2::Compression::fast()))
        } else {
            Box::new(file)
        };
        w.write_all(crate::align::modular_workflow::DUMP_HEADER.as_bytes()).unwrap();
        let _ = crate::align::modular_workflow::DUMP.set((std::sync::Mutex::new(Some(w)), options.args.dump_every.max(1)));
    }
    let _dump_finish = crate::align::modular_workflow::DumpFinish;
    crate::align::common::MAPQ_FIT_PENALTY.store((options.args.mapq_fit_penalty as f32).to_bits(), std::sync::atomic::Ordering::Relaxed);
    crate::align::common::MAPQ_MM_SCALE.store((options.args.mapq_mm_scale as f32).to_bits(), std::sync::atomic::Ordering::Relaxed);
    crate::align::common::UNPAIRED_PENALTY.store(options.args.unpaired_penalty.max(0) as u32, std::sync::atomic::Ordering::Relaxed);

    // Set before any database load -- the mapping happens below the option plumbing.
    {
        let spec = options.args.ref_budget.clone()
            .or_else(|| options.args.lazy_ref.then(|| "auto".to_string()));
        let kb = crate::database::common::resolve_ref_budget(spec.as_deref());
        crate::database::common::REF_BUDGET_KB.store(kb, std::sync::atomic::Ordering::Relaxed);
        if kb > 0 {
            eprintln!("--ref-budget: capping resident reference at {:.1} GB", kb as f64 / 1048576.0);
        }
    }

    // The whole database is read up front unless --lazy-ref says otherwise. Set before any load,
    // for the same reason as the budget above: the mapping happens below the option plumbing.
    crate::database::common::EAGER_REF
        .store(!options.args.lazy_ref, std::sync::atomic::Ordering::Relaxed);

    // The packed index with the line-hashed key table, for the whole index and for shards alike.
    crate::database::common::set_index_format(IndexFormat::Hashed);

    if !options.reference.exists() {
        eprintln!("Reference does not exist {:?}", options.reference);
        exit(9);
    }

    let db_paths = DBPaths::new(&options.reference);

    const K: usize = 31;
    const C: usize = 15;
    const F: usize = 16;
    const S: usize = 4; // 7 0.34 //8 0.37 //6 0.31 //5  0.29 //4 0.289 //3 0.413  //2  0.413
    const L: usize = C - S + 1; //1
    const CELLS_PER_BODY: u64 = 16;
    const HEADER_THRESHOLD: usize = 2;

    // Shards are windows of a packed index's bundle (`shard::slice`): their bounds are
    // discoverable from the reference path alone, so the sharded pipeline is selected here rather
    // than by reaching for a separate executable. --shards N picks a shard count (its bounds are
    // computed on first use); --no-shards ignores them. An index not built yet is built first, below.
    let format = IndexFormat::Hashed;
    let bundle_ready = format != IndexFormat::Blob
        && std::path::Path::new(&db_paths.bundle_file(format)).exists()
        && DB::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>::index_compatible(&db_paths)
        && !options.args.force_build;
    if bundle_ready {
        if let Some(n_shards) = options.args.shard_slice {
            shard_bounds::<C, F, HEADER_THRESHOLD>(&db_paths, format, n_shards, options.args.threads as usize);
            return;
        }
        if !options.args.no_shards {
            let available = crate::shard::slice::ShardBounds::available(&db_paths, format);
            let chosen = options.args.shards.or_else(|| match available.len() {
                0 => None,
                1 => Some(available[0]),
                _ => {
                    // Several shard counts have bounds and none was named. Picking one silently
                    // would make the run's memory profile depend on what happens to be on disk, so
                    // say what is there and let the caller choose.
                    eprintln!(
                        "Shard bounds exist for several counts ({}). Pick one with --shards N, or --no-shards.",
                        available.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(", ")
                    );
                    exit(3);
                }
            });
            // One shard is the whole index: the unsharded run does that in one pass over the reads,
            // the sharded pipeline in two plus an evidence file.
            if let Some(n_shards) = chosen.filter(|&n| n > 1) {
                run_sharded::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(&db_paths, &options, n_shards);
                return;
            }
        }
    }

    // Build if the index is missing, if forced, OR if the on-disk blob is incompatible with this
    // build (wrong blob version / const params). The last check is what stops a stale index -- e.g.
    // one written before a flexmap key-format change -- from being memory-mapped and silently
    // misread (which collapses sensitivity to ~zero), forcing a clean rebuild instead.
    let index_compatible = DB::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>::index_compatible(&db_paths);
    if db_paths.names_valid() && !index_compatible && !options.args.force_build {
        log::warn!("On-disk index is incompatible with this build (format/params changed); rebuilding it.");
    }
    let build = options.args.force_build || !index_compatible;

    // Build or load database
    let db: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = match build {
        true => {
            
            let (_duration, result) = 
                time(|| DB::build(&options));
            let _ = result.save(&db_paths, GLOBAL_VERSION);
            result

        },
        false => {
            eprintln!("Load index.");
            let (duration, result) =
                time(|| {
                    DB::load(&db_paths, GLOBAL_VERSION)
                });
            eprintln!("Loading index took: {:?}", duration);
            result
        },
    };

    // Sharding asked for on an index that had to be built first: now that its bundle is written,
    // compute the bounds (--shard-slice) or run on the shards (--shards).
    let shards = options.args.shards.filter(|&n| !options.args.no_shards && n > 1);
    if options.args.shard_slice.is_some() || shards.is_some() {
        drop(db);
        match options.args.shard_slice {
            Some(n_shards) => shard_bounds::<C, F, HEADER_THRESHOLD>(&db_paths, format, n_shards, options.args.threads as usize),
            None => run_sharded::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(&db_paths, &options, shards.unwrap()),
        }
        return;
    }

    // Check if all files exist
    for file in &options.fwd {
        if !file.exists() {
            panic!("File passed with --rev does not exist: \n{}", file.to_str().unwrap());
        }
    }

    for file_option in &options.rev {
        match file_option {
            Some(file) => if !file.exists() {
                panic!("File passed with --rev does not exist: \n{}", file.to_str().unwrap());
            },
            None => {},
        }
    }

    // Prepare the --time-log TSV (truncate + header) before the per-input loop appends sections.
    if let Some(path) = &options.args.time_log {
        if let Err(e) = crate::align::stats::Stats::init_time_log(path) {
            eprintln!("warning: could not create --time-log {}: {}", path, e);
        }
    }

    // Timed directly, not through `time()`: that one is opt-in now, and this single measurement
    // is what the benchmarks read.
    let started = Instant::now();
    process_fastq_wrapper_modular::<K, C, F, S, L, HEADER_THRESHOLD,DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>>(&options, &db);
    log::info!("Modular: Process reads: {:?}", started.elapsed());
}

/// --shard-slice N: computes (or confirms) the bounds of N shards of the bundle of `format`.
fn shard_bounds<const C: usize, const F: usize, const T: usize>(db_paths: &DBPaths, format: IndexFormat, n_shards: usize, threads: usize) {
    let (duration, (bounds, computed)) = time(|| {
        crate::shard::slice::ShardBounds::load_or_compute::<C, F, T>(db_paths, format, n_shards, threads.max(1))
            .unwrap_or_else(|e| panic!("shard bounds of {}: {e}", db_paths.bundle_file(format)))
    });
    eprintln!(
        "{} bounds of {} shard(s) in {:?}: {}",
        if computed { "Computed" } else { "Reusing" },
        bounds.ranges().len(),
        duration,
        crate::shard::slice::ShardBounds::path_for(db_paths, format, n_shards)
    );
}

/// --shards N (or the one shard count with bounds): aligns paired-end reads on N windows of the
/// bundle, one after another.
fn run_sharded<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(db_paths: &DBPaths, options: &Options, n_shards: usize) {
    let pairs = sharded_pairs(options);
    if pairs.is_empty() {
        eprintln!("Sharded alignment is paired-end only: give -1 and -2.");
        exit(4);
    }
    eprintln!("Using {n_shards} shards of {:?}", options.reference);
    crate::shard::run::run_sharded::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(db_paths, options, n_shards, pairs);
}
