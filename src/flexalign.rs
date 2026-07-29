use std::time::{Duration, Instant};
use std::process::exit;

use crate::align::process_fastq::process_fastq_wrapper_modular;
use crate::database::flexmap::DB;
use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::options::{Args, Options};
use crate::GLOBAL_VERSION;


pub fn time<F, T>(f: F) -> (Duration, T)
    where F: FnOnce() -> T {

    let start: Instant = Instant::now();
    let result = f();
    (start.elapsed(), result)
}

// pub fn run_worker(args: Args) {
    
// }

/// The (reads_1, reads_2, output) triples a sharded run should process.
///
/// Output path follows the same rules as the unsharded path: `--output` names a file for a single
/// input pair and a directory for several, and the format is taken from the extension (`.sam` with
/// --sam, `.paf` otherwise) so a caller does not have to restate it.
fn sharded_pairs(options: &Options) -> Vec<(String, String, String)> {
    let ext = if options.args.sam { "sam" } else { "paf" };
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

    // A sliced index is discoverable from the reference path alone, so the sharded pipeline is
    // selected here rather than by reaching for a separate executable. --shards N picks a specific
    // slicing (or slices one); --no-shards ignores them.
    if options.args.shard_slice.is_none() && !options.args.no_shards {
        let available = crate::shard::slice::SliceManifest::available(&db_paths);
        let chosen = options.args.shards.or_else(|| {
            match available.len() {
                0 => None,
                1 => Some(available[0]),
                _ => {
                    // Several slicings exist and none was named. Picking one silently would make the
                    // run's memory profile depend on what happens to be on disk, so say what is
                    // there and let the caller choose.
                    eprintln!(
                        "Several sliced indexes exist ({}). Pick one with --shards N, or --no-shards.",
                        available.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(", ")
                    );
                    exit(3);
                }
            }
        });
        if let Some(n_shards) = chosen {
            let pairs = sharded_pairs(&options);
            if pairs.is_empty() {
                eprintln!("Sharded alignment is paired-end only: give -1 and -2.");
                exit(4);
            }
            eprintln!("Using the {n_shards}-shard index for {:?}", options.reference);
            crate::shard::run::run_sharded::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(
                &db_paths, &options, n_shards, pairs,
            );
            return;
        }
    }

    // Build if the index is missing, if forced, OR if the on-disk blob is incompatible with this
    // build (wrong blob version / const params). The last check is what stops a stale index -- e.g.
    // one written before a flexmap key-format change -- from being memory-mapped and silently
    // misread (which collapses sensitivity to ~zero), forcing a clean rebuild instead.
    let index_compatible =
        db_paths.valid_paths() && DB::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>::index_compatible(&db_paths);
    if db_paths.valid_paths() && !index_compatible && !options.args.force_build {
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
                time(|| DB::load(&db_paths, GLOBAL_VERSION));
            eprintln!("Loading index took: {:?}", duration);
            result
        },
    };

    // Shard-slicing mode: cut the built index into N physical shard files + a manifest, then exit.
    // The index is already built/loaded above, so `.flex.index.keys` exists for the slicer to read.
    if let Some(n_shards) = options.args.shard_slice {
        drop(db); // the slicer re-reads the keys file; free the in-RAM copy first
        let (duration, manifest) = time(|| {
            crate::shard::slice::slice_index::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(&db_paths, n_shards)
                .expect("slicing the index")
        });
        eprintln!(
            "Sliced index into {} shard(s) in {:?}: {}",
            manifest.shards.len(),
            duration,
            crate::shard::slice::SliceManifest::path_for(&db_paths, n_shards)
        );
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

    let (duration, _result) = time(|| process_fastq_wrapper_modular::<K, C, F, S, L, HEADER_THRESHOLD,DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>>(&options, &db));
    log::info!("Modular: Process reads: {:?}", duration);
}
