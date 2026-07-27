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

pub fn run(args: Args) {
    let options = Options::from_args(args);


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
            crate::shard::slice::SliceManifest::path_for(&db_paths)
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
