use std::time::{Duration, Instant};
use std::process::exit;

use crate::align::process_fastq::process_fastq_wrapper_modular;
use crate::database::flexmap::DB;
use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::options::{Args, Options, Profile};
use crate::GLOBAL_VERSION;

/// Expands a runtime [`Profile`] into a call to a generic function with the
/// matching compile-time const-generic parameters. Each arm is a separate
/// monomorphization; this `match` is the single runtime dispatch point and runs
/// once per invocation, before the parallel read loop. Adding a variant = one
/// arm here plus one [`Profile`] variant. The `L = C - S + 1` invariant lives
/// only in these tuples.
macro_rules! with_profile {
    ($profile:expr, $f:ident, $($arg:expr),* $(,)?) => {
        match $profile {
            //                   K   C   F  S   L   CELLS HDR
            Profile::Default => $f::<31, 15, 16, 4, 12,  16,  2>($($arg),*),
            // Profile::Mid  => $f::<25, 13, 14, 3, 11,  16,  2>($($arg),*),
            // Profile::Long => $f::<19, 11, 12, 2, 10,  16,  2>($($arg),*),
        }
    };
}


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

    // Single runtime dispatch point: pick the const-generic pipeline for the
    // selected profile, then run the fully-monomorphized core.
    with_profile!(options.args.profile, run_with_params, &options);
}

/// The monomorphized core of `run`: build/load the index and process the reads
/// for one fixed set of const-generic parameters. `run` selects the parameter
/// set at runtime via [`with_profile!`].
fn run_with_params<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(options: &Options) {
    let db_paths = DBPaths::new(&options.reference);

    let build = !db_paths.valid_paths() || options.args.force_build;

    // Build or load database
    let db: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = match build {
        true => {
            let (_duration, result) =
                time(|| DB::build(options));
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

    let (duration, _result) = time(|| process_fastq_wrapper_modular::<K, C, F, S, L, HEADER_THRESHOLD, DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>>(options, &db));
    log::info!("Modular: Process reads: {:?}", duration);
}
