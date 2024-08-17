use std::time::{Duration, Instant};
use std::process::exit;
use log::info;

use crate::align::process_fastq::process_fastq_wrapper;
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

pub fn run_worker(args: Args) {
    
}

pub fn run(args: Args) {
    let options = Options::from_args(args);

    if !options.reference.exists() {
        eprintln!("Invalid reference");
        exit(9);
    }

    let db_paths = DBPaths::new(&options.reference);
    
    let build = !db_paths.valid_paths() || options.args.force_build;
    
    const K: usize = 31;
    const C: usize = 15;
    const F: usize = 16; 
    const S: usize = 7; // 7 0.34 //8 0.37 //6 0.31 //5  0.29 //4 0.289 //3 0.413  //2  0.413
    const L: usize = C - S + 1; //1
    const CELLS_PER_BODY: u64 = 16;
    const HEADER_THRESHOLD: usize = 2;
    
    let db: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = match build {
        true => {
            
            let (duration, result) = 
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

    let (duration, _result) = time(|| process_fastq_wrapper::<K, C, F, S, L, HEADER_THRESHOLD,DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>>(&options, &db));
    eprintln!("Process reads: {:?}", duration);
}

