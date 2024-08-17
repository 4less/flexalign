use std::{collections::HashMap, fs::File, io::BufWriter};

use flexmap::flexmap::{Flexmap, FlexmapHash};
use savefile::save;
use bincode::{self, config};

use crate::{database::common::DBPaths, flexalign::time, options::{self, Options}};

pub fn default<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize
>(options: &Options) -> Result<(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {
    
    let db_paths = DBPaths::new(&options.reference);

    let (duration, result) = time(|| {
        flexmap::build::default_build::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(
            &options.reference, options.args.max_range_size
        )
    });
    println!("Building database: {:?}", duration);

    let (flexmap, reference2id, id2reference) = result.expect("Building database works");

    const GLOBAL_VERSION: u32 = 1;
    let mut file = match File::create(&db_paths.index_path) {
        Err(why) => panic!("couldn't open {}: {}", db_paths.index_path.display(), why),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &flexmap);

    let mut file = match File::create(&db_paths.id2reference_path) {
        Err(why) => panic!(
            "couldn't open {}: {}",
            db_paths.id2reference_path.display(),
            why
        ),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &id2reference);


    let mut file = match File::create(&db_paths.reference2id_path) {
        Err(why) => panic!(
            "couldn't open {}: {}",
            db_paths.id2reference_path.display(),
            why
        ),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &reference2id);


    return Ok((flexmap, reference2id, id2reference));
}


pub fn hash<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize
>(options: &Options) -> Result<(FlexmapHash<C, F, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {
    
    let db_paths = DBPaths::new(&options.reference);

    let (duration, result) = time(|| {
        flexmap::build::hash_build::<K, C, F, S, L, HEADER_THRESHOLD>(
            &options.reference, options.args.max_range_size
        )
    });
    eprintln!("Building database: {:?}", duration);

    let (flexmap, reference2id, id2reference) = result.expect("Building database works");

    const GLOBAL_VERSION: u32 = 1;

    let mut file: File = match File::create(&db_paths.index_path) {
        Err(why) => panic!("couldn't open {}: {}", db_paths.index_path.display(), why),
        Ok(file) => file,
    };

    let mut writer = BufWriter::new(file);
    // let config = bincode::config::standard();
    // let _ = bincode::encode_into_std_write(&flexmap, &mut writer, config);

    let _ = save(&mut writer, GLOBAL_VERSION, &flexmap);


    eprintln!("Done saving..");

    let mut file = match File::create(&db_paths.id2reference_path) {
        Err(why) => panic!(
            "couldn't open {}: {}",
            db_paths.id2reference_path.display(),
            why
        ),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &id2reference);


    let mut file = match File::create(&db_paths.reference2id_path) {
        Err(why) => panic!(
            "couldn't open {}: {}",
            db_paths.id2reference_path.display(),
            why
        ),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &reference2id);


    return Ok((flexmap, reference2id, id2reference));
}
