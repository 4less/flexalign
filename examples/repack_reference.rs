//! Writes a bundle again with the current format: a reference of short sequences (a marker
//! database) at 2 bits a base (`refseq2`).
//!
//! Run:
//!   cargo +nightly run --release --example repack_reference -- \
//!       <reference.fna of the bundle> <reference.fna path to write next to> [hashed]
//!
//! The bundle of the first is read (as an alignment run reads it) and written as the bundle of the
//! second, so a built index can be converted without building it again.

use std::path::PathBuf;

use flexalign::database::common::{set_index_format, DBPaths, FlexalignDatabase, IndexFormat};
use flexalign::database::flexmap::DB;

type Db = DB<31, 15, 16, 4, 12, 16, 2>;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: {} <reference.fna> <out reference.fna> [hashed]", args[0]);
        std::process::exit(2);
    }
    let format = match args.get(3).map(String::as_str) {
        None | Some("hashed") => IndexFormat::Hashed,
        Some(other) => panic!("unknown index format {other}"),
    };
    set_index_format(format);
    let db = Db::load(&DBPaths::new(PathBuf::from(&args[1])), 1);
    let out = DBPaths::new(PathBuf::from(&args[2]));
    db.save(&out, 1).expect("writing the bundle");
    eprintln!("wrote {}", out.bundle_file(format));
}
