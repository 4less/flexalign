#![feature(test)]
#![feature(slice_internals)]
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![feature(map_try_insert)]
#![feature(fn_traits, unboxed_closures)]

extern crate test;

#[macro_use]
extern crate savefile_derive;

pub mod database;
pub mod align;
pub mod flexalign;
pub mod options;
pub mod io;
// `misc` was a KHashEntry (de)serialization micro-benchmark; flexmap dropped the hash key type.
// Dead and uncompiled — retained on disk. // pub mod misc;
pub mod utils;
pub mod distance;
pub mod shard;

const GLOBAL_VERSION: u32 = 1;

const GOLDSTD_EVAL_ENV_VAL: Option<&str> = option_env!("FLEXALIGN_GOLDSTD_EVAL");
pub const GOLDSTD_EVAL: bool = match GOLDSTD_EVAL_ENV_VAL {
    Some(val) => val.as_bytes().len() > 0 && val.as_bytes()[0] == b'1',
    None => false,
};

const STAGE_TRACE_ENV_VAL: Option<&str> = option_env!("FLEXALIGN_STAGE_TRACE");
/// Per-stage timing *inside* the alignment, plus which stage dropped an anchor.
///
/// Compile-time, like [`GOLDSTD_EVAL`], and for the same reason: the alignment stages are short
/// enough that reading the clock around each one is not free. Five `Instant::now()` pairs, over
/// `--align-top-y` anchors, over two mates, is ~40 clock reads per read — several seconds across a
/// 10 M-pair sample, charged to the very thing being measured. Off, every branch below folds away
/// and `AlignTrace` costs nothing.
///
/// Implied by `GOLDSTD_EVAL`, because the truth-survival funnel needs the drop stage to attribute a
/// loss to the filter that caused it — a gold-standard build gets the whole funnel with no second
/// flag to remember.
///
/// ```bash
/// FLEXALIGN_STAGE_TRACE=1 cargo +nightly build --release   # timings only
/// FLEXALIGN_GOLDSTD_EVAL=1 cargo +nightly build --release  # timings + survival funnel
/// ```
pub const STAGE_TRACE: bool = GOLDSTD_EVAL
    || match STAGE_TRACE_ENV_VAL {
        Some(val) => val.as_bytes().len() > 0 && val.as_bytes()[0] == b'1',
        None => false,
    };


pub fn add(left: usize, right: usize) -> usize {
    left + right
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
