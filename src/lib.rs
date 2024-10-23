#![feature(slice_internals)]
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![feature(map_try_insert)]
#![feature(fn_traits, unboxed_closures)]

#[macro_use]
extern crate savefile_derive;

pub mod database;
pub mod align;
pub mod flexalign;
pub mod options;
pub mod io;
pub mod misc;
pub mod utils;

const GLOBAL_VERSION: u32 = 1;

const GOLDSTD_EVAL_ENV_VAL: Option<&str> = option_env!("FLEXALIGN_GOLDSTD_EVAL");
pub const GOLDSTD_EVAL: bool = match GOLDSTD_EVAL_ENV_VAL {
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
