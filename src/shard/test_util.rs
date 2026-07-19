//! Shared test-only helpers for the shard modules (mock index + FASTQ generators).
//!
//! The mock reports every core-mer as present with a single, key-dependent headerless cell. That
//! is enough to exercise kmer extraction -> range lookup -> emit -> group -> replay without building
//! a real 2^(2C) flexmap, and the key-dependent cell means different core-mers yield different
//! seeds, so a seed-level parity check is meaningful rather than vacuous.

use std::io::Write;

use flexmap::values::{VCell, VRange};

use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::options::Options;

// Small scheme so the key space (2^(2C)) is tiny and splits many ways: K = C + 2F.
pub const K: usize = 8;
pub const C: usize = 4;
pub const F: usize = 2;
pub const S: usize = 2;
pub const L: usize = C - S + 1; // 3
pub const KEY_SPACE: u64 = 1 << (2 * C); // 256

/// A stand-in index: every core-mer is present, with one distinct headerless cell per key. Header
/// lookups (`get_rname` etc.) are never touched by the pipeline front-half, so they are stubbed.
#[derive(Clone)]
pub struct MockDb {
    cells: Vec<VCell>,
}

impl MockDb {
    pub fn new() -> Self {
        // One distinct, non-zero cell per core-mer.
        let cells = (0..KEY_SPACE).map(|k| VCell(((k + 1) << 8) | 0x5)).collect();
        Self { cells }
    }
}

impl FlexalignDatabase for MockDb {
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange<'_>> {
        let i = canonical_kmer as usize;
        // A shard restricts keys before delegating, and core-mers are < KEY_SPACE, so this is in
        // bounds; guard anyway rather than risk a panic corrupting a test.
        self.cells.get(i..i + 1).map(|slice| VRange::new(None, slice))
    }
    fn get_rid(&self, _reference: &str) -> Option<&usize> {
        None
    }
    fn get_rname(&self, _id: usize) -> Option<&str> {
        None
    }
    fn get_reference(&self, _id: usize) -> Option<&[u8]> {
        None
    }
    fn build(_options: &Options) -> Self {
        unimplemented!("mock is not built from a reference")
    }
    fn save(&self, _paths: &DBPaths, _version: u32) -> Result<(), std::io::Error> {
        unimplemented!()
    }
    fn load(_paths: &DBPaths, _version: u32) -> Self {
        unimplemented!()
    }
}

fn push_record(buf: &mut Vec<u8>, i: usize, mate: u8) {
    let len = 40 + (i * 3) % 20;
    let seq: Vec<u8> = (0..len).map(|j| b"ACGT"[(i * 5 + j * 3) % 4]).collect();
    let qual = vec![b'I'; len];
    write!(buf, "@read{}/{}\n", i, mate).unwrap();
    buf.extend_from_slice(&seq);
    buf.push(b'\n');
    buf.extend_from_slice(b"+\n");
    buf.extend_from_slice(&qual);
    buf.push(b'\n');
}

/// A paired FASTQ (R1, R2) with `n` records each.
pub fn make_pair(n: usize) -> (Vec<u8>, Vec<u8>) {
    let (mut d1, mut d2) = (Vec::new(), Vec::new());
    for i in 0..n {
        push_record(&mut d1, i, 1);
        push_record(&mut d2, i, 2);
    }
    (d1, d2)
}
