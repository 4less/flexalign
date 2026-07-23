//! Key-sharded prefilter with deferred rejoin. See ProjectShard.md.
//!
//! The flexmap is sliced by core-mer into N shards, each small enough to be resident. A shard
//! pass runs the pipeline down to the range lookup and writes what it found; the rejoin reads
//! every shard's evidence back, rebuilds the seed list per read, and runs the reference-based
//! stages once.
//!
//! Status: `record` and `chunkfile` are the intermediate format and its container; `db` is the
//! key-range view a shard pass aligns against; `manifest` is what every pass must agree on;
//! `emit` and `replay` are the two sides of the cut. The pass drivers -- the loop that runs a
//! shard over a FASTQ, and the rejoin that streams the FASTQ back alongside the evidence -- are
//! not implemented yet, and neither is the read-id plumbing they need from bioreader (§6).

pub mod acceptance;
pub mod chunkfile;
pub mod db;
pub mod disk;
pub mod emit;
pub mod manifest;
pub mod pass;
pub mod record;
pub mod rejoin;
pub mod rejoin_align;
pub mod replay;
pub mod slice;

#[cfg(test)]
pub(crate) mod test_util;
