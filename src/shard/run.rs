//! The sharded alignment run, as a library entry point.
//!
//! It lives here, in the library, so the MAIN binary can select it automatically: a sliced index is
//! discoverable from the reference path (`<index>.s<n>.shards.json`, see
//! [`SliceManifest::available`]), so a caller should not have to know which executable to reach for.
//! There is deliberately no separate sharded binary -- one entry point, one set of options.
//!
//! Shape of the run, and why: the SHARD loop is the outer one and samples are inner, so each shard
//! is memory-mapped once for the whole batch instead of once per sample, and is dropped before the
//! next loads. Evidence is streamed to per-(shard, sample) files rather than accumulated.
//!
//! The memory contract, which the whole design exists to provide:
//! **at any instant either ONE shard is mapped, or the reference is -- never two shards, and never a
//! shard together with the reference.** The shard loop drops each shard before the next loads, and the
//! reference is opened only afterwards, for the join. Reads are STREAMED from their files on every
//! pass and never held decompressed, so peak RSS is `max(one shard, reference)` and does not grow
//! with the number of samples or the depth of the reads.

use std::fs::File;
use std::io::{BufReader, Read};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use bioreader::utils::is_gzip;
use flate2::read::GzDecoder;

use crate::database::common::DBPaths;
use crate::database::flexmap::DB;
use crate::io::output_buffer::{OutputBuffer, OutputTarget};
use crate::options::Options;
use crate::shard::disk::{count_batches, run_disk_rejoin_align, run_shard_pass_to_disk};
use crate::shard::slice::{slice_index, PhysicalShardDB, SliceManifest};

const BUFFER_SIZE: usize = 1 << 24;
const MAX_HEADERED: usize = 15;

/// One sample of a batch: its read FILES and its output path. The reads are never held in memory.
pub struct ShardSample {
    pub r1: String,
    pub r2: String,
    pub out: String,
    n_batches: u64,
}

/// A fresh streaming reader over a read file, transparently gunzipping.
///
/// Every pass over a sample re-opens its files instead of reusing a decompressed copy. Holding the
/// reads was the single largest term in a sharded batch's memory: a 10 M-pair sample is ~6.5 GB
/// DECOMPRESSED (0.5 GB gzipped), so five samples held at once were 32.6 GB -- more than the 15.4 GB
/// reference and four times one 7.5 GB shard. That inverted the whole point of sharding, making the
/// sharded batch the largest-RSS contender of all.
///
/// The cost of streaming is re-inflating the gzip once per shard pass (CPU), which is the correct
/// trade here: the design exists to bound memory, and peak RSS must be
/// `max(one shard, reference)` -- never a shard plus a reference plus every sample's reads.
fn open_reads(path: &str) -> Box<dyn Read + Send> {
    let file = File::open(path).unwrap_or_else(|e| panic!("open {}: {}", path, e));
    if is_gzip(path).unwrap_or(false) {
        Box::new(GzDecoder::new(BufReader::with_capacity(1 << 22, file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 22, file))
    }
}

/// Align every `(reads_1, reads_2, output)` triple against the `n_shards` slicing of `paths`.
///
/// Slices the index first if that slicing does not exist yet -- the only step that loads the full
/// index, and a one-time cost per shard count.
pub fn run_sharded<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(
    paths: &DBPaths,
    options: &Options,
    n_shards: usize,
    pairs: Vec<(String, String, String)>,
) {
    type ShardOf<const C: usize, const F: usize, const CB: u64, const HT: usize> =
        PhysicalShardDB<C, F, CB, HT>;

    // -t, not a constant. This was hardcoded to 4, so sharded runs ignored --threads entirely and
    // were benchmarked at 4 workers against everything else's 16 -- most of an apparent slowdown.
    let threads: u32 = options.args.threads.max(1);

    let manifest = match SliceManifest::load(&SliceManifest::path_for(paths, n_shards)) {
        Ok(m) if m.shards.len() == n_shards => {
            eprintln!("Reusing existing slice ({n_shards} shards)");
            m
        }
        _ => {
            eprintln!("Slicing index into {n_shards} physical shards ...");
            slice_index::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(paths, n_shards).expect("slice")
        }
    };

    // Only the FILE NAMES are kept; each pass re-opens and streams them (see `open_reads`). The
    // batch amortises the shard loads, not the reads -- holding the reads is what blew peak RSS past
    // the unsharded contender it is supposed to undercut.
    let samples: Vec<ShardSample> = pairs
        .into_iter()
        .map(|(r1, r2, out)| {
            let n_batches = count_batches(open_reads(&r1), open_reads(&r2), BUFFER_SIZE, threads);
            ShardSample { r1, r2, out, n_batches }
        })
        .collect();
    eprintln!(
        "{} sample(s), batches: {}",
        samples.len(),
        samples.iter().map(|s| s.n_batches.to_string()).collect::<Vec<_>>().join(", ")
    );

    // Evidence is per (shard, sample), keyed on the sample's own output path so concurrent runs and
    // re-runs cannot collide.
    let evidence: Vec<Vec<String>> = (0..manifest.shards.len())
        .map(|i| samples.iter().map(|s| format!("{}.evidence.s{i}.bin", s.out)).collect())
        .collect();

    let t_pass = Instant::now();
    for i in 0..manifest.shards.len() {
        let shard: ShardOf<C, F, CELLS_PER_BODY, HEADER_THRESHOLD> =
            ShardOf::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::load(&manifest, i);
        for (j, s) in samples.iter().enumerate() {
            let emit = run_shard_pass_to_disk::<
                K, C, F, S, L,
                ShardOf<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
                Box<dyn Read + Send>,
            >(
                open_reads(&s.r1), open_reads(&s.r2), &shard, s.n_batches, &evidence[i][j],
                BUFFER_SIZE, threads, options.args.max_best_flex, MAX_HEADERED,
            )
            .expect("shard pass");
            eprintln!(
                "  shard {i} / sample {j} ({}): emit headered={} headerless={}",
                s.r1, emit.headered, emit.headerless
            );
        }
        // `shard` drops here (mmap unmaps), freeing it before the next shard loads.
    }
    let pass_secs = t_pass.elapsed().as_secs_f64();

    // Rejoin + align: streams evidence back per batch and REPLAYS the anchors the shard passes
    // already found (run_disk_rejoin_align never calls a range extractor), so the k-mer index is not
    // needed here at all -- only the reference bases to align against. Loading the full 30 GB index,
    // as this used to, read a file the join never queries. The join phase touches NO index: it loads
    // the reference only, mapping the index without reading it.
    eprintln!("Loading reference for the align half (no index) ...");
    let full: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = DB::load_reference_only(paths);

    let t_align = Instant::now();
    let mut total_alignments = 0usize;
    for (j, s) in samples.iter().enumerate() {
        let sink = Arc::new(Mutex::new(OutputTarget::File(
            File::create(&s.out).unwrap_or_else(|e| panic!("create {}: {}", s.out, e)),
        )));
        let out_buffer = OutputBuffer::new(Arc::clone(&sink), BUFFER_SIZE);
        let ev: Vec<String> = (0..manifest.shards.len()).map(|i| evidence[i][j].clone()).collect();
        let stats = run_disk_rejoin_align::<
            K, C, F, S, L,
            DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>,
            Box<dyn Read + Send>,
        >(
            open_reads(&s.r1), open_reads(&s.r2), &full, options, out_buffer, &ev,
            options.args.ranges, BUFFER_SIZE, threads,
        )
        .expect("rejoin");
        total_alignments += stats.alignments_successful as usize;
        eprintln!("  sample {j}: alignments={} -> {}", stats.alignments_successful, s.out);
        for p in &ev {
            let _ = std::fs::remove_file(p);
        }
    }
    let align_secs = t_align.elapsed().as_secs_f64();

    eprintln!(
        "\nsharded done: {} sample(s), {} shard(s); passes {:.2}s + rejoin/align {:.2}s ; alignments={}",
        samples.len(),
        manifest.shards.len(),
        pass_secs,
        align_secs,
        total_alignments
    );
}
