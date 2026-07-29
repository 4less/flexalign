//! The sharded alignment run, as a library entry point.
//!
//! This is the pipeline `examples/shard_align` used to own. It lives here so the main binary can
//! select it automatically: a sliced index is discoverable from the reference path
//! (`<index>.s<n>.shards.json`, see [`SliceManifest::available`]), so callers should not have to
//! know which executable to reach for. The example is now a thin wrapper over this.
//!
//! Shape of the run, and why: the SHARD loop is the outer one and samples are inner, so each shard
//! is memory-mapped once for the whole batch instead of once per sample, and is dropped before the
//! next loads. Evidence is streamed to per-(shard, sample) files rather than accumulated, so peak
//! memory is one shard plus one batch of evidence -- independent of how many reads there are.

use std::fs::File;
use std::io::{BufReader, Read};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use bioreader::utils::is_gzip;
use flate2::read::GzDecoder;

use crate::database::common::{DBPaths, FlexalignDatabase};
use crate::database::flexmap::DB;
use crate::io::output_buffer::{OutputBuffer, OutputTarget};
use crate::options::Options;
use crate::shard::disk::{count_batches, run_disk_rejoin_align, run_shard_pass_to_disk};
use crate::shard::slice::{slice_index, PhysicalShardDB, SliceManifest};

const BUFFER_SIZE: usize = 1 << 24;
const THREADS: u32 = 4;
const MAX_HEADERED: usize = 15;

/// One sample of a batch: its read files, its output path, and its reads once loaded.
pub struct ShardSample {
    pub r1: String,
    pub r2: String,
    pub out: String,
    reads1: Vec<u8>,
    reads2: Vec<u8>,
    n_batches: u64,
}

fn read_all(path: &str) -> Vec<u8> {
    let gz = is_gzip(path).unwrap_or(false);
    let file = File::open(path).unwrap_or_else(|e| panic!("open {}: {}", path, e));
    let mut buf = Vec::new();
    if gz {
        GzDecoder::new(file).read_to_end(&mut buf).expect("gunzip");
    } else {
        BufReader::new(file).read_to_end(&mut buf).expect("read");
    }
    buf
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

    // All samples' reads are resident at once: the shard loop is the OUTER one, so every sample has
    // to be streamable while a given shard is loaded. Reads are small next to a shard.
    let samples: Vec<ShardSample> = pairs
        .into_iter()
        .map(|(r1, r2, out)| {
            let reads1 = read_all(&r1);
            let reads2 = read_all(&r2);
            let n_batches = count_batches(reads1.as_slice(), reads2.as_slice(), BUFFER_SIZE, THREADS);
            ShardSample { r1, r2, out, reads1, reads2, n_batches }
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
                &[u8],
            >(
                &s.reads1, &s.reads2, &shard, s.n_batches, &evidence[i][j], BUFFER_SIZE, THREADS,
                options.args.max_best_flex, MAX_HEADERED,
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

    // Rejoin + align: streams evidence back per batch. The full index is loaded ONCE and reused for
    // every sample -- the second place a batch would otherwise repeat a multi-GB load.
    eprintln!("Loading full index for the align half ...");
    let full: DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD> = DB::load(paths, 1);

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
            &[u8],
        >(
            &s.reads1, &s.reads2, &full, options, out_buffer, &ev, options.args.ranges,
            BUFFER_SIZE, THREADS,
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
