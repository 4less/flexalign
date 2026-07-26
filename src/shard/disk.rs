//! Disk-streaming shard passes and rejoin (ProjectShard.md §8, §13.10).
//!
//! This is the RAM-bounded path the whole design is for. Evidence is **never accumulated in RAM**:
//! each shard streams over *all* reads and writes its evidence to a per-shard chunk file (one chunk
//! per bioreader batch), then is dropped before the next shard loads. The rejoin re-streams the
//! reads and, per batch, `read_at`s each shard's chunk, decodes just that batch's evidence, and
//! aligns -- so only one batch's evidence is ever resident.
//!
//! Peak RAM is therefore ~ one shard (memory-mapped, touched pages) during a pass, and ~ one
//! batch's evidence during the rejoin -- independent of the read count. This is what lets the
//! sharded run stay below the unsharded run instead of holding every read's ranges at once.

use std::collections::HashMap;
use std::io;
use std::sync::Arc;
use std::time::Duration;

use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::{read_fastq_paired_end_state_par_ids, Merge};
use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

use crate::align::common::{KmerExtractor, NoSAMOutput, Or, RangeExtractor};
use crate::align::modular_workflow::ModularPE;
use crate::align::process::alignment::LIBWFA2Alignment;
use crate::align::process::anchor_extender::StdPairedAnchorExtender;
use crate::align::process::anchor_extractor::StdPairedAnchorExtractor;
use crate::align::process::kmer_extractor::StdKmerExtractor;
use crate::align::process::output::StdPAFOutput;
use crate::align::process::range_extractor::StdRangeExtractor;
use crate::align::process::seed_extractor::StdSeedExtractor;
use crate::align::stats::Stats;
use crate::database::common::FlexalignDatabase;
use crate::io::output_buffer::OutputBuffer;
use crate::options::Options;

use super::chunkfile::{ChunkFileReader, ChunkFileWriter};
use super::emit::{emit_ranges, EmitStats};
use super::record::{GroupReader, GroupRecord, GroupWriter};
use super::replay::ReplaySeedExtractor;

fn add_emit(a: &mut EmitStats, b: &EmitStats) {
    a.headered += b.headered;
    a.headerless += b.headerless;
    a.discarded += b.discarded;
    a.capped += b.capped;
    a.cells += b.cells;
}

/// Counts the batches the reader will cut (= the number of chunks each evidence file needs). A quick
/// pass with no per-read work; batch ids are deterministic, so the count matches the real passes.
pub fn count_batches<R: std::io::Read + Send>(
    read1: R,
    read2: R,
    buffer_size: usize,
    threads: u32,
) -> u64 {
    #[derive(Default, Clone)]
    struct Max(u64);
    impl Merge for Max {
        fn merge_from(&mut self, other: &mut Self) {
            self.0 = self.0.max(other.0);
        }
    }
    read_fastq_paired_end_state_par_ids(
        read1,
        read2,
        buffer_size,
        threads,
        Max(0),
        |_r1: &RefFastqRecord, _r2: &RefFastqRecord, id: ReadId, st: &mut Max| {
            if id.batch_id + 1 > st.0 {
                st.0 = id.batch_id + 1;
            }
        },
    )
    .expect("count batches")
    .0
}

/// Runs one shard over *all* reads, writing its evidence to `evidence_path` as `n_chunks` chunks
/// (one per batch). Nothing is accumulated across batches: each worker holds only the batch it is
/// on and flushes it as a chunk when the batch changes (and its final batch on merge). Returns the
/// emit counters.
pub fn run_shard_pass_to_disk<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    D: FlexalignDatabase + Sync + Clone,
    R: std::io::Read + Send,
>(
    read1: R,
    read2: R,
    db: &D,
    n_chunks: u64,
    evidence_path: &str,
    buffer_size: usize,
    threads: u32,
    max_best_flex: usize,
    max_headered: usize,
) -> io::Result<EmitStats> {
    let writer = Arc::new(ChunkFileWriter::create(evidence_path, n_chunks)?);

    struct PassState {
        chunk: Option<Arc<ChunkFileWriter>>,
        gw: GroupWriter,
        current_batch: Option<u64>,
        emit: EmitStats,
    }
    impl PassState {
        fn flush_current(&mut self) {
            if let Some(bid) = self.current_batch.take() {
                let buf = self.gw.take();
                self.chunk.as_ref().unwrap().write_chunk(bid, &buf).expect("write chunk");
            }
        }
    }
    impl Default for PassState {
        fn default() -> Self {
            Self { chunk: None, gw: GroupWriter::new(), current_batch: None, emit: EmitStats::default() }
        }
    }
    impl Clone for PassState {
        // Per-worker clone shares the writer but starts a fresh, empty batch buffer.
        fn clone(&self) -> Self {
            Self { chunk: self.chunk.clone(), gw: GroupWriter::new(), current_batch: None, emit: EmitStats::default() }
        }
    }
    impl Merge for PassState {
        fn merge_from(&mut self, other: &mut Self) {
            other.flush_current(); // flush the worker's final batch
            add_emit(&mut self.emit, &other.emit);
        }
    }

    let mut kmer_ex = StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default();
    let mut range_ex = StdRangeExtractor::<K, C, F, D>::new(db);
    let mut scratch = Stats::default();

    let init = PassState {
        chunk: Some(Arc::clone(&writer)),
        gw: GroupWriter::new(),
        current_batch: None,
        emit: EmitStats::default(),
    };

    let mut dist_buf: Vec<u8> = Vec::new();
    let worker = move |rec1: &RefFastqRecord, rec2: &RefFastqRecord, id: ReadId, st: &mut PassState| {
        if st.current_batch != Some(id.batch_id) {
            st.flush_current();
            st.current_batch = Some(id.batch_id);
        }
        for (mate, rec) in [(0u8, rec1), (1u8, rec2)] {
            let kmers = kmer_ex.generate(rec, &mut scratch);
            let ranges = range_ex.generate(kmers, rec, &mut scratch);
            let mut records = Vec::new();
            let est = emit_ranges::<F>(ranges, max_best_flex, max_headered, &mut dist_buf, &mut records);
            add_emit(&mut st.emit, &est);
            if !records.is_empty() {
                st.gw
                    .push(&GroupRecord { idx: id.idx_in_batch as u32, mate, ranges: records })
                    .expect("push group (ids must ascend within a batch)");
            }
        }
    };

    let final_state =
        read_fastq_paired_end_state_par_ids(read1, read2, buffer_size, threads, init, worker)
            .expect("paired fastq read");

    let emit = final_state.emit;
    drop(final_state); // releases the extra writer Arc so we can finalize it below

    Arc::try_unwrap(writer)
        .unwrap_or_else(|_| panic!("chunk writer still shared at finish"))
        .finish()?;
    Ok(emit)
}

/// Rejoin over on-disk evidence: re-streams the reads and, per batch, decodes just that batch's
/// evidence from each shard file, rebuilds seeds, and aligns into PAF. Only one batch's evidence is
/// resident at a time. `evidence_paths` are the per-shard chunk files written by
/// [`run_shard_pass_to_disk`]; `db` is the full index (references + names).
pub fn run_disk_rejoin_align<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    D: FlexalignDatabase + Clone + Sync + Send,
    R: std::io::Read + Send,
>(
    read1: R,
    read2: R,
    db: &D,
    options: &Options,
    out_buffer: OutputBuffer,
    evidence_paths: &[String],
    max_ranges: usize,
    buffer_size: usize,
    threads: u32,
) -> io::Result<Stats> {
    let readers: Arc<Vec<ChunkFileReader>> = Arc::new(
        evidence_paths
            .iter()
            .map(|p| ChunkFileReader::open(p))
            .collect::<io::Result<_>>()?,
    );

    let output: Or<StdPAFOutput, NoSAMOutput> =
        Or { a: Some(StdPAFOutput::new(out_buffer)), b: None };

    let mut modular_pe = ModularPE {
        options,
        db,
        kmer_extractor_fwd: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
        kmer_extractor_rev: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
        range_extractor_fwd: StdRangeExtractor::<K, C, F, D>::new(db),
        range_extractor_rev: StdRangeExtractor::<K, C, F, D>::new(db),
        seed_extractor_fwd: StdSeedExtractor::<K, C, F>::new(
            options.args.max_best_flex,
            options.args.max_range_size,
            options.args.min_ranges,
            options.args.ranges,
            options.args.mask_flank_mult, options.args.flank_slack,
        ),
        seed_extractor_rev: StdSeedExtractor::<K, C, F>::new(
            options.args.max_best_flex,
            options.args.max_range_size,
            options.args.min_ranges,
            options.args.ranges,
            options.args.mask_flank_mult, options.args.flank_slack,
        ),
        anchor_extractor: StdPairedAnchorExtractor::new(),
        anchor_extender: StdPairedAnchorExtender::new(db),
        align: LIBWFA2Alignment::default(),
        output,
        rec_fwd_revc: OwnedFastqRecord::new(),
        rec_rev_revc: OwnedFastqRecord::new(),
    };

    let mut replay_fwd = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut replay_rev = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut sort_time = Duration::ZERO;
    // Per-batch decoded evidence (idx_in_batch -> gathered groups across shards); only one batch resident.
    let mut current_batch: Option<u64> = None;
    let mut batch_map: HashMap<u32, Vec<GroupRecord>> = HashMap::new();
    let mut chunk_buf: Vec<u8> = Vec::new();
    let empty: Vec<GroupRecord> = Vec::new();

    let worker = move |rec_fwd: &RefFastqRecord,
                       rec_rev: &RefFastqRecord,
                       id: ReadId,
                       stats: &mut Stats| {
        if current_batch != Some(id.batch_id) {
            batch_map.clear();
            for reader in readers.iter() {
                chunk_buf.clear();
                reader.read_chunk(id.batch_id as usize, &mut chunk_buf).expect("read chunk");
                if chunk_buf.is_empty() {
                    continue;
                }
                let mut gr = GroupReader::new(&chunk_buf).expect("group reader");
                while let Some(g) = gr.next_group().expect("decode group") {
                    batch_map.entry(g.idx).or_default().push(g);
                }
            }
            current_batch = Some(id.batch_id);
        }

        let groups = batch_map.get(&(id.idx_in_batch as u32)).unwrap_or(&empty);
        replay_fwd.gather(groups, 0);
        replay_rev.gather(groups, 1);
        let seeds_fwd = replay_fwd.generate(&mut sort_time);
        let seeds_rev = replay_rev.generate(&mut sort_time);
        modular_pe.align_from_seeds(rec_fwd, rec_rev, seeds_fwd, seeds_rev, stats);
    };

    Ok(
        read_fastq_paired_end_state_par_ids(read1, read2, buffer_size, threads, Stats::default(), worker)
            .expect("paired fastq read"),
    )
}
