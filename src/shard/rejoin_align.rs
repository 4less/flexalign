//! The rejoin's alignment half: replayed seeds -> PAF (ProjectShard.md §13.8, alignment wiring).
//!
//! [`run_rejoin_seeds`](super::rejoin::run_rejoin_seeds) rebuilds each read's seeds from the shard
//! evidence; this driver takes it the rest of the way. It streams the FASTQ back (with read ids),
//! looks each read's seeds up by ordinal, and runs them through `ModularPE::align_from_seeds` -- the
//! exact reference-comparison half the unsharded pipeline uses -- so the sharded run emits the same
//! PAF the unsharded run would, modulo the §4 seed differences.
//!
//! The database passed here is the *full* index: the align half needs `get_reference` (Hamming
//! extension) and `get_rname` (output), not `get_vrange`. A references-only loader that avoids the
//! resident keys is a later RAM refinement; correctness and speed do not depend on it.


use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::read_fastq_paired_end_state_par_ids;
use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;


use crate::align::modular_workflow::ModularPE;
use crate::align::process::alignment::LIBWFA2Alignment;
use crate::align::process::anchor_extender::StdPairedAnchorExtender;
use crate::align::process::anchor_extractor::StdPairedAnchorExtractor;
use crate::align::process::kmer_extractor::StdKmerExtractor;
use crate::align::process::output::make_output;
use crate::align::process::range_extractor::StdRangeExtractor;
use crate::align::process::seed_extractor::StdSeedExtractor;
use crate::align::stats::Stats;
use crate::database::common::FlexalignDatabase;
use crate::io::output_buffer::OutputBuffer;
use crate::options::Options;

use super::pass::ShardEvidence;
use super::record::GroupRecord;
use super::replay::ReplaySeedExtractor;

/// Rebuilds seeds from `evidence`, then streams the reads and aligns each read's seeds into PAF via
/// the shared `align_from_seeds`. Returns the run's `Stats`. `out_buffer` is where PAF lands
/// (file or stdout); `db` is the full index (references + names).
pub fn run_rejoin_align<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    D: FlexalignDatabase + Clone + Sync + Send,
>(
    evidence: &[ShardEvidence],
    read1: &[u8],
    read2: &[u8],
    db: &D,
    options: &Options,
    out_buffer: OutputBuffer,
    max_ranges: usize,
    buffer_size: usize,
    threads: u32,
) -> Stats {
    let output = make_output(options, out_buffer, db);

    // The front-half stages are unused by `align_from_seeds` but are structural fields of
    // `ModularPE`; constructing them is cheap. The align half uses the anchor extractor/extender,
    // the WFA2 aligner, the database, and the output.
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

    // Per-worker replay state (cloned per thread). Rather than materialising every read's seeds up
    // front, each read gathers its evidence from the shards and rebuilds its seeds on the fly, so
    // only one read's seeds are ever live -- the full seed set never sits in RAM at once.
    let mut replay_fwd = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut replay_rev = ReplaySeedExtractor::<K, C, F>::new(max_ranges);
    let mut groups: Vec<GroupRecord> = Vec::new();
    let mut sort_time = std::time::Duration::ZERO;

    let worker = move |rec_fwd: &RefFastqRecord,
                       rec_rev: &RefFastqRecord,
                       id: ReadId,
                       stats: &mut Stats| {
        groups.clear();
        for shard in evidence.iter() {
            for mate in [0u8, 1] {
                if let Some(ranges) = shard.ranges_for(id.ordinal, mate) {
                    groups.push(GroupRecord { idx: id.ordinal as u32, mate, ranges: ranges.to_vec() });
                }
            }
        }
        replay_fwd.gather(&groups, 0);
        replay_rev.gather(&groups, 1);
        let seeds_fwd = replay_fwd.generate(&mut sort_time);
        let seeds_rev = replay_rev.generate(&mut sort_time);
        // A read with no evidence yields empty seeds -> no anchors -> no PAF line, as unsharded.
        modular_pe.align_from_seeds(rec_fwd, rec_rev, seeds_fwd, seeds_rev, stats);
    };

    read_fastq_paired_end_state_par_ids(
        read1,
        read2,
        buffer_size,
        threads,
        Stats::default(),
        worker,
    )
    .expect("paired fastq read")
}
