use std::{path::PathBuf, str::FromStr};

use clap::Parser;

use crate::utils::infer_output_prefix;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)] // term_width sets it fixed, max term_width can be smaller
pub struct Args {
    /// Forward read of pair (.fastq, .fq)
    #[arg(num_args(0..), short = '1', long, default_values_t = ["".to_string()], action = clap::ArgAction::Append)]
    pub fwd: Vec<String>,

    /// Reverse read of pair (.fastq, .fq)
    #[arg(num_args(0..), short = '2', long, default_values_t = ["".to_string()], action = clap::ArgAction::Append)]
    pub rev: Vec<String>,

    /// Output file
    #[arg(short = '0', long)] // String::default()
    pub output: Option<String>,

    /// Database reference
    #[arg(short = 'r', long = "reference", default_value_t = String::default())]
    pub reference: String,

    /// Input map file 
    #[arg(short, long, default_value_t = String::default())]
    pub map: String,

    /// threads 
    #[arg(short, long, default_value_t = 1)]
    pub threads: u32,

    /// How many minimizers should be looked at
    #[arg(short = 'a', long = "ranges", default_value_t = 15)]
    pub ranges: usize,

    /// For a single minimizer, how many occurances may there be at max.
    #[arg(short = 'b', long = "max-range-size", default_value_t = 256)]
    pub max_range_size: usize,

    /// For all occurrences of a key, flexalign only takes the seeds with the highest matching flanking region.
    /// This limits the number of values to be retrieved in this scenario.
    #[arg(short = 'f', long = "max-best-flex", default_value_t = 16)]
    pub max_best_flex: usize,

    /// Hard ceiling on flank multiplicity. If more than this many flank headers tie at the best
    /// (minimum) Hamming distance to the read's flank, the range emits no seeds — a repetitive
    /// flank context that cannot be localized. Unlike --max-best-flex, this ceiling is NOT relaxed
    /// by the min-ranges recovery pass, so masked ranges stay masked. 0 disables it (default).
    #[arg(long = "mask-flank-mult", default_value_t = 0)]
    pub mask_flank_mult: usize,

    /// Extend at most this many anchor pairs (the top Z after the post-pairing sort).
    ///
    /// Extension is the first stage that touches reference sequence, so it is the first that costs
    /// real time; without a bound it runs over EVERY anchor a read produced (15.7 per read measured
    /// on the protal marker DB). Anchors are already sorted by anchor score when they arrive, so
    /// the top Z are the ones worth the Hamming pass.
    #[arg(short = 'z', long = "extend-top-z", default_value_t = 32)]
    pub extend_top_z: usize,

    /// Ceiling on the tie run past `--extend-top-z`, as a MULTIPLE of it.
    ///
    /// The Z cut is a count applied to a score order, so anchors tying the Zth slot are kept or
    /// dropped by arrival order -- index-lookup order, not evidence. Extension therefore continues
    /// while the next anchor scores exactly equal to the last one taken, letting reference sequence
    /// settle the tie instead of list position.
    ///
    /// 1 disables the tie run (the plain top-Z cut); 2 lets the window double; larger values allow
    /// longer ties. The cap is a SAFETY BELT WITH NO MEASURED COST OR BENEFIT: the concern was that
    /// a conserved marker hit by a whole reference family could tie hundreds deep, and each extra
    /// anchor costs a reference fetch plus a full-read Hamming pass. That pathology has not been
    /// observed -- on protal, cap 64 (effectively uncapped) ran no slower than cap 2 (38.95-39.16 s
    /// vs 39.10-39.27 s, two reps each, warm) and added 4 successful alignments out of ~165,000.
    /// Do not describe the cap as fixing a measured problem; it bounds one that has not appeared.
    ///
    /// Default 1 (off) is PROVISIONAL and rests on incomplete evidence, not on a measured
    /// regression. What is actually measured, cap 2 vs cap 1, same binary, warm, `-t 16`:
    ///
    ///   bacteroides (whole-genome DB)  recall +0.0001 pp, precision +0.0001 pp, 0 record delta
    ///
    /// i.e. a wash on that dataset. The marker-DB arms (protal, markersim) do not yet have a
    /// same-binary paired accuracy measurement -- earlier figures compared against `results/`
    /// baselines produced by an OLDER COMMIT (pre-`c92c8a8`, some at `-t 8`), so they measured that
    /// commit plus this flag. Settle the default with paired runs before quoting any per-dataset
    /// delta here.
    #[arg(long = "extend-tie-cap", default_value_t = 1)]
    pub extend_tie_cap: usize,

    /// align the top y anchors. This happens after anchor extension
    #[arg(short = 'y', long = "align-top-y", default_value_t = 4)]
    pub align_top_y: usize,

    /// Minimum average nucleotide identity, as a fraction, for an alignment to be pursued and
    /// reported. Serves double duty: it is converted into the aligner's abort score, so a candidate
    /// whose alignment cost passes below this identity is dropped mid-alignment rather than
    /// finished, and it gates the output, so no record is emitted below it. The identity is
    /// approximated from the alignment cost (gap cost is priced as mismatch cost), so indel-rich
    /// alignments read slightly lower than their true identity. The default is the usual species
    /// boundary: on a marker database most alignments below it are reads of a species the database
    /// lacks, placed on a relative's marker.
    #[arg(long = "min-ani", default_value_t = 0.95)]
    pub min_ani: f64,

    /// MAPQ subtracted from a mate whose partner did not align to the same reference. Such lone
    /// mates are far less often right than paired ones of the same MAPQ (on protal-mini, 63-83% vs
    /// 93-96% at MAPQ 20-39), so without it a MAPQ cutoff drops good pairs before bad lone mates.
    #[arg(long = "lone-mate-mapq-penalty", default_value_t = 25)]
    pub lone_mate_mapq_penalty: u8,

    /// Alignment cost charged, when the aligned candidates are ranked, for a mate a pair does not
    /// place -- no anchor, or dropped by the aligner. bwa-mem's unpaired penalty plays this role (17
    /// score points, about 34 in this cost space). 0 ranks a one-sided pair by its one mate alone.
    #[arg(long = "unpaired-penalty", default_value_t = 34)]
    pub unpaired_penalty: i32,

    /// A pair whose best and second-best anchor pairs tie on the anchor score (MAPQ 0) takes the
    /// margin of their gapped alignments instead, up to this MAPQ. Lifts reads the alignments tell
    /// apart out of MAPQ 0, at a small cost in precision at MAPQ >= 1 (human: recall at MAPQ >= 1
    /// 96.79 -> 97.08 %, precision 99.88 -> 99.81 %). 0 (default) keeps MAPQ 0.
    #[arg(long = "tie-mapq-cap", default_value_t = 0)]
    pub tie_mapq_cap: u8,

    /// With `--repeat-index`: MAPQ ceiling for a pair with at least `--unresolved-min-cores` cores
    /// the repeat index could not resolve and fewer than `--unresolved-max-seeds` seeds behind its
    /// placement. Such a pair was placed among the copies seeding saw, and the copy it came from
    /// may be one it did not (human: 83 % of these right, against 99-100 % for other pairs).
    #[arg(long = "unresolved-mapq-cap", default_value_t = 3)]
    pub unresolved_mapq_cap: u8,

    /// See `--unresolved-mapq-cap`.
    #[arg(long = "unresolved-min-cores", default_value_t = 6)]
    pub unresolved_min_cores: usize,

    /// See `--unresolved-mapq-cap`.
    #[arg(long = "unresolved-max-seeds", default_value_t = 10)]
    pub unresolved_max_seeds: u32,

    /// Do not write a mate the winning pair has no usable alignment for from its own best aligned
    /// candidate (the lone mate fallback); leave it unmapped.
    #[arg(long = "no-lone-mate-fallback", action)]
    pub no_lone_mate_fallback: bool,

    /// Do not look for a mate the index never placed next to its placed partner (the overhang
    /// rescue: the partner's marker, within the insert size, the mate allowed to run off its end).
    #[arg(long = "no-overhang-rescue", action)]
    pub no_overhang_rescue: bool,

    /// Seed the cores the index drops for having more than --max-range-size positions (young
    /// repeats: LINE-1, satellites) from a repeat index kept next to the index as `<index>.repeat`,
    /// built on the first run that asks for it. A read's core then finds the positions where its
    /// whole K-mer (core and flanks) occurs.
    #[arg(long = "repeat-index", action)]
    pub repeat_index: bool,

    /// Most positions a repeat core's exact K-mer may have to be used as seeds (at most 128);
    /// more leave the core unresolved (`--repeat-index`).
    #[arg(long = "repeat-bucket-max", default_value_t = 64, value_parser = clap::value_parser!(u16).range(1..=128))]
    pub repeat_bucket_max: u16,

    /// Seed reads from regions dense in dropped cores (young repeats) with an FM-index of those
    /// regions, kept next to the index as `<index>.rbwt` and built on the first run that asks for
    /// it. A read's match is extended along the read until it tells the copies apart, so copies
    /// identical over 31 bases but not over the read are resolved. With `--repeat-index`: the BWT
    /// is asked only about the cores the repeat index leaves unresolved.
    #[arg(long = "repeat-bwt", action)]
    pub repeat_bwt: bool,

    /// `--repeat-bwt`: a read is searched when it carries at least this many dropped cores, and the
    /// index holds the regions where a read would (fewer are placed by the read's other cores).
    #[arg(long = "bwt-min-dropped", default_value_t = 8)]
    pub bwt_min_dropped: u32,

    /// `--repeat-bwt`: the read length the indexed regions are sized for (a K-mer's region reaches
    /// this far less K either side). Longer reads still work; their matches stop at a region's end.
    #[arg(long = "bwt-read-len", default_value_t = 150)]
    pub bwt_read_len: u32,

    /// `--repeat-bwt`: a match is turned into positions when it has at most this many occurrences;
    /// more leave its cores unresolved (counted for the MAPQ cap, as with `--repeat-index`).
    #[arg(long = "bwt-max-locate", default_value_t = 8)]
    pub bwt_max_locate: u32,

    /// `--repeat-bwt`: a pair whose mates both align at this identity without the BWT's seeds is
    /// not tried with them (else the result with them is kept when its pair score is better).
    #[arg(long = "bwt-redo-ani", default_value_t = 0.99)]
    pub bwt_redo_ani: f64,

    /// `--repeat-bwt`: the unresolved-repeat MAPQ ceiling (`--unresolved-mapq-cap`) for a pair
    /// written from the BWT's seeds.
    #[arg(long = "bwt-unresolved-cap", default_value_t = 1)]
    pub bwt_unresolved_cap: u8,

    /// `--repeat-bwt`: MAPQ ceiling of a pair written from the BWT's seeds.
    #[arg(long = "bwt-mapq-cap", default_value_t = 60)]
    pub bwt_mapq_cap: u8,

    /// Minimum number of read bases an alignment must actually place for the record to be emitted.
    ///
    /// This is a BASE COUNT, not a fraction. A fraction is the wrong shape on a marker DB: its
    /// denominator has to be "what COULD align", which shrinks as a read hangs further off the end
    /// of a gene -- so the same 40 aligned bases read as 27% on a read sitting inside a gene and
    /// 100% on one dovetailing off it. The gate then varies with geometry rather than with evidence.
    /// A base count does not move: 35 bases is 35 bases of sequence agreement wherever it sits.
    ///
    /// 35 is a little over one syncmer window (k=31) plus slack -- below that a match is short
    /// enough to occur by chance across 14.5 M markers. Set to 0 to disable.
    #[arg(long = "min-query-coverage", default_value_t = 35)]
    pub min_query_coverage: usize,

    /// Rescue a mate the index never reported, by scanning where the pair geometry says it must be.
    ///
    /// Fires only when one mate anchored and the other produced nothing on that reference -- ~1% of
    /// pairs, but the dominant recall loss: 194k of the 224k mates missed on markersim had their
    /// partner emitted. Relaxing `--max-best-flex`/`--ranges` recovers only 16% of them, because the
    /// partner's k-mers are absent or suppressed rather than merely deprioritised; scanning the
    /// reference directly needs no index lookup at all.
    ///
    /// The scan only locates. What it finds is judged by the ordinary alignment and output gates.
    #[arg(long = "mate-rescue", action)]
    pub mate_rescue: bool,

    /// Resolve a pair from its unique k-mers alone when they allow it (SAM output only).
    ///
    /// Every position of a mate's k-mers whose keys have at most `--fast-hit-max-positions`
    /// positions votes for a (reference, diagonal, strand); the winner needs at least
    /// `--fast-hit-min-hits` votes, and the mate is compared base by base against that window.
    /// At most `--fast-hit-max-mismatch` mismatches, the whole read inside the reference, and a
    /// partner resolved the same way on the same reference, the opposite strand and within 1 kb,
    /// and the pair is written as two gapless records without seeding, anchoring or alignment,
    /// with MAPQ from the vote margin over the runner-up locus as the full path scores anchors.
    /// Everything else takes the full path. Modelled on protal 2.0's sparse-first gate.
    /// On by default; `--no-fast-hit` turns it off. Each worker probes it on its first 4096 pairs
    /// and switches it off for the rest of its run when fewer than one in five succeed (a
    /// marker-gene reference, where few reads have any home), so it costs nothing where it does
    /// not pay. (`--fast-hit` is accepted and means the default.)
    #[arg(long = "fast-hit", action)]
    pub fast_hit: bool,

    /// The anchor gate is on by default; this turns it off. An anchor pair is kept only if one
    /// of its anchors has two seeds, an exact full-k-mer seed, or a lone core-only seed whose
    /// flanks mismatch the read in at most `--anchor-gate-flank` bases. Modelled on protal 2.0's
    /// threshold-score gate; the flank rule is what keeps it lossless -- a chance 15-mer match
    /// has ~12 of 16 flank bases wrong, a real read with errors a few. Measured: the reduced
    /// marker DB keeps every correct mate and runs 11% fewer instructions; bacteroides is
    /// identical; the 5 GB mixed set places 8 more correct and 2 more wrong of 4 M mates.
    #[arg(long = "no-anchor-gate", action)]
    pub no_anchor_gate: bool,

    /// Accepted and ignored: the gate is the default now (`--no-anchor-gate` turns it off).
    #[arg(long = "anchor-gate", action, hide = true)]
    pub anchor_gate: bool,

    /// When no anchor pair of a read passes the gate, keep this many of the best weak ones
    /// instead of none (0 drops the read).
    #[arg(long = "anchor-gate-keep", default_value_t = 0)]
    pub anchor_gate_keep: usize,

    /// A lone core-only seed passes the anchor gate when its flanks mismatch the read in at
    /// most this many bases (a chance 15-mer match has ~12 of 16 wrong; a real read with errors
    /// has a few). 16 lets every lone seed through.
    #[arg(long = "anchor-gate-flank", default_value_t = 6)]
    pub anchor_gate_flank: u32,

    /// Map the packed index (compact or hashed) into memory instead of reading it: nothing to
    /// copy at start-up, pages fault in as they are touched. Worth it when the file is in the
    /// page cache and the run is short; the read-in path gets huge pages, this does not.
    #[arg(long = "index-mmap", action)]
    pub index_mmap: bool,

    /// Never take the fast-hit path: every pair goes through seeding, anchoring and alignment.
    #[arg(long = "no-fast-hit", action, conflicts_with = "fast_hit")]
    pub no_fast_hit: bool,

    /// Most mismatches a fast-hit mate may carry (`--fast-hit`).
    #[arg(long = "fast-hit-max-mismatch", default_value_t = 6)]
    pub fast_hit_max_mismatch: usize,

    /// Fewest k-mer hits on the winning diagonal a fast-hit mate needs (`--fast-hit`). 2, not 3:
    /// on bacteroides 3 left 31% of pairs to the full path and 2 leaves 21%, with MORE mates
    /// placed correctly and fewer wrongly (the Hamming check does the real deciding).
    #[arg(long = "fast-hit-min-hits", default_value_t = 2)]
    pub fast_hit_min_hits: usize,

    /// First round of the fast-hit path looks up every N-th k-mer only; the rest are fetched when
    /// that round does not place the pair. 1 looks everything up at once (`--fast-hit`).
    #[arg(long = "fast-hit-stride", default_value_t = 2)]
    pub fast_hit_stride: usize,

    /// Most index positions a k-mer's key may have for the fast-hit path to use it; keys with more
    /// carry flank headers and are left to the full path (`--fast-hit`).
    #[arg(long = "fast-hit-max-positions", default_value_t = 2)]
    pub fast_hit_max_positions: usize,

    /// HARD CAP on rescue attempts per read. Rank is a poor selector on its own -- single-sided
    /// candidates are the normal contents of the ranking tail, so `top 2` already costs 53x `top 1`
    /// for 17% more yield. This only stops a pathological read; `--mate-rescue-score-frac` does the
    /// actual selecting.
    #[arg(long = "mate-rescue-top", default_value_t = 16)]
    pub mate_rescue_top: usize,

    /// Attempt rescue only on candidates scoring at least this fraction of the read's BEST
    /// candidate. Selecting by score rather than rank adapts per read: one clear winner yields one
    /// attempt, several genuinely competitive candidates yield a few, and a long junk tail yields
    /// none of it -- which is what makes looking past rank 1 affordable at all.
    #[arg(long = "mate-rescue-score-frac", default_value_t = 0.85)]
    pub mate_rescue_score_frac: f64,

    /// Minimum `core_matches` (summed seed length) on the ANCHORED side before its position is
    /// trusted as a prior. The whole scan rests on that anchor being real; rescuing the partner of a
    /// weak single-end hit just manufactures a well-formed pair in the wrong place.
    #[arg(long = "mate-rescue-min-core", default_value_t = 30)]
    pub mate_rescue_min_core: usize,

    /// Mismatch margin the best offset must beat the runner-up by. A partner matching about equally
    /// well at several offsets is in a repeat and its position carries no information -- and since a
    /// 1-base shift destroys a Hamming match entirely, a genuine hit normally wins by a mile.
    #[arg(long = "mate-rescue-margin", default_value_t = 8)]
    pub mate_rescue_margin: u32,

    /// Score added to a candidate that has BOTH mates, when ranking candidates.
    ///
    /// Without it the key is just score(mate1)+score(mate2) with a missing side worth zero -- and
    /// since `core_matches` sums overlapping seed lengths it is not bounded by read length, so a
    /// seed-rich repeat hit on ONE mate can outrank a genuine full pair. Two mates agreeing on a
    /// locus at the right insert distance is stronger evidence than one mate scoring well. 0 = off.
    #[arg(long = "proper-pair-bonus", default_value_t = 0)]
    pub proper_pair_bonus: i32,

    /// Pieces the partner is split into for the rescue scan.
    ///
    /// A full-length Hamming is destroyed by one indel -- every base after it mismatches, so the read
    /// exceeds the budget at every offset and is rejected outright. Per segment, an indel can only
    /// ruin the segment containing it; the rest still locate the read. Costs the same byte
    /// comparisons (3 x 50 == 1 x 150). Segments disagreeing about the read start ARE the indel, and
    /// each is carried into the anchor as its own seed so the gapped aligner sees it. 1 = old
    /// full-length behaviour.
    #[arg(long = "mate-rescue-segments", default_value_t = 1)]
    pub mate_rescue_segments: usize,

    /// Base floor for a mate whose PARTNER already aligned cleanly -- the pair-aware counterpart of
    /// `--min-query-coverage`.
    ///
    /// That flag asks whether a mate carries enough sequence to be placed on its own, which is the
    /// wrong question for a pair: 20 bases cannot establish a location, but they can be consistent
    /// with one the partner has already established at the correct insert distance. Judging the mates
    /// independently throws that joint evidence away. Only one mate may lean on the other, and
    /// identity and end-to-end still apply in full -- this relaxes how MUCH sequence is required,
    /// never whether it agrees. Set equal to `--min-query-coverage` to disable.
    #[arg(long = "min-query-coverage-mate", default_value_t = 20)]
    pub min_query_coverage_mate: usize,

    /// Emit PAF instead of SAM. SAM is the default: it carries the CIGAR and NM that the aligner
    /// already computed, so PAF is the lossy option and should be the one you ask for.
    ///
    /// Either way only aligned records are written -- reads below --min-ani are omitted entirely
    /// rather than emitted as unmapped, so the output is already filtered.
    #[arg(long = "paf", action)]
    pub paf: bool,

    /// Minimum number of ranges for lookup. With max-best-flex defines, none of the ranges might actually yield any seeds.
    #[arg(long = "min-ranges", default_value_t = 4)]
    pub min_ranges: usize,

    /// force_build
    #[arg(long = "force-build", action)]
    pub force_build: bool,

    /// debug
    #[arg(long = "debug", action)]
    pub debug: bool,

    /// Verbosity: 0 Off, 1 Error, 2 Warn, 3 Info, 4 Debug, 5 Trace. All logging goes to stderr --
    /// stdout carries the alignments.
    #[arg(long = "log-level", default_value_t = 3)]
    pub log_level: usize,

    /// Output reads that align to the wrong reference sequence. Does not check position of alignment. Requires env variable to be set.
    #[arg(long = "output-prefix-fp-reads")] // String::default()
    pub output_prefix_fp_reads: Option<String>,

    /// Slice the (built) index into N physical shard files + a manifest, then exit. Enables
    /// RAM-bounded shard passes: each shard loads only its slice of the keys array.
    #[arg(long = "shard-slice", value_name = "N")]
    pub shard_slice: Option<usize>,

    /// Dump the candidate anchor list for every read to stderr as `ANCHOR` TSV lines.
    ///
    /// One line per candidate per mate: rank, reference, position, anchor score, aligned score,
    /// whether it reached the aligner. Answers the question a wrong call always raises -- was the
    /// correct reference absent from the candidates, or present and outscored? Those need opposite
    /// fixes (seeding vs scoring), and no amount of output-gate tuning distinguishes them.
    #[arg(long = "dump-anchors", action)]
    pub dump_anchors: bool,

    /// Accept PARTIAL alignments -- a read matching only part way into the middle of a reference.
    ///
    /// ON by default. Partial alignments -- a read matching only part way into the middle of a
    /// reference -- do not explain the read and are rejected.
    ///
    /// It was off historically because it "removed 0.156 pp of alignments to gain 0.20 pp precision,
    /// 87% of what it rejected was correct". That measurement was taken while `to_sam` could only
    /// clip at the 3' end, so a read hanging off the START of a reference was indistinguishable from
    /// a mid-reference partial hit and was rejected with it. Leading soft-clips are now emitted, so
    /// the two cases are separable and the rejection lands only on genuine partials.
    ///
    /// `--allow-partial` restores the old permissive behaviour.
    #[arg(long = "allow-partial", action)]
    pub allow_partial: bool,

    /// Compute MAPQ from the gapped-alignment scores of the best and second-best candidate,
    /// instead of from their pre-alignment anchor scores.
    ///
    /// The default MAPQ is computed before alignment runs, from `core_matches - mismatches`. Two
    /// candidates that tie as seed chains but separate once aligned therefore get a low MAPQ they
    /// no longer deserve, and vice versa -- and since the anchors are finally *ranked* by aligned
    /// score, the reported alignment may not even be the one the MAPQ describes.
    #[arg(long = "mapq-from-alignment", action)]
    pub mapq_from_alignment: bool,

    /// With `--mapq-from-alignment`: MAPQ per mismatch's worth of alignment cost between the best
    /// pair and the runner-up (capped at 60); 0 uses the relative margin instead.
    #[arg(long = "mapq-per-mismatch", default_value_t = 0)]
    pub mapq_per_mismatch: u32,

    /// When the best candidate has one mate and its partner is rescued next to it, rescue the
    /// partners of the runner-ups too and take the MAPQ from the completed pairs.
    #[arg(long = "rescue-mapq", action)]
    pub rescue_mapq: bool,

    /// With `--repeat-bwt`: a pair left under the unresolved-repeat MAPQ cap is resolved as a pair --
    /// every copy of each mate's long exact stretches listed by the BWT, kept where the mates'
    /// copies make proper pairs -- and, when those lists are complete, given the aligned margin as
    /// MAPQ instead of the cap.
    #[arg(long = "pair-resolve", action)]
    pub pair_resolve: bool,

    /// `--pair-resolve`: MAPQ ceiling of a pair resolved on one mate's copies only, some of them left
    /// without a partner copy (about 9 in 10 right on human).
    #[arg(long = "pair-loose-mapq", default_value_t = 20)]
    pub pair_loose_mapq: u8,

    /// `--pair-resolve`: the shortest exact stretch of a mate whose copies are listed.
    #[arg(long = "pair-min-match", default_value_t = 40)]
    pub pair_min_match: u32,

    /// `--pair-resolve`: most copies of a stretch listed; one with more leaves the pair incomplete.
    #[arg(long = "pair-max-copies", default_value_t = 128)]
    pub pair_max_copies: u32,

    /// MAPQ ceiling of a mate whose partner is unmapped or not within 1 kb of it (255: off). On a
    /// genome such mates are mostly wrong at any MAPQ; on a marker database a mate hanging off the
    /// end of a gene often is not, hence off by default.
    #[arg(long = "discordant-mapq-cap", default_value_t = 255)]
    pub discordant_mapq_cap: u8,

    /// Before the aligned candidates are ranked, rescue the partner next to each one-mate candidate
    /// that could still beat the best full pair, and rank it as a pair.
    #[arg(long = "complete-singles", action)]
    pub complete_singles: bool,

    /// MAPQ ceiling of a pair neither of whose mates got a seed from the index or the repeat index
    /// (placed by the repeat BWT or a rescue alone).
    #[arg(long = "seedless-pair-mapq", default_value_t = 60)]
    pub seedless_pair_mapq: u8,

    /// A mate also counts as seedless (for --seedless-pair-mapq / --seedless-mate-mapq) when none of
    /// its index ranges has at most this many positions (0: only a mate without any seed).
    #[arg(long = "seedless-multiplicity", default_value_t = 0)]
    pub seedless_multiplicity: usize,

    /// MAPQ ceiling of a pair with one such mate.
    #[arg(long = "seedless-mate-mapq", default_value_t = 60)]
    pub seedless_mate_mapq: u8,

    /// Spend more on pairs seeding can see are hard (a mate without an index range of at most
    /// --adaptive-multiplicity positions): --adaptive-flex flank ties per core, --adaptive-top-y
    /// candidates aligned.
    #[arg(long = "adaptive", action)]
    pub adaptive: bool,

    #[arg(long = "adaptive-multiplicity", default_value_t = 256)]
    pub adaptive_multiplicity: usize,

    #[arg(long = "adaptive-flex", default_value_t = 128)]
    pub adaptive_flex: usize,

    #[arg(long = "adaptive-top-y", default_value_t = 16)]
    pub adaptive_top_y: usize,

    /// MAPQ ceiling from the pair's core copy counts: this minus --mapq-copies-lo x log2 of the
    /// fewer copies (the mate whose least repeated core has the fewest) minus --mapq-copies-hi x
    /// log2 of the other mate's. 0: off.
    #[arg(long = "mapq-copies", default_value_t = 0)]
    pub mapq_copies: usize,

    #[arg(long = "mapq-copies-lo", default_value_t = 2.0)]
    pub mapq_copies_lo: f64,

    #[arg(long = "mapq-copies-hi", default_value_t = 3.0)]
    pub mapq_copies_hi: f64,

    /// With --adaptive: a mate whose every core was dropped for too many flank ties is seeded
    /// again keeping up to this many ties per core (0: off).
    #[arg(long = "adaptive-ties", default_value_t = 0)]
    pub adaptive_ties: usize,

    /// MAPQ of a pair with no runner-up at another place.
    #[arg(long = "mapq-unique", default_value_t = 30)]
    pub mapq_unique: u32,

    /// Align on N windows of the index, one in memory at a time (RAM-bounded). 1 is the unsharded
    /// run: the whole index in memory, the reads read once.
    ///
    /// Normally unnecessary -- shard bounds are discovered from the reference path
    /// (`<bundle>.s<n>.bounds.json`) and used automatically. Give this when bounds exist for several
    /// counts and you want a specific one, or for a new count (its bounds are computed on first use).
    #[arg(long = "shards", value_name = "N")]
    pub shards: Option<usize>,

    /// Ignore any sliced index and align against the whole one.
    #[arg(long = "no-shards", action)]
    pub no_shards: bool,

    /// Cap the resident reference: `auto`, or a size like `8G` / `512M`. Absent = unbounded.
    ///
    /// Alignment touches ~0.6 MB of distinct reference pages per read, so on a large reference the
    /// whole thing goes resident and sets peak RSS. Unbounded is right when RAM is plentiful --
    /// resident pages cost nothing there. Under a memory limit it is what gets the run OOM-killed,
    /// so with a budget the rejoin drops reference pages once residency exceeds it and re-faults
    /// what it needs again (minor faults; the pages stay in the page cache).
    ///
    /// `auto` reads the cgroup limit if there is one -- container, slurm, `systemd-run -p
    /// MemoryMax=` -- and MemAvailable otherwise, then takes 60% to leave room for reads, evidence
    /// and output buffers.
    #[arg(long = "ref-budget", value_name = "SIZE|auto")]
    pub ref_budget: Option<String>,

    /// Do NOT read the reference sequences up front; page them in as alignment needs them, and
    /// bound what stays resident (implies `--ref-budget auto`).
    ///
    /// By default the WHOLE database -- index and reference -- is read into memory during the load.
    /// `mmap` transfers nothing by itself: the multi-GB read happens later, one 4K fault at a time
    /// in whatever order queries touch it, charged to read processing where it looks like compute
    /// instead of like the load it is. One sequential pass up front moves the same bytes, lets the
    /// kernel read ahead, and makes the reported load figure a real one.
    ///
    /// Use this when RAM is the binding constraint rather than time: alignment touches only the few
    /// hundred bases around each anchor, so most of the reference is never needed, and skipping the
    /// up-front read keeps residency down at the cost of faulting during alignment.
    #[arg(long = "lazy-ref")]
    pub lazy_ref: bool,

    /// Write the per-stage timing block as a machine-readable TSV to this path (one row per
    /// metric: `input\tmetric\tvalue\tunit`), in addition to the human-readable stderr block.
    /// The file is truncated on open; one section is appended per input file processed.
    #[arg(long = "time-log", value_name = "FILE")]
    pub time_log: Option<String>,

    /// Widen flank-seed selection: keep every flank cell whose distance to the read's flank is
    /// within this many mismatches of the *best* (minimum) distance, not only the exact minimum.
    /// 0 (default) reproduces the old behaviour (best-distance ties only); 1 also keeps cells one
    /// mismatch worse, etc. Widening still respects `--max-best-flex` on the resulting set. Applies
    /// to the unsharded path; the sharded wire format carries a single per-range distance, so the
    /// sharded path stays at 0.
    #[arg(long = "flank-slack", value_name = "N", default_value_t = 0)]
    pub flank_slack: u32,
}

impl Args {
    /// True when no reverse mate was given: the single-end pipeline stops at anchors.
    pub fn single_end(&self) -> bool {
        self.rev.iter().all(|r| r.is_empty())
    }

    /// Should this run emit SAM? SAM is the default, but the single-end pipeline never runs gapped
    /// alignment and so has no CIGAR to put in a record -- there PAF is the only possible output,
    /// not a preference. Resolved here so the format is decided once rather than at each writer.
    /// Whether the fast-hit path is on (default yes; `--no-fast-hit` turns it off).
    pub fn fast_hit_on(&self) -> bool {
        !self.no_fast_hit
    }

    pub fn emit_sam(&self) -> bool {
        !self.paf && !self.single_end()
    }
}

#[derive(Debug)]
pub struct Options {
    pub fwd: Vec<PathBuf>,
    pub rev: Vec<Option<PathBuf>>,
    pub output_prefix: Option<Vec<PathBuf>>,
    pub reference: PathBuf,
    
    pub args: Args,
}


impl Options {
    pub fn from_args(args: Args) -> Self {
        let mut options = Options {
            fwd: vec![PathBuf::default(); 0],
            rev: vec![None; 0],
            reference: PathBuf::default(),
            output_prefix: None,
            args: args,
        };
        Self::init(&mut options);
        options
    }

    pub fn init(&mut self) {
        self.fwd.extend(self.args.fwd.iter().map(|x| x.into()));
        // `-2` defaults to one empty string: no reverse mate, the single-end pipeline.
        self.rev.extend(self.args.rev.iter().map(|x| (!x.is_empty()).then(|| x.into())));

        if self.fwd.len() > 1 {
            if self.args.output.is_none() {
                panic!("When processing multiple files in one run you need to provide an output folder for the results to be stored. (--output FOLDER)")
            }

            let inputs = self.fwd.iter().map(|x| x.to_string_lossy().into_owned()).collect::<Vec<String>>();
            self.output_prefix = Some(infer_output_prefix(&inputs)
                .iter()
                .map(|s| { 
                    let mut p = PathBuf::from_str(&self.args.output.as_ref().unwrap()).expect("Cannot turn string into path");
                    p.push(s);
                    p
                })
                .collect::<Vec<_>>());
        } else if self.fwd.len() == 1 && self.args.output.is_some() {
            // ONE input pair: `--output` names the output FILE itself, as the README and the
            // sharded path both specify ("a file for a single input pair, a directory for several").
            //
            // This used to derive a name from the INPUT path and then build it with `vec![p; 0]` --
            // a zero-length vector -- so the per-input lookup in process_fastq found nothing and
            // panicked with "There is no output for input", and `--output` was ignored either way.
            self.output_prefix =
                Some(vec![PathBuf::from(self.args.output.as_ref().unwrap())]);
        }

        if self.output_prefix.is_some() {
            for s in self.output_prefix.as_ref().unwrap() {
                // STDERR, not stdout: stdout carries the alignments (PAF/SAM), so anything else
                // written there corrupts the output stream.
                log::debug!("output -> {:?}", s);
            }
        }

        self.reference.push(self.args.reference.clone());
    }
}

pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

/// The k-mer index a run builds, converts to and loads. Runs always use `Hashed`; the other
/// formats remain for the examples and tools that read older indexes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum IndexFormat {
    /// Packed values, compact key table: the least memory.
    Compact,
    /// Packed values, line-hashed key table: one cache line a lookup, ~10 bytes a used key more.
    Hashed,
    /// The Flexmap blob (memory-mapped).
    Blob,
}
