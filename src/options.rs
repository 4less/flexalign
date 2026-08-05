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
    /// alignments read slightly lower than their true identity.
    #[arg(long = "min-ani", default_value_t = 0.85)]
    pub min_ani: f64,

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

    /// Use the N-shard slicing of the index (RAM-bounded: one shard resident at a time).
    ///
    /// Normally unnecessary -- a sliced index is discovered from the reference path
    /// (`<index>.s<n>.shards.json`) and used automatically. Give this when several slicings exist
    /// and you want a specific one, or to slice a new shard count on the fly.
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
        self.rev.extend(self.args.rev.iter().map(|x| Some(x.into())));

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