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

    /// After the seeds are grouped into anchors, the top x will be extended with the use of hamming distance.
    /// This affects speed negatively but sensitivity and precision positively
    #[arg(short = 'x', long = "extend-top-x", default_value_t = 4)]
    pub extend_top_x: usize,

    /// Extend at most this many anchor pairs (the top Z after the post-pairing sort).
    ///
    /// Extension is the first stage that touches reference sequence, so it is the first that costs
    /// real time; without a bound it runs over EVERY anchor a read produced (15.7 per read measured
    /// on the protal marker DB). Anchors are already sorted by anchor score when they arrive, so
    /// the top Z are the ones worth the Hamming pass. `--extend-top-x` was intended for this and
    /// was never wired up.
    #[arg(short = 'z', long = "extend-top-z", default_value_t = 32)]
    pub extend_top_z: usize,

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

    /// Minimum fraction of the read that must be covered by the alignment (aligned query bases /
    /// read length) for a record to be emitted. `--min-ani` prices only the *aligned* window, so a
    /// short local match to the wrong reference -- with the rest of the read soft-clipped away for
    /// free -- can pass the identity gate; this coverage floor rejects those partial hits. A real
    /// short read maps end-to-end (coverage ~1.0). Set to 0.0 to disable. Default 0.70: on a marker
    /// DB this lifts precision-of-mapped from ~26% to ~92% (matching protal) while keeping per-read
    /// recall above protal's; wrong-reference junk is concentrated below ~0.5 coverage and precision
    /// is flat above the knee, so 0.70 clears it with margin.
    #[arg(long = "min-query-coverage", default_value_t = 0.70)]
    pub min_query_coverage: f64,

    /// Emit SAM instead of PAF. Only aligned records are written -- reads below --min-ani are
    /// omitted entirely rather than emitted as unmapped, so the output is already filtered.
    #[arg(long = "sam", action)]
    pub sam: bool,

    /// Minimum number of ranges for lookup. With max-best-flex defines, none of the ranges might actually yield any seeds.
    #[arg(long = "min-ranges", default_value_t = 4)]
    pub min_ranges: usize,

    /// force_build
    #[arg(long = "force-build", action)]
    pub force_build: bool,

    /// debug
    #[arg(long = "debug", action)]
    pub debug: bool,

    /// align the top y anchors. This happens after anchor extension 
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

    /// Require alignments to be end-to-end on the read, treating soft-clipping as legitimate only
    /// where the read runs off the reference.
    ///
    /// OFF by default, and measured as a net loss until leading soft-clips are emitted: `to_sam`
    /// clips only at the 3' end, so a read hanging off the START of a reference is indistinguishable
    /// from a mid-reference partial hit and gets rejected with it. On the protal marker DB this
    /// removed 0.156 pp of alignments to gain 0.20 pp precision -- 87% of what it rejected was
    /// correct. Revisit once left overhangs are representable.
    #[arg(long = "end-to-end", action)]
    pub end_to_end: bool,

    /// Compute MAPQ from the gapped-alignment scores of the best and second-best candidate,
    /// instead of from their pre-alignment anchor scores.
    ///
    /// The default MAPQ is computed before alignment runs, from `core_matches - mismatches`. Two
    /// candidates that tie as seed chains but separate once aligned therefore get a low MAPQ they
    /// no longer deserve, and vice versa -- and since the anchors are finally *ranked* by aligned
    /// score, the reported alignment may not even be the one the MAPQ describes.
    #[arg(long = "mapq-from-alignment", action)]
    pub mapq_from_alignment: bool,

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

    /// Alias for `--ref-budget auto`.
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

#[derive(Debug)]
pub struct Options {
    pub fwd: Vec<PathBuf>,
    pub rev: Vec<Option<PathBuf>>,
    pub output_prefix: Option<Vec<PathBuf>>,
    pub reference: PathBuf,
    pub reference_database: PathBuf,
    
    pub args: Args,
}


impl Options {
    pub fn from_args(args: Args) -> Self {
        let mut options = Options {
            fwd: vec![PathBuf::default(); 0],
            rev: vec![None; 0],
            reference: PathBuf::default(),
            reference_database: PathBuf::default(),
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
            let s = self.fwd.first().unwrap().to_str().unwrap();
            let s = s.strip_suffix(".gz").unwrap_or(s);
            let s = s.strip_suffix(".bz").unwrap_or(s);
            let s = s.strip_suffix(".bz2").unwrap_or(s);     // Remove .gz if present
            let s = s.rsplit_once('.').map_or(s, |(left, _)| left);
            
            self.output_prefix = Some(vec![PathBuf::from(s); 0]);
        }
        
        if self.output_prefix.is_some() {
            for s in self.output_prefix.as_ref().unwrap() {
                println!("{:?}", s);
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