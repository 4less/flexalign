use std::{fmt::Display, time::Duration};
use rgb::RGB8;
use bioreader::parallel::fastq::Merge;
use textplots::{Chart, ColorPlot, Shape};

use crate::GOLDSTD_EVAL;

use super::common::{AlignStage, AlignTrace};
use super::eval::MapqEvaluation; 
#[derive(Clone, Debug)]
pub struct Stats {
    pub reads_processed: usize,
    pub kmers_processed: usize,
    pub minimizer: usize,
    pub ranges: usize,
    pub seeds: usize,
    pub anchors: usize,
    pub alignments: usize,
    pub alignments_successful: usize,
    pub alignments_partial: usize,
    pub alignments_dropped: usize,
    pub retrieved_ranges: usize,
    /// Total flank-header elements scanned during seed extraction (sum of N over every headered or
    /// discarded range walked, including the recovery pass). Lets us report seed-extraction time
    /// per range and per header element.
    pub header_elements: usize,

    // --debug diagnostics: for simulated reads whose true reference is encoded in the header,
    // whether the true anchor is present among all anchors / is the best anchor / the best anchor
    // is the true one after alignment.
    pub dbg_checked: usize,
    pub dbg_true_anchor_present: usize,
    pub dbg_true_anchor_best: usize,

    /// Records dropped by the `--min-ani` output filter, and records that passed the filter but
    /// could not be written as SAM because the anchor carried no usable CIGAR / reference span
    /// (an anchor that went through `whole_align`, which does not set `reference_cigar_range`).
    /// Both are silent losses otherwise, so they are counted here.
    pub filtered_below_min_ani: usize,
    pub sam_records_skipped: usize,

    pub time_get_kmers: Duration,
    pub time_get_minimizer: Duration,
    pub time_get_vranges: Duration,
    pub time_get_ranges: Duration,
    pub time_range_sorting: Duration,
    pub time_seed_sorting: Duration,
    pub time_anchor_sorting: Duration,
    pub time_reverse_complement: Duration,
    pub time_hamming_screen: Duration,
    pub time_get_anchors: Duration,
    pub time_range_header: Duration,
    pub time_offset: Duration,
    pub time_checking_anchors: Duration,
    pub time_alignment: Duration,

    /// Per-stage timings *inside* the alignment, and the histogram of how far anchors got before
    /// they were dropped. Both are written only under [`STAGE_TRACE`](crate::STAGE_TRACE); see
    /// [`AlignTrace`]. `align_stage_reached` is indexed by [`AlignStage`] and counts every aligned
    /// anchor; `align_stage_reached_true` counts only anchors on the read's true reference, so the
    /// two together say which filter discards true alignments rather than merely how many die.
    pub trace: AlignTrace,
    pub align_stage_reached: [usize; 6],
    pub align_stage_reached_true: [usize; 6],

    /// Truth-survival funnel outside the alignment itself. Gated on `GOLDSTD_EVAL`, so zero cost
    /// when off. Each counts reads whose true anchor was still alive at that point.
    pub funnel_true_in_screen: usize,
    /// Same test over the window extension ACTUALLY ran (`screen_z` plus any tie run). Paired with
    /// `funnel_true_in_screen` on purpose: the screen counter is fixed-width so it stays comparable
    /// across runs, but with `--extend-tie-cap > 1` it no longer describes what was screened, and
    /// alone it would report "truth lost at the screen" for reads whose truth the tie run rescued.
    pub funnel_true_in_extended: usize,
    pub funnel_true_in_align_window: usize,
    pub funnel_true_reported: usize,
    /// Same bottom-of-funnel test as `funnel_true_reported`, but counting MATES rather than pairs.
    ///
    /// The pair counter marks a pair reported when EITHER mate carries the truth, so it cannot see a
    /// pair that emitted one mate and dropped the other -- which measurement showed to be the
    /// dominant recall loss (194k of 224k missing mates on markersim had their partner emitted).
    /// Against a denominator of `2 * dbg_checked`, this counter does see it.
    pub funnel_true_reported_mates: usize,
    /// The rest of the funnel, per MATE. The pair-level counters above are satisfied when EITHER
    /// mate carries the truth, so all of them -- not just `reported` -- are blind to a pair that
    /// keeps one mate and loses the other. Only these say at WHICH stage a mate disappeared.
    pub funnel_true_anchor_mates: usize,
    pub funnel_true_in_screen_mates: usize,
    pub funnel_true_in_align_window_mates: usize,

    /// Shape of what was emitted per pair, counted on EVERY run (not just gold builds): a pair that
    /// emits one mate and drops the other is a silent half of a placement, invisible in the mate
    /// totals and in the pair-level funnel alike.
    pub pairs_half_emitted: usize,
    /// A half-pair split by WHY its dropped mate is missing (gold builds only): the mate either had
    /// no anchor on the true reference at all -- seeding never offered it, which a scan-based rescue
    /// can still reach -- or it had one and the ranking chose a pair without it, which rescue cannot
    /// help with and which needs the selection fixed instead.
    pub halfpair_mate_never_anchored: usize,
    pub halfpair_mate_lost_selection: usize,
    pub pairs_none_emitted: usize,

    /// Where a half/empty pair is BORN, in the paired anchor extractor. All three describe a
    /// reference on which BOTH mates had seeds -- so the pairing had the material for a full pair
    /// and produced less than one:
    ///   `empty_side`   one mate's seeds grouped into zero anchors; the other is paired with None.
    ///   `multi_empty`  same, but the surviving side has SEVERAL anchors -- the cross-product loop
    ///                  then iterates zero times and discards every one of them.
    ///   `insert_size`  both sides had anchors, but no fwd/rev combination came within the insert
    ///                  size cutoff, and there is no fallback -- so both mates are dropped.
    /// Mate rescue (`--mate-rescue`): pairs where one side anchored and the other did not, so a scan
    /// was attempted, and how many of those located a partner within the mismatch budget. `attempted`
    /// tracks the opportunity, `succeeded` the yield -- the gap between them is how often the pair
    /// geometry was right but the sequence was not there.
    /// Mates emitted only because their partner aligned cleanly -- they cleared
    /// `--min-query-coverage-mate` but not `--min-query-coverage`. Isolates what the pair-aware
    /// floor is actually contributing, so it can be judged rather than assumed.
    pub mate_backed_emitted: usize,

    pub rescue_attempted: usize,
    pub rescue_succeeded: usize,
    /// Why a rescue did not happen or was refused. `skipped_weak` is the anchored side failing the
    /// quality gate (its position is not worth trusting); `rejected_ambiguous` is the scan finding a
    /// partner that matched about equally well elsewhere in the window, i.e. a repeat.
    pub rescue_skipped_weak: usize,
    pub rescue_rejected_ambiguous: usize,
    /// Rescues whose segments disagreed about the read start -- i.e. the scan located the partner AND
    /// detected an indel in it. Previously these were simply rejected.
    pub rescue_indel_suspected: usize,
    /// Reads where rescue fired AND a two-sided candidate existed somewhere in the list. Settles
    /// whether rescue can pre-empt a full pair that scoring would have chosen anyway -- the pair
    /// score sums both mates, so this should be rare, but it was asserted once without evidence.
    pub rescue_with_full_pair_available: usize,

    pub pair_leak_empty_side: usize,
    pub pair_leak_multi_empty: usize,
    pub pair_leak_insert_size: usize,

    pub threads: usize,

    pub gold_std_evaluation: Option<MapqEvaluation>,
}

pub trait EDisplay {
    fn edisplay(&mut self);
}

impl<'a> EDisplay for Chart<'a> {
    fn edisplay(&mut self) {
        self.axis();
        self.figures();

        eprintln!("{}", self);
    }
}

impl Stats {
    /// Per-thread divisor used by the human-readable `Display` (durations are summed across threads
    /// during `Merge`, so the reported per-stage time is the sum divided by the thread count).
    fn thread_divisor(&self) -> u32 {
        (self.threads as u32).max(1)
    }

    /// Truncate `path` and write the `--time-log` TSV header. Call once before processing; each
    /// input then appends a section via [`append_time_log`](Self::append_time_log).
    pub fn init_time_log(path: &str) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "input\tmetric\tvalue\tunit")
    }

    /// Append this input's timing block to a TSV `--time-log` file as tidy rows
    /// (`input\tmetric\tvalue\tunit`), so a multi-input run yields one file with one section per
    /// input. Times are the same per-thread-averaged seconds shown in the stderr block; counts and
    /// per-read rates are also emitted.
    pub fn append_time_log(&self, path: &str, input: &str) -> std::io::Result<()> {
        use std::io::Write;
        let td = self.thread_divisor();
        let secs = |d: Duration| (d / td).as_secs_f64();
        let reads = self.reads_processed.max(1) as f64;

        // (metric, value, unit). Order mirrors the stderr block, then counts, then per-read rates.
        let rows: Vec<(&str, f64, &str)> = vec![
            ("reverse_complement", secs(self.time_reverse_complement), "s"),
            ("kmers",              secs(self.time_get_kmers), "s"),
            ("minimizers",         secs(self.time_get_minimizer), "s"),
            ("ranges",             secs(self.time_get_ranges), "s"),
            ("vranges",            secs(self.time_get_vranges), "s"),
            ("range_sorting",      secs(self.time_range_sorting), "s"),
            ("range_headers",      secs(self.time_range_header), "s"),
            ("seed_sorting",       secs(self.time_seed_sorting), "s"),
            ("anchors",            secs(self.time_get_anchors), "s"),
            ("anchor_sorting",     secs(self.time_anchor_sorting), "s"),
            ("hamming_screen",     secs(self.time_hamming_screen), "s"),
            ("offsets",            secs(self.time_offset), "s"),
            ("checking_anchors",   secs(self.time_checking_anchors), "s"),
            ("alignment",          secs(self.time_alignment), "s"),
            // Sub-stages of `alignment`, from AlignTrace. All zero unless STAGE_TRACE is compiled
            // in; they are components of the line above and must not be added to it.
            ("align_seed_extension",   secs(self.trace.t_seed_extension), "s"),
            ("align_gapless_prefilter", secs(self.trace.t_prefilter), "s"),
            ("align_left_flank",       secs(self.trace.t_left_flank), "s"),
            ("align_interior",         secs(self.trace.t_middle), "s"),
            ("align_right_flank",      secs(self.trace.t_right_flank), "s"),
            ("header_elements",         self.header_elements as f64, "count"),
            // Total (thread-summed) seed-extraction time over total counts -- CPU time per unit.
            ("seedext_s_per_range",     self.time_range_header.as_secs_f64() / self.ranges.max(1) as f64, "s"),
            ("seedext_s_per_header",    self.time_range_header.as_secs_f64() / self.header_elements.max(1) as f64, "s"),
            ("reads",                   self.reads_processed as f64, "count"),
            ("alignments",              self.alignments as f64, "count"),
            ("alignments_successful",   self.alignments_successful as f64, "count"),
            ("alignments_partial",      self.alignments_partial as f64, "count"),
            ("alignments_dropped",      self.alignments_dropped as f64, "count"),
            ("filtered_below_min_ani",  self.filtered_below_min_ani as f64, "count"),
            ("sam_records_skipped",     self.sam_records_skipped as f64, "count"),
            ("threads",                 self.threads as f64, "count"),
            // Truth-survival funnel. Zero unless GOLDSTD_EVAL; `funnel_checked` is the denominator
            // (reads whose true reference could be read off the header), and each later row is how
            // many of those still had their true anchor alive at that point.
            ("funnel_checked",              self.dbg_checked as f64, "count"),
            ("funnel_anchor_formed",        self.dbg_true_anchor_present as f64, "count"),
            ("funnel_in_screen",            self.funnel_true_in_screen as f64, "count"),
            ("funnel_in_extended",          self.funnel_true_in_extended as f64, "count"),
            ("funnel_screen_ranked_best",   self.dbg_true_anchor_best as f64, "count"),
            ("funnel_in_align_window",      self.funnel_true_in_align_window as f64, "count"),
            ("funnel_reported",             self.funnel_true_reported as f64, "count"),
            ("funnel_reported_mates",       self.funnel_true_reported_mates as f64, "count"),
            ("funnel_anchor_mates",         self.funnel_true_anchor_mates as f64, "count"),
            ("funnel_in_screen_mates",      self.funnel_true_in_screen_mates as f64, "count"),
            ("funnel_in_align_window_mates", self.funnel_true_in_align_window_mates as f64, "count"),
            ("pairs_half_emitted",          self.pairs_half_emitted as f64, "count"),
            ("halfpair_mate_never_anchored", self.halfpair_mate_never_anchored as f64, "count"),
            ("halfpair_mate_lost_selection", self.halfpair_mate_lost_selection as f64, "count"),
            ("pairs_none_emitted",          self.pairs_none_emitted as f64, "count"),
            ("mate_backed_emitted",         self.mate_backed_emitted as f64, "count"),
            ("rescue_attempted",            self.rescue_attempted as f64, "count"),
            ("rescue_succeeded",            self.rescue_succeeded as f64, "count"),
            ("rescue_skipped_weak",         self.rescue_skipped_weak as f64, "count"),
            ("rescue_rejected_ambiguous",   self.rescue_rejected_ambiguous as f64, "count"),
            ("rescue_indel_suspected",      self.rescue_indel_suspected as f64, "count"),
            ("rescue_with_full_pair_available", self.rescue_with_full_pair_available as f64, "count"),
            ("pair_leak_empty_side",        self.pair_leak_empty_side as f64, "count"),
            ("pair_leak_multi_empty",       self.pair_leak_multi_empty as f64, "count"),
            ("pair_leak_insert_size",       self.pair_leak_insert_size as f64, "count"),
            ("minimizers_per_read",     self.minimizer as f64 / reads, "rate"),
            ("ranges_per_read",         self.ranges as f64 / reads, "rate"),
            ("retrieved_ranges_per_read", self.retrieved_ranges as f64 / reads, "rate"),
            ("seeds_per_read",          self.seeds as f64 / reads, "rate"),
            ("anchors_per_read",        self.anchors as f64 / reads, "rate"),
        ];

        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;

        // How far anchors got before being dropped, per AlignStage: once over every aligned anchor,
        // once over anchors on the read's true reference. Written directly rather than pushed into
        // `rows`, whose metric names are `&'static str`. A stage that kills many anchors and no true
        // ones is a good filter; one that kills few but takes true alignments with them is not, and
        // that distinction is the whole point of counting both.
        for (i, st) in AlignStage::ALL.iter().enumerate() {
            let label = st.label();
            writeln!(f, "{input}\talign_stage_{label}\t{}\tcount", self.align_stage_reached[i])?;
            writeln!(f, "{input}\talign_stage_true_{label}\t{}\tcount", self.align_stage_reached_true[i])?;
        }

        for (metric, value, unit) in rows {
            if unit == "s" || unit == "rate" {
                writeln!(f, "{input}\t{metric}\t{value:.9}\t{unit}")?;
            } else {
                writeln!(f, "{input}\t{metric}\t{value:.0}\t{unit}")?;
            }
        }
        Ok(())
    }

    pub fn plot_mapq(&self) {
        if self.gold_std_evaluation.is_none() { return };
        
        let gse = self.gold_std_evaluation.as_ref().unwrap();
        Chart::new(300, 120, 0f32, 100f32)
            .linecolorplot(&Shape::Continuous(Box::new(|x| {
                let be = gse.binary_evaluator(x as usize);
                be.true_positive_rate() as f32
            })), 
            RGB8 {
                r: 255_u8,
                g: 0,
                b: 0,
            },)
            .linecolorplot(&Shape::Continuous(Box::new(|x| {
                let be: super::eval::BinaryEvaluator = gse.binary_evaluator(x as usize);
                be.false_positive_rate() as f32
            })), 
            RGB8 {
                r: 255_u8,
                g: 255_u8,
                b: 0,
            },).edisplay();
    }
}

impl Merge for Stats {
    fn merge_from(&mut self, other: &mut Self) {
        self.reads_processed += other.reads_processed;
        self.kmers_processed += other.kmers_processed;
        self.minimizer += other.minimizer;

        self.time_reverse_complement += other.time_reverse_complement;
        self.time_hamming_screen += other.time_hamming_screen;

        self.time_range_sorting += other.time_range_sorting;
        self.time_seed_sorting += other.time_seed_sorting;
        self.time_anchor_sorting += other.time_anchor_sorting;

        self.time_get_kmers += other.time_get_kmers;
        self.time_get_minimizer += other.time_get_minimizer;
        self.time_get_ranges += other.time_get_ranges;
        self.time_get_vranges += other.time_get_vranges;
        self.time_range_header += other.time_range_header;
        self.time_get_anchors += other.time_get_anchors;
        self.time_offset += other.time_offset;
        self.time_checking_anchors += other.time_checking_anchors;
        self.time_alignment += other.time_alignment;

        self.ranges += other.ranges;
        self.retrieved_ranges += other.retrieved_ranges;
        self.header_elements += other.header_elements;
        self.dbg_checked += other.dbg_checked;
        self.dbg_true_anchor_present += other.dbg_true_anchor_present;
        self.dbg_true_anchor_best += other.dbg_true_anchor_best;
        self.filtered_below_min_ani += other.filtered_below_min_ani;
        self.trace.merge_from(&other.trace);
        for i in 0..self.align_stage_reached.len() {
            self.align_stage_reached[i] += other.align_stage_reached[i];
            self.align_stage_reached_true[i] += other.align_stage_reached_true[i];
        }
        self.funnel_true_in_screen += other.funnel_true_in_screen;
        self.funnel_true_in_extended += other.funnel_true_in_extended;
        self.funnel_true_in_align_window += other.funnel_true_in_align_window;
        self.funnel_true_reported += other.funnel_true_reported;
        self.funnel_true_reported_mates += other.funnel_true_reported_mates;
        self.funnel_true_anchor_mates += other.funnel_true_anchor_mates;
        self.funnel_true_in_screen_mates += other.funnel_true_in_screen_mates;
        self.funnel_true_in_align_window_mates += other.funnel_true_in_align_window_mates;
        self.pairs_half_emitted += other.pairs_half_emitted;
        self.halfpair_mate_never_anchored += other.halfpair_mate_never_anchored;
        self.halfpair_mate_lost_selection += other.halfpair_mate_lost_selection;
        self.pairs_none_emitted += other.pairs_none_emitted;
        self.mate_backed_emitted += other.mate_backed_emitted;
        self.rescue_attempted += other.rescue_attempted;
        self.rescue_succeeded += other.rescue_succeeded;
        self.rescue_skipped_weak += other.rescue_skipped_weak;
        self.rescue_rejected_ambiguous += other.rescue_rejected_ambiguous;
        self.rescue_indel_suspected += other.rescue_indel_suspected;
        self.rescue_with_full_pair_available += other.rescue_with_full_pair_available;
        self.pair_leak_empty_side += other.pair_leak_empty_side;
        self.pair_leak_multi_empty += other.pair_leak_multi_empty;
        self.pair_leak_insert_size += other.pair_leak_insert_size;
        self.sam_records_skipped += other.sam_records_skipped;
        self.seeds += other.seeds;
        self.anchors += other.anchors;
        self.alignments += other.alignments;
        self.alignments_successful += other.alignments_successful;
        self.alignments_partial += other.alignments_partial;
        self.alignments_dropped += other.alignments_dropped;
        self.threads += 1;

        if self.gold_std_evaluation.is_some() && other.gold_std_evaluation.is_some() {
            self.gold_std_evaluation.as_mut().unwrap().merge_from(&mut other.gold_std_evaluation.as_mut().unwrap());
        }
    }
}

impl Display for Stats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, 
            "\
            Time for getting reverse complement.........{:?}\n\
            Time for getting kmers......................{:?}\n\
            ....Time for getting minimizers.............{:?}\n\
            Time for getting ranges.....................{:?}\n\
            ....Time for getting vranges................{:?}\n\
            Time for sorting ranges.....................{:?}\n\
            Time for getting range headers..............{:?}\n\
            ....Header elements scanned.................{}\n\
            ....Seed-extraction time per range..........{:?}\n\
            ....Seed-extraction time per header elem....{:?}\n\
            Time for sorting seeds......................{:?}\n\
            Time for getting anchors....................{:?}\n\
            Time for sorting anchors....................{:?}\n\
            Time for the hamming screen.................{:?}\n\
            Time for calculating offsets................{:?}\n\
            Time for checking anchors...................{:?}\n\
            Time for alignment..........................{:?}\n\
            ....Time for seed extension.................{:?}\n\
            ....Time for gapless prefilter..............{:?}\n\
            ....Time for left flank.....................{:?}\n\
            ....Time for interior.......................{:?}\n\
            ....Time for right flank....................{:?}\n\
            Anchor stops gapless prefilter..............{} all, {} true\n\
            Anchor stops left flank.....................{} all, {} true\n\
            Anchor stops interior.......................{} all, {} true\n\
            Anchor stops right flank....................{} all, {} true\n\
            Anchor stops complete.......................{} all, {} true\n\
            Truth funnel reads checked..................{}\n\
            ....Truth anchor formed.....................{}\n\
            ....Truth in hamming screen.................{}\n\
            ....Truth ranked best after screen..........{}\n\
            ....Truth in alignment window...............{}\n\
            ....Truth reported..........................{}\n\n\
            MATE-level funnel (denominator 2x checked)\n\
            ....Truth anchor formed (mates).............{}\n\
            ....Truth in hamming screen (mates).........{}\n\
            ....Truth in alignment window (mates).......{}\n\
            ....Truth reported (mates)..................{}\n\n\
            Pairs emitting only ONE mate................{:?}\n\
            ....dropped mate never anchored (seeding)...{:?}\n\
            ....dropped mate lost ranking (selection)...{:?}\n\
            Pairs emitting NEITHER mate.................{:?}\n\
            Mates emitted on partner evidence...........{:?}\n\
            Mate rescue attempted.......................{:?}\n\
            Mate rescue succeeded.......................{:?}\n\
            ....skipped, anchored side too weak.........{:?}\n\
            ....rejected, ambiguous in window...........{:?}\n\
            ....located WITH an indel (segments split)..{:?}\n\
            ....fired though a full pair existed........{:?}\n\
            ....pairing: one side grouped to 0 anchors..{:?}\n\
            ....pairing: ^ and other side had several...{:?}\n\
            ....pairing: no pair within insert size.....{:?}\n\n\
            Total Reads.................................{:?}\n\
            Total Alignments............................{:?}\n\
            Total Alignments successful.................{:?}\n\
            Total Alignments partial....................{:?}\n\
            Total Alignments dropped....................{:?}\n\
            Filtered below min-ani/gate.................{:?}\n\
            SAM records skipped (no CIGAR)..............{:?}\n\
            Total Minimizers per read...................{:.2}x\n\
            Total Ranges per read.......................{:.2}x\n\
            Total Retrieved Ranges per read.............{:.2}x\n\
            Total Seeds per read........................{:.2}x\n\
            Total Anchors per read......................{:.2}x\n\
            Total Alignments per read...................{:.2}x\n\
            Total Alignments success per read...........{:.2}x\n\
            Total Alignments partial per read...........{:.2}x\n\
            Total Alignments dropped per read...........{:.2}x\
            {}",
            self.time_reverse_complement / self.threads as u32,
            self.time_get_kmers / self.threads as u32,
            self.time_get_minimizer / self.threads as u32,
            self.time_get_ranges / self.threads as u32,
            self.time_get_vranges / self.threads as u32,
            self.time_range_sorting / self.threads as u32,
            self.time_range_header / self.threads as u32,
            self.header_elements,
            // Per-range / per-header: total (thread-summed) seed-extraction time over total counts.
            self.time_range_header / (self.ranges.min(u32::MAX as usize).max(1) as u32),
            self.time_range_header / (self.header_elements.min(u32::MAX as usize).max(1) as u32),
            self.time_seed_sorting / self.threads as u32,
            self.time_get_anchors / self.threads as u32,
            self.time_anchor_sorting / self.threads as u32,
            self.time_hamming_screen / self.threads as u32,
            self.time_offset / self.threads as u32,
            self.time_checking_anchors / self.threads as u32,
            self.time_alignment / self.threads as u32,
            self.trace.t_seed_extension / self.threads as u32,
            self.trace.t_prefilter / self.threads as u32,
            self.trace.t_left_flank / self.threads as u32,
            self.trace.t_middle / self.threads as u32,
            self.trace.t_right_flank / self.threads as u32,
            self.align_stage_reached[AlignStage::Prefilter as usize],
            self.align_stage_reached_true[AlignStage::Prefilter as usize],
            self.align_stage_reached[AlignStage::LeftFlank as usize],
            self.align_stage_reached_true[AlignStage::LeftFlank as usize],
            self.align_stage_reached[AlignStage::Middle as usize],
            self.align_stage_reached_true[AlignStage::Middle as usize],
            self.align_stage_reached[AlignStage::RightFlank as usize],
            self.align_stage_reached_true[AlignStage::RightFlank as usize],
            self.align_stage_reached[AlignStage::Complete as usize],
            self.align_stage_reached_true[AlignStage::Complete as usize],
            self.dbg_checked,
            self.dbg_true_anchor_present,
            self.funnel_true_in_screen,
            self.dbg_true_anchor_best,
            self.funnel_true_in_align_window,
            self.funnel_true_reported,
            self.funnel_true_anchor_mates,
            self.funnel_true_in_screen_mates,
            self.funnel_true_in_align_window_mates,
            self.funnel_true_reported_mates,
            self.pairs_half_emitted,
            self.halfpair_mate_never_anchored,
            self.halfpair_mate_lost_selection,
            self.pairs_none_emitted,
            self.mate_backed_emitted,
            self.rescue_attempted,
            self.rescue_succeeded,
            self.rescue_skipped_weak,
            self.rescue_rejected_ambiguous,
            self.rescue_indel_suspected,
            self.rescue_with_full_pair_available,
            self.pair_leak_empty_side,
            self.pair_leak_multi_empty,
            self.pair_leak_insert_size,
            self.reads_processed,
            self.alignments,
            self.alignments_successful,
            self.alignments_partial,
            self.alignments_dropped,
            self.filtered_below_min_ani,
            self.sam_records_skipped,
            self.minimizer as f64 / self.reads_processed as f64,
            self.ranges as f64 / self.reads_processed as f64,
            self.retrieved_ranges as f64 / self.reads_processed as f64,
            self.seeds as f64 / self.reads_processed as f64,
            self.anchors as f64 / self.reads_processed as f64,
            self.alignments as f64 / self.reads_processed as f64,
            self.alignments_successful as f64 / self.reads_processed as f64,
            self.alignments_partial as f64 / self.reads_processed as f64,
            self.alignments_dropped as f64 / self.reads_processed as f64,
            if self.gold_std_evaluation.is_some() {
                "\n\n".to_string() + &self.gold_std_evaluation.as_ref().unwrap().to_string()
            } else { 
                "".to_string() 
            })
    }
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            reads_processed: 0,
            kmers_processed: 0,
            minimizer: 0,
            ranges: 0,
            seeds: 0,
            anchors: 0,
            alignments: 0,
            alignments_successful: 0,
            alignments_partial: 0,
            alignments_dropped: 0,
            retrieved_ranges: 0,
            header_elements: 0,
            dbg_checked: 0,
            dbg_true_anchor_present: 0,
            dbg_true_anchor_best: 0,
            filtered_below_min_ani: 0,
            sam_records_skipped: 0,

            time_reverse_complement: Duration::default(),
            time_hamming_screen: Duration::default(),
            time_get_kmers: Duration::default(),
            time_get_minimizer: Duration::default(),
            time_get_ranges: Duration::default(),
            time_get_vranges: Duration::default(),
            time_range_sorting: Duration::default(),
            time_seed_sorting: Duration::default(),
            time_anchor_sorting: Duration::default(),
            time_range_header: Duration::default(),
            time_offset: Duration::default(),
            time_checking_anchors: Duration::default(),
            time_get_anchors: Duration::default(),
            time_alignment: Duration::default(),
            trace: AlignTrace::default(),
            align_stage_reached: [0; 6],
            align_stage_reached_true: [0; 6],
            funnel_true_in_screen: 0,
            funnel_true_in_extended: 0,
            funnel_true_in_align_window: 0,
            funnel_true_reported: 0,
            funnel_true_reported_mates: 0,
            funnel_true_anchor_mates: 0,
            funnel_true_in_screen_mates: 0,
            funnel_true_in_align_window_mates: 0,
            pairs_half_emitted: 0,
            halfpair_mate_never_anchored: 0,
            halfpair_mate_lost_selection: 0,
            pairs_none_emitted: 0,
            mate_backed_emitted: 0,
            rescue_attempted: 0,
            rescue_succeeded: 0,
            rescue_skipped_weak: 0,
            rescue_rejected_ambiguous: 0,
            rescue_indel_suspected: 0,
            rescue_with_full_pair_available: 0,
            pair_leak_empty_side: 0,
            pair_leak_multi_empty: 0,
            pair_leak_insert_size: 0,

            threads: 0,

            gold_std_evaluation: if GOLDSTD_EVAL { Some(MapqEvaluation::default()) } else { None },
        }
    }
}
