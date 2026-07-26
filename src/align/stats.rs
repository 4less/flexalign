use std::{fmt::Display, time::Duration};
use rgb::RGB8;
use bioreader::parallel::fastq::Merge;
use textplots::{Chart, ColorPlot, Shape};

use crate::GOLDSTD_EVAL;

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

    pub time_get_kmers: Duration,
    pub time_get_minimizer: Duration,
    pub time_get_vranges: Duration,
    pub time_get_ranges: Duration,
    pub time_range_sorting: Duration,
    pub time_seed_sorting: Duration,
    pub time_anchor_sorting: Duration,
    pub time_reverse_complement: Duration,
    pub time_extend_anchors: Duration,
    pub time_get_anchors: Duration,
    pub time_range_header: Duration,
    pub time_offset: Duration,
    pub time_checking_anchors: Duration,
    pub time_alignment: Duration,

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
            ("extend_anchors",     secs(self.time_extend_anchors), "s"),
            ("offsets",            secs(self.time_offset), "s"),
            ("checking_anchors",   secs(self.time_checking_anchors), "s"),
            ("alignment",          secs(self.time_alignment), "s"),
            ("header_elements",         self.header_elements as f64, "count"),
            // Total (thread-summed) seed-extraction time over total counts -- CPU time per unit.
            ("seedext_s_per_range",     self.time_range_header.as_secs_f64() / self.ranges.max(1) as f64, "s"),
            ("seedext_s_per_header",    self.time_range_header.as_secs_f64() / self.header_elements.max(1) as f64, "s"),
            ("reads",                   self.reads_processed as f64, "count"),
            ("alignments",              self.alignments as f64, "count"),
            ("alignments_successful",   self.alignments_successful as f64, "count"),
            ("alignments_partial",      self.alignments_partial as f64, "count"),
            ("alignments_dropped",      self.alignments_dropped as f64, "count"),
            ("threads",                 self.threads as f64, "count"),
            ("minimizers_per_read",     self.minimizer as f64 / reads, "rate"),
            ("ranges_per_read",         self.ranges as f64 / reads, "rate"),
            ("retrieved_ranges_per_read", self.retrieved_ranges as f64 / reads, "rate"),
            ("seeds_per_read",          self.seeds as f64 / reads, "rate"),
            ("anchors_per_read",        self.anchors as f64 / reads, "rate"),
        ];

        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
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
        self.time_extend_anchors += other.time_extend_anchors;

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
            Time for extending anchors..................{:?}\n\
            Time for calculating offsets................{:?}\n\
            Time for checking anchors...................{:?}\n\
            Time for alignment..........................{:?}\n\n\
            Total Reads.................................{:?}\n\
            Total Alignments............................{:?}\n\
            Total Alignments successful.................{:?}\n\
            Total Alignments partial....................{:?}\n\
            Total Alignments dropped....................{:?}\n\
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
            self.time_extend_anchors / self.threads as u32,
            self.time_offset / self.threads as u32,
            self.time_checking_anchors / self.threads as u32,
            self.time_alignment / self.threads as u32,
            self.reads_processed,
            self.alignments,
            self.alignments_successful,
            self.alignments_partial,
            self.alignments_dropped,
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

            time_reverse_complement: Duration::default(),
            time_extend_anchors: Duration::default(),
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
            
            threads: 0,

            gold_std_evaluation: if GOLDSTD_EVAL { Some(MapqEvaluation::default()) } else { None },
        }
    }
}
