use std::{cmp::max, collections::HashMap, fmt::Display, time::Duration};
use rgb::RGB8;
use bioreader::parallel::fastq::Merge;
use textplots::{Chart, ColorPlot, Plot, Shape};

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