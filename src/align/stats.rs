use std::{cmp::max, collections::HashMap, fmt::Display, time::Duration};

use bioreader::parallel::fastq::Merge;

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

    pub threads: usize,

    pub gold_std_evaluation: Option<MapqEvaluation>,
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

        self.time_get_ranges += other.time_get_ranges;
        self.time_range_header += other.time_range_header;
        self.time_get_anchors += other.time_get_anchors;
        self.time_offset += other.time_offset;
        self.time_checking_anchors += other.time_checking_anchors;

        self.ranges += other.ranges;
        self.seeds += other.seeds;
        self.anchors += other.anchors;
        self.threads += 1;

        if self.gold_std_evaluation.is_some() && other.gold_std_evaluation.is_some() {
            self.gold_std_evaluation.as_mut().unwrap().merge_from(&mut other.gold_std_evaluation.as_mut().unwrap());
        }
    }
}

impl Display for Stats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, 
            "Time for sorting ranges.....................{:?}\n\
            Time for sorting seeds......................{:?}\n\
            Time for sorting anchors....................{:?}\n\
            Time for getting ranges.....................{:?}\n\
            Time for getting range headers..............{:?}\n\
            Time for getting anchors....................{:?}\n\
            Time for getting reverse complement.........{:?}\n\
            Time for extending anchors..................{:?}\n\
            Time for calculating offsets................{:?}\n\
            Time for checking anchors...................{:?}\n\n\
            Total Reads.................................{:?}\n\
            Total Ranges per read.......................{:.2}x\n\
            Total Seeds per read........................{:.2}x\n\
            Total Anchors per read......................{:.2}x\
            {}",
            self.time_range_sorting / self.threads as u32,
            self.time_seed_sorting / self.threads as u32,
            self.time_anchor_sorting / self.threads as u32,
            self.time_get_ranges / self.threads as u32,
            self.time_range_header / self.threads as u32,
            self.time_get_anchors / self.threads as u32,
            self.time_reverse_complement / self.threads as u32,
            self.time_extend_anchors / self.threads as u32,
            self.time_offset / self.threads as u32,
            self.time_checking_anchors / self.threads as u32,
            self.reads_processed,
            self.ranges as f64 / self.reads_processed as f64,
            self.seeds as f64 / self.reads_processed as f64,
            self.anchors as f64 / self.reads_processed as f64,
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

            time_reverse_complement: Duration::default(),
            time_extend_anchors: Duration::default(),
            time_get_ranges: Duration::default(),
            time_range_sorting: Duration::default(),
            time_seed_sorting: Duration::default(),
            time_anchor_sorting: Duration::default(),
            time_range_header: Duration::default(),
            time_offset: Duration::default(),
            time_checking_anchors: Duration::default(),
            time_get_anchors: Duration::default(),
            
            threads: 0,

            gold_std_evaluation: if GOLDSTD_EVAL { Some(MapqEvaluation::default()) } else { None },
        }
    }
}