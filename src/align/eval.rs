use std::{cmp::max, fmt::Display};

use bioreader::parallel::fastq::Merge;

#[derive(Clone, Debug)]
pub struct BinaryEvaluator {
    pub tps: u64,
    pub fps: u64,
    pub tns: u64,
    pub fns: u64,
}

impl BinaryEvaluator {

    fn actual_positives(&self) -> u64 {
        self.tps + self.fns
    }

    fn actual_negatives(&self) -> u64 {
        self.fps + self.tns
    }

    fn predicted_positives(&self) -> u64 {
        self.tps + self.fps
    }

    fn predicted_negatives(&self) -> u64 {
        self.tns + self.fns
    }

    fn total(&self) -> u64 {
        self.tps + self.tns + self.fps + self.fns
    }

    fn sensitivity(&self) -> f64 {
        self.tps as f64 / self.actual_positives() as f64
    }

    fn recall(&self) -> f64 {
        self.sensitivity()
    }

    fn true_positive_rate(&self) -> f64 {
        self.sensitivity()
    }

    fn false_positive_rate(&self) -> f64 {
        self.fps as f64 / self.actual_negatives() as f64
    }

    fn true_negative_rate(&self) -> f64 {
        self.tns as f64 / self.actual_negatives() as f64
    }

    fn false_negative_rate(&self) -> f64 {
        self.fns as f64 / self.actual_negatives() as f64
    }
    
    fn negative_predictive_value(&self) -> f64 {
        self.tns as f64 / self.predicted_negatives() as f64
    }

    fn specificity(&self) -> f64 {
        self.true_negative_rate()
    }

    fn precision(&self) -> f64 {
        self.tps as f64 / self.predicted_positives() as f64
    }

    fn positive_predictive_value(&self) -> f64 {
        self.precision()
    }

    fn f1_score(&self) -> f64 {
        (2 * self.tps) as f64 / (2 * self.tps + self.fps + self.fns) as f64
    }

    fn accuracy(&self) -> f64 {
        (self.tps + self.fns) as f64 / self.total() as f64
    }
}

impl Default for BinaryEvaluator {
    fn default() -> Self {
        Self { tps: 0, fps: 0, tns: 0, fns: 0 }
    }
}

#[derive(Clone, Debug)]
pub struct MapqEvaluation {
    pub mapq_correct: Vec<u64>,
    pub mapq_incorrect: Vec<u64>,
}

impl Display for MapqEvaluation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let binary_eval = self.binary_evaluator(0);

        let mut str = String::default();
        str.push_str("MAPQ\tTP\tFP\tFN\tTN\tSensitivity\tPrecision\tF1\tSpecificity\tAccuracy\n");
        for mapq_threshold in 0..max(self.mapq_correct.len(), self.mapq_incorrect.len()) {
            let binary_eval = self.binary_evaluator(mapq_threshold);
            str.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n",
                mapq_threshold,
                binary_eval.tps,
                binary_eval.fps,
                binary_eval.fns,
                binary_eval.tns,
                binary_eval.sensitivity(),
                binary_eval.precision(),
                binary_eval.f1_score(),
                binary_eval.specificity(),
                binary_eval.accuracy(),
                binary_eval.true_negative_rate(),
                binary_eval.negative_predictive_value(),
            ));
        }

        write!(f, "{}", str)

        // write!(f,
        //     "TP:           {}\n\
        //     FP:            {}\n\
        //     FN:            {}\n\
        //     TN:            {}\n\
        //     Sensitivity:   {}\n\
        //     Precision:     {}\n\
        //     F1-Score:      {}",
        //     binary_eval.tps, 
        //     binary_eval.fps, 
        //     binary_eval.fns, 
        //     binary_eval.tns,
        //     binary_eval.sensitivity(),
        //     binary_eval.precision(),
        //     binary_eval.f1_score())
    }
}

impl MapqEvaluation {
    pub fn binary_evaluator(&self, mapq_threshold: usize) -> BinaryEvaluator {
        BinaryEvaluator {
            tps: self.mapq_correct.iter().skip(mapq_threshold).sum(),
            fps: self.mapq_incorrect.iter().skip(mapq_threshold).sum(),
            tns: self.mapq_incorrect.iter().take(mapq_threshold).sum(),
            fns: self.mapq_correct.iter().take(mapq_threshold).sum(),
        }
    }

    pub fn add(&mut self, correct: bool, mapq: u64) {
        match correct {
            true => {
                if mapq >= self.mapq_correct.len() as u64 {
                    self.mapq_correct.resize(mapq as usize + 1, 0);
                }
                self.mapq_correct[mapq as usize] += 1;
            },
            false => {
                if mapq >= self.mapq_incorrect.len() as u64 {
                    self.mapq_incorrect.resize(mapq as usize + 1, 0);
                }
                self.mapq_incorrect[mapq as usize] += 1;
            },
        }
    }
}

impl Default for MapqEvaluation {
    fn default() -> Self {
        Self { mapq_correct: Vec::new(), mapq_incorrect: Vec::new() }
    }
}

impl Merge for MapqEvaluation {
    fn merge_from(self: &mut Self, other: &mut Self) {
        if self.mapq_correct.len() < other.mapq_correct.len() {
            self.mapq_correct.resize(other.mapq_correct.len(), 0);
        }
        if self.mapq_incorrect.len() < other.mapq_incorrect.len() {
            self.mapq_incorrect.resize(other.mapq_incorrect.len(), 0);
        }
        
        for i in 0..other.mapq_correct.len() {
            self.mapq_correct[i] += other.mapq_correct[i];
        }
        for i in 0..other.mapq_incorrect.len() {
            self.mapq_incorrect[i] += other.mapq_incorrect[i];
        }
    }
}

