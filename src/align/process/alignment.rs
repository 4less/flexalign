
use libwfa2::{affine_wavefront::{AffineWavefronts, AlignmentSpan, AlignmentStatus, HeuristicStrategy}, bindings::wfa::wavefront_aligner_set_heuristic_xdrop};

use crate::align::{common::{Align, Heuristic, Status}, sam::{Cigar, CigarRef}};


// pub struct FastAlignment {
    
// }



pub struct LIBWFA2Alignment {
    pub aligner: AffineWavefronts,
    pub cigar: Cigar,
}


pub fn ani_abort_score(min_ani: f64, mismatch: i32, overlap_length: i32) -> i32 {
    let score = (1.0 - min_ani) * overlap_length as f64 * mismatch as f64;
    score.ceil() as i32
}


unsafe impl Send for LIBWFA2Alignment{}

impl Clone for LIBWFA2Alignment {
    fn clone(&self) -> Self {
        let mut new_aligner = AffineWavefronts::default();

        new_aligner.set_alignment_scope(self.aligner.get_alignment_scope());
        new_aligner.set_alignment_span(self.aligner.get_alignment_span());
        new_aligner.set_memory_mode(self.aligner.get_memory_mode());
        
        self.aligner.get_heuristics().iter().for_each(|h| {
            new_aligner.set_heuristic(h);
        });

        new_aligner.set_max_alignment_score(self.aligner.get_max_alignment_steps());

        Self { 
            aligner: new_aligner,
            cigar: self.cigar.clone(),
        }
    }
}

impl Align for LIBWFA2Alignment {
    fn align(&mut self, q: &[u8], r: &[u8]) -> (i32, &Cigar, Status) {
        self.cigar.0.clear();

        // Perform alignment
        match self.aligner.align(q, r) {
            AlignmentStatus::Completed => {
                self.cigar.0.extend_from_slice(self.aligner.cigar());
                (self.aligner.score(), &self.cigar, Status::OK)
            },
            AlignmentStatus::Partial => {
                self.cigar.0.extend_from_slice(self.aligner.cigar());
                if self.aligner.score() != std::i32::MIN {
                    eprintln!("Yooo score {}", self.aligner.score());
                    panic!("Yoo");
                }
                if !self.cigar.0.is_empty() {
                    eprintln!("{} {} {}", self.aligner.score(), self.cigar.0.len(), String::from_utf8_lossy(&self.cigar.0))
                }

                (self.aligner.score(), &self.cigar, Status::Partial)
            },
            AlignmentStatus::MaxStepsReached => {
                (std::i32::MIN, &self.cigar, Status::Dropped)
            },
            AlignmentStatus::OOM => panic!("Out of memory error."),
            AlignmentStatus::Unattainable => panic!("Alignment status unattainable"),
            AlignmentStatus::Undefined => panic!("Undefined alignment status"),
        }
    }
    
    fn align_into(&mut self, q: &[u8], r: &[u8], cigar: &mut Cigar) -> (i32, Status) {
        // Perform alignment
        match self.aligner.align(q, r) {
            AlignmentStatus::Completed => {
                cigar.0.extend_from_slice(self.aligner.cigar());
                (self.aligner.score(), Status::OK)
            },
            AlignmentStatus::Partial => {
                cigar.0.extend_from_slice(self.aligner.cigar());
                (self.aligner.score(), Status::Partial)
            },
            AlignmentStatus::MaxStepsReached => {
                (std::i32::MIN, Status::Dropped)
            },
            AlignmentStatus::OOM => panic!("Out of memory error."),
            AlignmentStatus::Unattainable => panic!("Alignment status unattainable"),
            AlignmentStatus::Undefined => panic!("Undefined alignment status"),
        }
    }
    
    fn set_ends_free(&mut self, qstart: i32, qend: i32, rstart: i32, rend: i32) {
        self.aligner.set_alignment_span(AlignmentSpan::EndsFree { pattern_begin_free: qstart, pattern_end_free: qend, text_begin_free: rstart, text_end_free: rend });
    }


}

impl LIBWFA2Alignment {
    pub fn set_penalties(&mut self, match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) {
        self.aligner.set_penalties(match_, mismatch, gap_opening, gap_extension);
    }

    pub fn with_penalties(match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) -> Self {
        Self {
            aligner: AffineWavefronts::with_penalties(match_, mismatch, gap_opening, gap_extension),
            cigar: Cigar(Vec::new()),
        }
    }

    pub fn set_below_ani_abort(&mut self, min_ani: f64, overlap_length: usize) {
        // (std::ceil((1 - min_ani) * static_cast<double>(overlap_length)) * mismatch_penalty) + 1;
        let mismatch = unsafe{ *self.aligner.aligner() }.penalties.mismatch;
        let score = ani_abort_score(min_ani, mismatch, overlap_length as i32);
        self.aligner.set_max_alignment_score(score);
    }
}

impl Heuristic for LIBWFA2Alignment {
    fn set_max_alignment_score(&mut self, score: i32) {
        self.aligner.set_max_alignment_score(score);
    }
}

impl Default for LIBWFA2Alignment {
    fn default() -> Self {
        let mut aligner = AffineWavefronts::with_penalties(0, 4, 6, 2);
        // aligner.set_heuristic(&HeuristicStrategy::XDrop { xdrop: std::i32::MIN, score_steps: 2 });
        // aligner.set_heuristic(&HeuristicStrategy::BandedStatic { band_min_k: -1, band_max_k: 1 });
        aligner.set_alignment_scope(libwfa2::affine_wavefront::AlignmentScope::Alignment);
        aligner.set_alignment_span(libwfa2::affine_wavefront::AlignmentSpan::End2End);

        // unsafe { wavefront_aligner_set_heuristic_xdrop(aligner.aligner_mut(), std::i32::MIN, 2) };

        Self { 
            aligner: aligner,
            cigar: Cigar(Vec::new()),
        }
    }
}