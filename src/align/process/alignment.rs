
use libwfa2::affine_wavefront::{AffineWavefronts, AlignmentSpan, AlignmentStatus};

use crate::align::{common::{penalties, Align, Heuristic, Status}, sam::Cigar};


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


/// The approximate ANI an alignment cost implies -- the inverse of [`ani_abort_score`].
///
/// `score` is the aligner's score as stored on an anchor, i.e. negative or zero; `length` is the
/// query length the alignment covers. The approximation prices gap cost as though it were mismatch
/// cost, so an indel-rich alignment reads lower than its true identity. That is the same
/// simplification `ani_abort_score` makes in the other direction, which is what keeps the abort
/// threshold and the output filter consistent: an alignment that survives the abort budget derived
/// from `min_ani` also passes `approx_ani(..) >= min_ani`.
///
/// Returns 0.0 for a zero-length query and clamps at 0.0, so a dropped alignment (whose score is
/// hugely negative) reports as 0.0 identity rather than a negative number.
pub fn approx_ani(score: i32, mismatch: i32, length: usize) -> f64 {
    if length == 0 || mismatch == 0 {
        return 0.0;
    }
    let cost = (-score) as f64;
    (1.0 - cost / (length as f64 * mismatch as f64)).max(0.0)
}


/// The approximate ANI implied by a Hamming distance over `length` query bases -- the fallback for
/// an anchor that never went through gapped alignment (so has no cost to invert).
pub fn hamming_ani(hamming: u64, length: usize) -> f64 {
    if length == 0 {
        return 0.0;
    }
    1.0 - (hamming as f64 / length as f64)
}


unsafe impl Send for LIBWFA2Alignment{}

impl Clone for LIBWFA2Alignment {
    fn clone(&self) -> Self {
        // Penalties must be carried over explicitly. `AffineWavefronts::default()` comes up with
        // WFA2-lib's built-in scheme (mismatch 4, gap 6/2), *not* ours -- and since the whole
        // pipeline struct is cloned once per worker thread, a clone that forgets them leaves every
        // thread but one scoring on a different scale than `align_middle`, which charges
        // `penalties::MISMATCH` in Rust. That produced per-thread-varying scores for identical
        // alignments. It went unnoticed while the configured penalties happened to equal WFA2-lib's
        // defaults; it stops being invisible the moment they differ.
        let p = unsafe { *self.aligner.aligner() }.penalties;
        let mut new_aligner = AffineWavefronts::with_penalties(
            p.match_,
            p.mismatch,
            p.gap_opening1,
            p.gap_extension1,
        );

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
        let mut aligner = AffineWavefronts::with_penalties(
            penalties::MATCH,
            penalties::MISMATCH,
            penalties::GAP_OPENING,
            penalties::GAP_EXTENSION,
        );
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::common::penalties;

    /// The penalties the pipeline scores with must be the penalties the aligner actually charges.
    /// `libwfa2` 0.1.1 builds its aligner from WFA2-lib's compiled-in defaults and then writes into
    /// the derived `penalties` struct, which is not necessarily what the alignment kernel reads --
    /// so this asserts on observed cost, not on the field we set.
    #[test]
    fn aligner_charges_the_configured_mismatch_penalty() {
        let mut a = LIBWFA2Alignment::default();
        //                    v one mismatch, nothing else
        let q = b"ACGTACGTACGTACGT";
        let r = b"ACGTACGTTCGTACGT";
        let (score, _cigar, _status) = a.align(q, r);
        assert_eq!(
            -score,
            penalties::MISMATCH,
            "one mismatch should cost exactly penalties::MISMATCH"
        );
    }

    /// Regression: `clone()` used to start from `AffineWavefronts::default()` and drop the
    /// penalties, so worker threads (which each get a clone of the pipeline) scored mismatches at
    /// WFA2-lib's built-in 4 while the original scored at ours. Identical alignments then got
    /// different scores depending on which thread happened to handle the read.
    #[test]
    fn clone_preserves_penalties() {
        let original = LIBWFA2Alignment::default();
        let mut cloned = original.clone();

        let q = b"ACGTACGTACGTACGT";
        let r = b"ACGTACGTTCGTACGT";
        let (score, _cigar, _status) = cloned.align(q, r);
        assert_eq!(-score, penalties::MISMATCH, "a cloned aligner must charge the same mismatch penalty");

        let mut twice_cloned = cloned.clone();
        let (score, _cigar, _status) = twice_cloned.align(q, r);
        assert_eq!(-score, penalties::MISMATCH, "penalties must survive repeated cloning");
    }

    #[test]
    fn aligner_charges_the_configured_gap_penalty() {
        let mut a = LIBWFA2Alignment::default();
        // One single-base deletion.
        let q = b"ACGTACGTACGTACGT";
        let r = b"ACGTACGTCGTACGT";
        let (score, _cigar, _status) = a.align(q, r);
        assert_eq!(
            -score,
            penalties::GAP_OPENING + penalties::GAP_EXTENSION,
            "a 1bp gap should cost gap_opening + gap_extension"
        );
    }
}
