
use wfa2lib_rs::aligner::{AffineAligner, AlignStatus, AlignmentScope};
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

use crate::align::{common::{penalties, Align, Heuristic, Status}, sam::Cigar};


/// The gapped-alignment heuristic. WF-adaptive prunes wavefronts far behind the best diagonal;
/// alignments are optimal unless the optimum strays that far.
pub const HEURISTIC: HeuristicStrategy = HeuristicStrategy::WfAdaptive {
    min_wavefront_length: 10,
    max_distance_threshold: 50,
    steps_between_cutoffs: 1,
};

/// Gapped alignment with `wfa2lib-rs`, the Rust port of WFA2-lib.
///
/// (The name is kept from when this wrapped the C library through `libwfa2`; the behaviour is the
/// same: gap-affine penalties, CIGAR scope, WF-adaptive, `set_max_alignment_score` as the step
/// limit, scores in WFA2-lib's convention -- the negated penalty with a zero match score.)
pub struct LIBWFA2Alignment {
    pub aligner: AffineAligner,
    pub cigar: Cigar,
    /// The penalties the aligner was built with (the aligner has no `Clone`; see `rebuilt`).
    penalties: AffinePenalties,
    /// The step limit last set through `set_max_alignment_score`.
    max_steps: i32,
    /// The ends-free span requested by `set_ends_free`, applied (clamped) at align time.
    ///
    /// WFA2 rejects an ends-free span whose free-length exceeds the sequence it applies to
    /// (`P0<=|P|`, `Tf<=|T|`, ...). The dovetail offsets the callers pass can exceed a slice that
    /// turns out shorter than the offset -- which never happened on uniform 150 bp simulated reads
    /// but does on real data (short fragments, adapters). Storing the request and clamping it
    /// against the actual `q`/`r` lengths inside `align` makes every caller safe at once.
    span: Option<(i32, i32, i32, i32)>,
}

impl LIBWFA2Alignment {
    fn new_aligner(penalties: &AffinePenalties, max_steps: i32) -> AffineAligner {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(penalties.clone()));
        aligner.set_alignment_scope(AlignmentScope::ComputeAlignment);
        // WF-adaptive (10, 50, 1): what flexalign has always run (WFA2-lib's default).
        aligner.set_heuristic(HEURISTIC);
        aligner.set_max_alignment_steps(max_steps);
        aligner
    }

    /// Align `q` (pattern) against `r` (text) with the stored span, as `wavefront_align` did.
    fn run(&mut self, q: &[u8], r: &[u8]) -> (i32, AlignStatus) {
        let wfa_score = match self.span {
            Some((qs, qe, rs, re)) => {
                let (ql, rl) = (q.len() as i32, r.len() as i32);
                self.aligner.set_alignment_free_ends(qs.min(ql), qe.min(ql), rs.min(rl), re.min(rl));
                self.aligner.align_endsfree(q, r)
            }
            None => self.aligner.align_end2end(q, r),
        };
        let status = self.aligner.status();
        // WFA2-lib's classic score: the negated penalty without a match reward; with one, the
        // score of the CIGAR in the same convention.
        let score = if self.penalties.match_ == 0 {
            -wfa_score
        } else {
            self.aligner.cigar().score_gap_affine(&self.penalties)
        };
        (score, status)
    }

    /// Map an alignment outcome to flexalign's (score, status), as the C-backed wrapper did.
    fn outcome(score: i32, status: AlignStatus) -> (i32, Status) {
        match status {
            AlignStatus::Completed => (score, Status::OK),
            // A heuristic made the end unreachable: no usable score (WFA2-lib reported i32::MIN).
            AlignStatus::Partial => (std::i32::MIN, Status::Partial),
            AlignStatus::MaxStepsReached => (std::i32::MIN, Status::Dropped),
            AlignStatus::Oom => panic!("Out of memory error."),
            other => panic!("Unexpected alignment status {:?}", other),
        }
    }
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


// The aligner owns plain heap buffers (raw pointers into its own arena) and is never shared:
// each worker thread gets its own clone.
unsafe impl Send for LIBWFA2Alignment{}

impl Clone for LIBWFA2Alignment {
    fn clone(&self) -> Self {
        // A fresh aligner with every setting (penalties, scope, heuristic, step limit, span). The
        // pipeline struct is cloned once per worker thread, so every thread must score exactly
        // like the original.
        Self {
            aligner: Self::new_aligner(&self.penalties, self.max_steps),
            cigar: self.cigar.clone(),
            penalties: self.penalties.clone(),
            max_steps: self.max_steps,
            span: self.span,
        }
    }
}

impl Align for LIBWFA2Alignment {
    fn align(&mut self, q: &[u8], r: &[u8]) -> (i32, &Cigar, Status) {
        self.cigar.0.clear();
        let (score, status) = self.run(q, r);
        if matches!(status, AlignStatus::Completed | AlignStatus::Partial) {
            self.cigar.0.extend_from_slice(self.aligner.cigar().operations_slice());
        }
        let (score, status) = Self::outcome(score, status);
        (score, &self.cigar, status)
    }

    fn align_into(&mut self, q: &[u8], r: &[u8], cigar: &mut Cigar) -> (i32, Status) {
        let (score, status) = self.run(q, r);
        if matches!(status, AlignStatus::Completed | AlignStatus::Partial) {
            cigar.0.extend_from_slice(self.aligner.cigar().operations_slice());
        }
        Self::outcome(score, status)
    }

    fn set_ends_free(&mut self, qstart: i32, qend: i32, rstart: i32, rend: i32) {
        // Stored, not applied: the clamp needs the sequence lengths, which are known only at align
        // time (see `span` / `run`).
        self.span = Some((qstart, qend, rstart, rend));
    }
}

impl LIBWFA2Alignment {
    pub fn set_penalties(&mut self, match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) {
        self.penalties = AffinePenalties { match_, mismatch, gap_opening, gap_extension };
        self.aligner = Self::new_aligner(&self.penalties, self.max_steps);
    }

    pub fn with_penalties(match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) -> Self {
        let penalties = AffinePenalties { match_, mismatch, gap_opening, gap_extension };
        Self {
            aligner: Self::new_aligner(&penalties, i32::MAX),
            cigar: Cigar(Vec::new()),
            penalties,
            max_steps: i32::MAX,
            span: None,
        }
    }

    pub fn set_below_ani_abort(&mut self, min_ani: f64, overlap_length: usize) {
        let score = ani_abort_score(min_ani, self.penalties.mismatch, overlap_length as i32);
        self.set_max_alignment_score(score);
    }
}

impl Heuristic for LIBWFA2Alignment {
    fn set_max_alignment_score(&mut self, score: i32) {
        // WFA2-lib's max alignment score is its step limit (`max_alignment_steps`).
        self.max_steps = score;
        self.aligner.set_max_alignment_steps(score);
    }
}

impl Default for LIBWFA2Alignment {
    fn default() -> Self {
        Self::with_penalties(
            penalties::MATCH,
            penalties::MISMATCH,
            penalties::GAP_OPENING,
            penalties::GAP_EXTENSION,
        )
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
