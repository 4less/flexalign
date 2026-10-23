//! Gapped alignment with `wfa2lib-rs`, the pure-Rust port of WFA2-lib, as a drop-in for
//! [`LIBWFA2Alignment`](super::alignment::LIBWFA2Alignment) (the C library through `libwfa2`).
//!
//! It is configured to do exactly what the libwfa2 aligner does today, so a run with either differs
//! only in the implementation: the same penalties (`penalties::*`), alignment scope (score and
//! CIGAR), end-to-end span unless `set_ends_free` asks otherwise (clamped per call, as there), the
//! abort score as WFA2's max alignment steps, and the same heuristic,
//! [`HEURISTIC`](super::alignment::HEURISTIC): WF-adaptive (min wavefront length 10, max distance
//! 50, every step).
//!
//! Scores follow the libwfa2 convention the pipeline is written against: the negated cost (0 for an
//! exact match), `i32::MIN` for an alignment dropped by the abort score.

use wfa2lib_rs::aligner::{AffineAligner, AlignStatus, AlignmentScope, AlignmentSpan};
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

use super::alignment::ani_abort_score;
use crate::align::{common::{penalties, Align, Heuristic, Status}, sam::Cigar};

/// The libwfa2 aligner's [`HEURISTIC`](super::alignment::HEURISTIC), in wfa2lib-rs terms.
const DEFAULT_HEURISTIC: HeuristicStrategy = HeuristicStrategy::WfAdaptive {
    min_wavefront_length: 10,
    max_distance_threshold: 50,
    steps_between_cutoffs: 1,
};

pub struct WFA2RSAlignment {
    aligner: AffineAligner,
    penalties: AffinePenalties,
    max_score: i32,
    cigar: Cigar,
    /// The ends-free span requested by `set_ends_free`, applied (clamped) at align time; see
    /// `LIBWFA2Alignment::span`. Like there, once set it applies to every later alignment.
    span: Option<(i32, i32, i32, i32)>,
}

impl WFA2RSAlignment {
    pub fn with_penalties(match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) -> Self {
        let penalties = AffinePenalties { match_, mismatch, gap_opening, gap_extension };
        Self::configured(penalties, i32::MAX, None)
    }

    fn configured(penalties: AffinePenalties, max_score: i32, span: Option<(i32, i32, i32, i32)>) -> Self {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(penalties));
        aligner.set_alignment_scope(AlignmentScope::ComputeAlignment);
        aligner.set_heuristic(DEFAULT_HEURISTIC);
        aligner.set_max_alignment_steps(max_score);
        Self { aligner, penalties, max_score, cigar: Cigar(Vec::new()), span }
    }

    pub fn set_penalties(&mut self, match_: i32, mismatch: i32, gap_opening: i32, gap_extension: i32) {
        *self = Self::configured(
            AffinePenalties { match_, mismatch, gap_opening, gap_extension },
            self.max_score,
            self.span,
        );
    }

    pub fn set_below_ani_abort(&mut self, min_ani: f64, overlap_length: usize) {
        let score = ani_abort_score(min_ani, self.penalties.mismatch, overlap_length as i32);
        self.set_max_alignment_score(score);
    }

    /// Runs the alignment; returns the libwfa2-convention score and status. The CIGAR is left in
    /// the aligner.
    fn run(&mut self, q: &[u8], r: &[u8]) -> (i32, Status) {
        let cost = match self.span {
            Some((qs, qe, rs, re)) => {
                let (ql, rl) = (q.len() as i32, r.len() as i32);
                self.aligner.set_alignment_free_ends(qs.min(ql), qe.min(ql), rs.min(rl), re.min(rl));
                self.aligner.align_endsfree(q, r)
            }
            None => {
                self.aligner.alignment_span = AlignmentSpan::End2End;
                self.aligner.align_end2end(q, r)
            }
        };
        match self.aligner.status() {
            AlignStatus::Completed => (-cost, Status::OK),
            AlignStatus::Partial => (-cost, Status::Partial),
            AlignStatus::MaxStepsReached => (i32::MIN, Status::Dropped),
            AlignStatus::Oom => panic!("Out of memory error."),
            AlignStatus::Unattainable => panic!("Alignment status unattainable"),
            other => panic!("Unexpected alignment status {other:?}"),
        }
    }
}

// SAFETY: wfa2lib-rs's aligner is not `Send` only because of raw sequence pointers (valid during
// one `align` call, which never spans a thread move) and the boxed match function of
// `align_lambda`, which is never set here. Each worker thread owns its own clone and nothing is
// shared, as with `LIBWFA2Alignment`.
unsafe impl Send for WFA2RSAlignment {}

impl Clone for WFA2RSAlignment {
    /// A fresh aligner with the same configuration (the wavefront memory is per aligner, and the
    /// pipeline clones one per worker thread).
    fn clone(&self) -> Self {
        let mut c = Self::configured(self.penalties, self.max_score, self.span);
        c.cigar = self.cigar.clone();
        c
    }
}

impl Align for WFA2RSAlignment {
    fn align(&mut self, q: &[u8], r: &[u8]) -> (i32, &Cigar, Status) {
        self.cigar.0.clear();
        let (score, status) = self.run(q, r);
        if !matches!(status, Status::Dropped) {
            self.cigar.0.extend_from_slice(self.aligner.cigar().operations_slice());
        }
        (score, &self.cigar, status)
    }

    fn align_into(&mut self, q: &[u8], r: &[u8], cigar: &mut Cigar) -> (i32, Status) {
        let (score, status) = self.run(q, r);
        if !matches!(status, Status::Dropped) {
            cigar.0.extend_from_slice(self.aligner.cigar().operations_slice());
        }
        (score, status)
    }

    fn set_ends_free(&mut self, qstart: i32, qend: i32, rstart: i32, rend: i32) {
        self.span = Some((qstart, qend, rstart, rend));
    }
}

impl Heuristic for WFA2RSAlignment {
    fn set_max_alignment_score(&mut self, score: i32) {
        self.max_score = score;
        self.aligner.set_max_alignment_steps(score);
    }
}

impl Default for WFA2RSAlignment {
    fn default() -> Self {
        Self::with_penalties(penalties::MATCH, penalties::MISMATCH, penalties::GAP_OPENING, penalties::GAP_EXTENSION)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::process::alignment::LIBWFA2Alignment;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn base(&mut self) -> u8 {
            b"ACGT"[(self.next() % 4) as usize]
        }
    }

    /// A random reference window and a read from it with substitutions and short indels.
    fn pair(rng: &mut Rng, len: usize, edits: usize) -> (Vec<u8>, Vec<u8>) {
        let r: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        let mut q = r.clone();
        for _ in 0..edits {
            let i = (rng.next() as usize) % q.len().max(1);
            match rng.next() % 4 {
                0 | 1 => q[i] = rng.base(),
                2 => q.insert(i, rng.base()),
                _ => {
                    q.remove(i);
                }
            }
        }
        (q, r)
    }

    /// The point of the type: on the same inputs it gives libwfa2's scores, statuses and CIGARs.
    #[test]
    fn matches_libwfa2_end2end_and_ends_free() {
        let mut rng = Rng(0x9E3779B97F4A7C15);
        let (mut c, mut rs) = (LIBWFA2Alignment::default(), WFA2RSAlignment::default());
        let mut checked = 0;
        for case in 0..2000 {
            let (q, r) = pair(&mut rng, 30 + case % 140, case % 9);
            let (sc, cc, stc) = c.align(&q, &r);
            let cc = cc.0.clone();
            let (sr, cr, str_) = rs.align(&q, &r);
            assert_eq!(sc, sr, "case {case}: end2end score");
            assert_eq!(std::mem::discriminant(&stc), std::mem::discriminant(&str_), "case {case}: status");
            assert_eq!(cc, cr.0, "case {case}: end2end cigar");
            checked += 1;
        }
        // Ends-free, as the flank alignment calls it.
        c.set_ends_free(0, 8, 0, 8);
        rs.set_ends_free(0, 8, 0, 8);
        for case in 0..2000 {
            let (q, mut r) = pair(&mut rng, 40 + case % 100, case % 7);
            r.extend((0..case % 6).map(|_| rng.base()));
            let (sc, cc, _) = c.align(&q, &r);
            let cc = cc.0.clone();
            let (sr, cr, _) = rs.align(&q, &r);
            assert_eq!(sc, sr, "case {case}: ends-free score");
            assert_eq!(cc, cr.0, "case {case}: ends-free cigar");
            checked += 1;
        }
        assert_eq!(checked, 4000);
    }

    #[test]
    fn the_abort_score_drops_like_libwfa2() {
        let mut rng = Rng(7);
        for case in 0..500 {
            let (q, r) = pair(&mut rng, 150, case % 25);
            let (mut c, mut rs) = (LIBWFA2Alignment::default(), WFA2RSAlignment::default());
            c.set_below_ani_abort(0.95, 150);
            rs.set_below_ani_abort(0.95, 150);
            let (sc, _, stc) = c.align(&q, &r);
            let (sr, _, str_) = rs.align(&q, &r);
            assert_eq!(std::mem::discriminant(&stc), std::mem::discriminant(&str_), "case {case}: status");
            assert_eq!(sc, sr, "case {case}: score");
        }
    }

    #[test]
    fn clone_keeps_the_configuration() {
        let mut a = WFA2RSAlignment::default();
        a.set_max_alignment_score(5);
        let mut b = a.clone();
        let q = b"ACGTACGTACGTACGTACGT";
        let r = b"ACGAACGTTCGTACCTACGT"; // 3 mismatches: cost 30 > 5
        assert!(matches!(b.align(q, r).2, Status::Dropped));
        let (score, _, _) = WFA2RSAlignment::default().align(q, b"ACGTACGTTCGTACGTACGT");
        assert_eq!(-score, penalties::MISMATCH);
    }
}
