use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};
use kmerrs::consecutive::kmer::Kmer;

use super::{data_structures::{anchor::{Anchor, SeedGroupPair, SeedGroupPaired}, common::{Alignments, Seed}}, process::{anchor_extender::FixAnchorError, range_extractor::Range}, sam::Cigar, stats::Stats};


#[derive(Debug)]
pub enum Status {
    OK, Partial, Dropped
}


/// Alignment penalties, in WFA2's *cost* space -- costs are minimised and a match always costs 0
/// (the wavefront recurrence requires it), so bwa's scoring scheme cannot be passed through as-is.
///
/// These are bwa-mem's defaults -- match `+1`, mismatch `-4`, gap open `-6`, gap extend `-1` -- put
/// through the standard score-to-cost transform. For an alignment with `M` matches, `X` mismatches,
/// `G` gaps spanning `I+D` bases, bwa maximises
///
/// ```text
/// S = a*M - b*X - G*o - e*(I+D)
/// ```
///
/// and since `2M + 2X + I + D = |query| + |reference|`, substituting `M` out gives
///
/// ```text
/// S = a*(|query|+|reference|)/2 - (a+b)*X - o*G - (a/2 + e)*(I+D)
/// ```
///
/// Maximising that is the same as minimising `2(a+b)*X + 2o*G + (a+2e)*(I+D)`; the leading term is a
/// constant for a fixed pair of sequences and drops out, and the factor 2 keeps the values integral
/// (`a + 2e` is odd). So the ranking of alignments here is identical to bwa-mem's, on a scale
/// stretched by 2.
pub mod penalties {
    /// Always 0 -- WFA is a cost model, not a score model.
    pub const MATCH: i32 = 0;
    /// `2 * (a + b)` = `2 * (1 + 4)`.
    pub const MISMATCH: i32 = 10;
    /// `2 * o` = `2 * 6`.
    pub const GAP_OPENING: i32 = 12;
    /// `a + 2 * e` = `1 + 2 * 1`.
    pub const GAP_EXTENSION: i32 = 3;
}


pub trait Align {
    fn align(&mut self, q: &[u8], r: &[u8]) -> (i32, &Cigar, Status);
    fn align_into(&mut self, q: &[u8], r: &[u8], cigar: &mut Cigar) -> (i32, Status);
    fn set_ends_free(&mut self, qstart: i32, qend: i32, rstart: i32, rend: i32);
}

pub fn print_alignment(query: &[u8], reference: &[u8], cigar: &[u8]) {
    let mut qi = 0;
    let mut ri = 0;

    let mut q_str = String::default();
    let mut r_str = String::default();
    let mut m_str = String::default();

    cigar.iter().for_each(|i| {
        match i {
            c if *c == b'M' => {
                assert!(qi < query.len());
                assert!(ri < reference.len());
                q_str.push(query[qi] as char);
                m_str.push('|');
                r_str.push(reference[ri] as char);
                qi += 1;
                ri += 1;
            }
            c if *c == b'X' => {
                assert!(qi < query.len());
                assert!(ri < reference.len());
                q_str.push(query[qi] as char);
                m_str.push('.');
                r_str.push(reference[ri] as char);
                qi += 1;
                ri += 1;
            }
            c if *c == b'D' || *c == b'S' => {
                assert!(qi < query.len());
                q_str.push(query[qi] as char);
                m_str.push(' ');
                r_str.push('-');
                qi += 1;
            }
            c if *c == b'I' => {
                assert!(ri < reference.len());
                q_str.push('-');
                m_str.push(' ');
                r_str.push(reference[ri] as char);
                ri += 1;
            }
            c => panic!("Unknown cigar char {}", *c as char),
        }
    });

    eprintln!("{}\n{}\n{}", q_str, m_str, r_str);
    
}


pub fn is_alignment_valid(query: &[u8], reference: &[u8], cigar: &[u8]) -> bool {
    let mut qi = 0;
    let mut ri = 0;
    let mut _soft_counter = 0;

    for i in cigar {
        match i {
            c if *c == b'M' => {
                assert!(qi < query.len());
                assert!(ri < reference.len());
                if query[qi] != reference[ri] {
                    return false;
                }
                qi += 1;
                ri += 1;
            }
            c if *c == b'X' => {
                assert!(qi < query.len());
                assert!(ri < reference.len());
                qi += 1;
                ri += 1;
            }
            c if *c == b'D' || *c == b'S' => {
                _soft_counter += (*c == b'S') as u32;
                assert!(qi < query.len());
                qi += 1;
            }
            c if *c == b'I' => {
                if ri >= reference.len() {
                    eprintln!("qi: {}/{}, ri: {}/{}", qi, query.len(), ri, reference.len());
                    eprintln!("Cigar: {}", String::from_utf8_lossy(cigar));
                    eprintln!("rest q: {}", String::from_utf8_lossy(&query[qi..]));
                }
                assert!(ri < reference.len());
                ri += 1;
            }
            c => panic!("Unknown cigar char {}", *c as char),
        }
    }

    return true
}

pub trait Heuristic {
    fn set_max_alignment_score(&mut self, score: i32);
}


pub trait KmerExtractor<const K: usize> {
    fn generate(&mut self, rec: &RefFastqRecord, stats: &mut Stats) -> &[(usize, Kmer<K>)];
    fn retrieve(&self) -> &[(usize, Kmer<K>)];
}

pub trait RangeExtractor<const C: usize, const F: usize> {
    fn generate(&mut self, kmers: &[(usize, Kmer<C>)], rec: &RefFastqRecord, stats: &mut Stats) -> &[Range<'_, F>];
    fn retrieve(&self) -> &[Range<'_, F>];
}

pub trait SeedExtractor<const F: usize> {
    fn generate(&mut self, ranges: &[Range<F>], stats: &mut Stats) -> &[Seed];
    fn retrieve(&self) -> &[Seed];
}

pub trait AnchorExtractor {
    fn generate(&mut self, seeds: &[Seed], read_length: usize, stats: &mut Stats) -> &mut [Anchor];
    fn retrieve(&self) -> &[Anchor];
    fn retrieve_mut(&mut self) -> &mut [Anchor];
}


// pub type AnchorPair = (Option<Anchor>, Option<Anchor>);
#[derive(Clone, Debug)]
pub struct AnchorPair(pub Option<Anchor>, pub Option<Anchor>);

impl AnchorPair {
    pub fn resolve_orientation(&mut self, read_length_fwd: usize, read_length_rev: usize) {
        if self.0.as_ref().is_some_and(|a| !a.orientation_set) && self.1.as_ref().is_some_and(|a| a.orientation_set) {
            self.0.as_mut().unwrap().set_forward(!self.1.as_mut().unwrap().forward, read_length_fwd);
        }

        if self.1.as_ref().is_some_and(|a| !a.orientation_set) && self.0.as_ref().is_some_and(|a| a.orientation_set) {
            self.1.as_mut().unwrap().set_forward(!self.0.as_mut().unwrap().forward, read_length_rev);
        }
    }

    pub fn reference(&self) -> u64 {
        if self.0.is_some() { return self.0.as_ref().unwrap().reference } else { return self.1.as_ref().unwrap().reference }
    }

    pub fn validate(&self) -> Result<(), FixAnchorError> {
        let ok1 = self.0.as_ref().map_or(Ok(()), |a| a.validate());
        let ok2 = self.1.as_ref().map_or(Ok(()), |a| a.validate());
        ok1.and(ok2)
    }
}


pub trait PairedAnchorExtractor {
    fn generate(&mut self, seeds_fwd: &[Seed], seeds_rev: &[Seed], read_length_fwd: usize, read_length_rev: usize, stats: &mut Stats) -> &mut [AnchorPair];
    fn retrieve(&self) -> &[AnchorPair];
    fn retrieve_mut(&mut self) -> &mut [AnchorPair];
}
pub trait PairedAnchorSorter {
    fn sort(&self, anchors: &mut [AnchorPair], rec_fwd: &RefFastqRecord, rec_fwd_revc: &OwnedFastqRecord,
        rec_rev: &RefFastqRecord, rec_rev_revc: &OwnedFastqRecord, stats: &mut Stats);
}

pub trait PairedAnchorExtender {
    fn extend(&self, anchors: &mut [AnchorPair], rec_fwd: &RefFastqRecord, rec_fwd_revc: &OwnedFastqRecord,
        rec_rev: &RefFastqRecord, rec_rev_revc: &OwnedFastqRecord, stats: &mut Stats) -> usize;
}


pub trait AnchorAligner {
    fn align(&mut self, anchor: &Anchor) -> Alignments<'_>;
}

pub trait PairedAnchorMAPQ {
    fn anchor_mapq(anchors: &mut [AnchorPair]) -> u8;
}

pub trait AnchorScore {
    fn score(a: &Anchor) -> i32;
}

pub struct StdAnchorScore;
impl AnchorScore for StdAnchorScore {
    fn score(a: &Anchor) -> i32 {
        a.core_matches() as i32 - a.mismatches as i32
    }
}

pub struct StdPairedAnchorMAPQ;
impl StdPairedAnchorMAPQ {
    fn score(a: &Anchor) -> i32 {
        a.core_matches() as i32 - a.mismatches as i32
    }

    fn score_paired(a: &AnchorPair) -> i32 {
        (match &a.0 {
            Some(a) => Self::score(&a),
            None => 0,
        }) + (match &a.1 {
            Some(a) => Self::score(&a),
            None => 0,
        })
    }

    fn score_paired_ext(a: &AnchorPair) -> i32 {
        (match &a.0 {
            Some(a) => a.score,
            None => 0,
        }) + (match &a.1 {
            Some(a) => a.score,
            None => 0,
        })
    }
}
impl PairedAnchorMAPQ for StdPairedAnchorMAPQ {
    /// Phred-style pseudo-MAPQ in the SAM range `[0, 60]`.
    ///
    /// Confidence is the *relative* margin between the best and second-best anchor pair --
    /// `1 - score2/score1` -- rather than a raw score gap. The relative form is scale-invariant
    /// (independent of read length / penalties), so it spreads reads smoothly across the range
    /// instead of piling into the bimodal clusters a raw gap produces, and the explicit clamp means
    /// it can never wrap a `u8` (the old `(gap) as u8` folded any gap > 255 back down to a small
    /// value, scrambling the most confident reads).
    ///
    /// Requires `anchors` sorted best-first. A read with no competing pair has no margin evidence at
    /// all -- neither a decisive win nor a tie -- so it is given a fixed *middle* confidence rather
    /// than the maximum: on a marker DB a lone hit is empirically less reliable than a decisive
    /// margin over a real competitor, so it must not sit at the top of the MAPQ sweep.
    fn anchor_mapq(anchors: &mut [AnchorPair]) -> u8 {
        assert!(!anchors.is_empty());

        const MAPQ_MAX: f64 = 60.0;
        const MAPQ_UNIQUE: u8 = 30;

        let s1 = Self::score_paired(&anchors[0]).max(0);
        if s1 == 0 {
            return 0;
        }
        if anchors.len() <= 1 {
            return MAPQ_UNIQUE;
        }
        let s2 = Self::score_paired(&anchors[1]).max(0);
        let rel_margin = 1.0 - (s2 as f64 / s1 as f64); // in [0, 1]
        (MAPQ_MAX * rel_margin).round().clamp(0.0, MAPQ_MAX) as u8
    }
}

#[derive(Clone)]
pub struct Or<A, B> {
    pub a: Option<A>,
    pub b: Option<B>,
}

impl<A,B> Or<A,B> {
    pub fn new_a(a: A) -> Self {
        Self {
            a: Some(a),
            b: None,
        }
    }

    pub fn new_b(b: B) -> Self {
        Self {
            a: None,
            b: Some(b),
        }
    }

    pub fn a(&self) -> &A {
        &self.a.as_ref().unwrap()
    }

    pub fn b(&self) -> &B {
        &self.b.as_ref().unwrap()
    }

    pub fn has_a(&self) -> bool {
        self.a.is_some()
    }

    pub fn has_b(&self) -> bool {
        self.b.is_some()
    }
}

 pub trait PAFOutput {
    fn write(
        &mut self,
        query_name: &str,
        query_length: usize,
        query_start: i32,
        query_end: i32,
        fwd: bool,
        reference_name: &str,
        reference_length: usize,
        reference_start: i32,
        reference_end: i32,
        residue_matches: u32,
        alignment_block_length: usize,
        mapping_quality: u8,
    );
 }

 /// SAM record sink, shaped like [`PAFOutput`] so the two can sit side by side in an [`Or`].
 ///
 /// Only mapped records are written -- see `--min-ani`. `reference_start` is 0-based (the
 /// implementation converts to SAM's 1-based POS), and `cigar` carries the internal
 /// WFA-convention CIGAR, which the implementation translates and run-length encodes.
 pub trait SAMOutput {
    /// Emits `@HD` and one `@SQ` per reference. Must be called once, before any record.
    fn write_header(&mut self, references: &[(&str, usize)]);

    fn write(
        &mut self,
        query_name: &str,
        flag: u16,
        reference_name: &str,
        reference_start: usize,
        mapping_quality: u8,
        cigar: &Cigar,
        mate_reference_name: Option<&str>,
        mate_reference_start: usize,
        template_length: i64,
        seq: &[u8],
        qual: &[u8],
    );
 }

 #[derive(Clone)]
 pub struct NoSAMOutput;

 impl SAMOutput for NoSAMOutput {
    fn write_header(&mut self, _references: &[(&str, usize)]) {}

    fn write(
        &mut self,
        _query_name: &str,
        _flag: u16,
        _reference_name: &str,
        _reference_start: usize,
        _mapping_quality: u8,
        _cigar: &Cigar,
        _mate_reference_name: Option<&str>,
        _mate_reference_start: usize,
        _template_length: i64,
        _seq: &[u8],
        _qual: &[u8],
    ) {
    }
 }
 

// Weird naming. Think of something better
pub type SeedGroupPairedList = Vec<SeedGroupPaired>;
pub type SeedGroupPairList = Vec<SeedGroupPair>;

pub trait Print {
    fn print(&self);
}

impl Print for Vec<Anchor> {
    fn print(&self) {
        eprintln!("Anchor print -----");
        for a in self {
            eprintln!("\t{}", a);
        }
        eprintln!("----- Anchor print");
    }
}

impl Print for Vec<AnchorPair> {
    fn print(&self) {
        eprintln!("Anchor pair print -----");
        for AnchorPair(a1, a2) in self {
            eprintln!("\t---");
            eprintln!("\t\t{:?}", a1);
            eprintln!("\t\t{:?}", a2);
        }
        eprintln!("----- Anchor print");
    }
}

impl Print for &mut [AnchorPair] {
    fn print(&self) {
        eprintln!("Anchor pair print -----");
        for AnchorPair(a1, a2) in self.iter() {
            eprintln!("\t---");
            eprintln!("\t\t{:?}", a1);
            eprintln!("\t\t{:?}", a2);
        }
        eprintln!("----- Anchor print");
    }
}

impl Print for SeedGroupPairedList {
    fn print(&self) {
        eprintln!("Seed group print -----");
        for s in self {
            eprintln!("\t{}", s);
        }
        eprintln!("----- Seed group print");
    }
}
