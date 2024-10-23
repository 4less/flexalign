use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};
use flexmap::values::{VData, VRange};
use kmerrs::consecutive::kmer::Kmer;

use super::{data_structures::{Alignment, Alignments, Anchor, Seed}, process::{anchor_extractor::{SeedGroupPair, SeedGroupPaired}, range_extractor::Range}, sam::{Cigar, CigarRef}, stats::Stats};

#[derive(Debug)]
pub enum Status {
    OK, Partial, Dropped
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
    let mut soft_counter = 0;

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
                soft_counter += (*c == b'S') as u32;
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
    fn generate(&mut self, kmers: &[(usize, Kmer<C>)], stats: &mut Stats) -> &[Range<F>];
    fn retrieve(&self) -> &[Range<F>];
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


pub trait AnchorAligner {
    fn align(&mut self, anchor: &Anchor) -> Alignments;
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
}
impl PairedAnchorMAPQ for StdPairedAnchorMAPQ {
    fn anchor_mapq(anchors: &mut [AnchorPair]) -> u8 {
        assert!(!anchors.is_empty());
        if anchors.len() <= 1 { return 0 };

        // Requires anchors being sorted from best to worst anchor
        let best = &anchors[0];
        let second = &anchors[1];

        (Self::score_paired(&best) - Self::score_paired(&second)) as u8
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

 pub trait SAMOutput {
    fn write();
 }

 #[derive(Clone)]
 pub struct NoSAMOutput;

 impl SAMOutput for NoSAMOutput {
    fn write() {
        todo!()
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

