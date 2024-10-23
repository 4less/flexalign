use std::{cmp::{max, min}, fmt::Display, iter::{zip, Zip}, ops::Range, slice::Iter, thread::current};

use bioreader::sequence::{fastq_record::{OwnedFastqRecord, RefFastqRecord}, utils::reverse_complement_into_vec};
use colored::{Color, Colorize};
use thiserror::Error;
use triple_accel::hamming as triple_hamming;

use crate::align::common::Status;

use super::{common::{print_alignment, Align, Heuristic}, errors::{AlignmentError, AlignmentResult}, sam::Cigar};


#[derive(Debug, Clone, Error)]
struct ResolveOrientationError;
impl Display for ResolveOrientationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error while resolving the orientation for an anchor")
    }
}

#[derive(Clone, Debug)]
#[repr(C)]
pub struct Seed {
    pub rpos: u64,
    pub rval: u64,
    pub qpos: u32,
    pub mismatch: u8,
    pub length: u8,
    pub flag: u8,
}

impl Seed {
    #[inline(always)]
    pub fn from_flexmer<const K: usize, const C: usize, const F: usize>(qpos: usize, rpos: u64, reference: u64, dist: u32) -> Self {
        Self {
            qpos: if dist == 0 { qpos as u32 } else { qpos as u32 + (F as u32/2) },
            rpos: if dist == 0 { rpos as u64 } else { rpos + (F as u64/2) },
            rval: reference,
            mismatch: dist as u8,
            length: if dist == 0 {K as u8} else {C as u8},
            flag: 0,
        }
    }

    #[inline(always)]
    pub fn from_coremer<const K: usize, const C: usize, const F: usize>(qpos: usize, rpos: u64, reference: u64) -> Self {
        Self {
            qpos: qpos as u32 + (F as u32/2),
            rpos:  rpos + (F as u64/2),
            rval: reference,
            mismatch: 0,
            length: C as u8,
            flag: 0,
        }
    }

    pub fn offset(&self) -> u64 {
        self.rpos as u64 - self.qpos as u64
    }

    pub fn offsets(&self,  read_length: usize) -> (i64, i64) {
        // ((self.rpos as u64 - self.qpos as u64) | (1 << 62), self.rpos as u64 + self.qpos as u64)

        (self.rpos as i64 - self.qpos as i64, self.rpos as i64 - (read_length as i64 - self.length as i64 - self.qpos as i64))
    }

    pub fn offset_dist(&self, other: &Self, read_length: usize) -> u64 {
        let (oa1, oa2) = self.offsets(read_length);
        let (ob1, ob2) = other.offsets(read_length);
        let min1 = min(oa1.abs_diff(ob1), oa1.abs_diff(ob2));
        let min2 = min(oa2.abs_diff(ob1), oa2.abs_diff(ob2));
        min(min1, min2)
    }

    pub fn closest_offset(&self, other: &Self, read_length: usize) -> (i64, bool, u64) {
        let (self_fwd, self_rev) = self.offsets(read_length);
        let (other_fwd, other_rev) = other.offsets(read_length);

        let diff_fwd = self_fwd.abs_diff(other_fwd);
        let diff_rev = self_rev.abs_diff(other_rev);

        if diff_fwd < diff_rev {
            (self_fwd, true, diff_fwd)
        } else {
            (self_rev, false, diff_rev)
        }
    }

    pub fn reverse(&self, read_length: usize) -> Seed {
        Seed {
            qpos: read_length as u32 - self.length as u32 - self.qpos as u32,
            rpos: self.rpos,
            rval: self.rval,
            mismatch: self.mismatch,
            length: self.length,
            flag: 0,
        }
    }

    pub fn to_visual_string_x(&self, read_length: Option<usize>) -> String {
        let mut output = String::new();

        match read_length {
            Some(read_length) => {
                let spaces = String::from_utf8(vec![b' '; read_length - self.length as usize - self.qpos as usize]).unwrap();
                let xes = String::from_utf8(vec![b'X'; self.length as usize]).unwrap();
                output += &spaces;
                output += &xes;
                output += "          ";
                output += &self.to_string();
            },
            None => {
                let spaces = String::from_utf8(vec![b' '; (self.qpos) as usize]).unwrap();
                let xes = String::from_utf8(vec![b'X'; self.length as usize]).unwrap();
                output += &spaces;
                output += &xes;
            },
        };

        output
    }


    pub fn to_visual_string(&self, read: &[u8]) -> String {
        let mut output = String::new();

        let spaces = String::from_utf8(vec![b' '; (self.qpos) as usize]).unwrap();
        let xes = String::from_utf8_lossy(&read[self.qpos as usize..self.qpos as usize + self.length as usize]);
        output += &spaces;
        output += &xes;

        output
    }
}

impl Display for Seed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "reference: {}  rpos: {},  qpos: {}, mismatch: {}, length: {}, offsets: {:?}", 
            self.rval,
            self.rpos,
            self.qpos,
            self.mismatch, 
            self.length,
            self.offsets(150)) // CHANGE!!!
    }
}


//                       overlap                      containment  
//        OFFSET_FWD_OTHER     OFFSET_FWD_SELF   CONTAINED_OTHER    CONTAINED_SELF
//  self:   .........              ..........       .........            ...
//  other:     ............     .....                  ....           ...........
//  Ignore c, as other is fully contained in self) 
pub enum SeedOverlap {
    OffsetFwdOther,
    OffsetFwdSelf,
    ContainedOther,
    ContainedSelf,
    NoOverlap
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct AnchorSeed {
    pub qpos: u32,
    pub rpos: u64,
    pub length: u32,
}

impl Display for AnchorSeed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "qpos: {}, rpos: {}, length: {}",
            self.qpos, self.rpos, self.length)
    }
}

impl AnchorSeed {
    pub fn qbegin(&self) -> usize {
        self.qpos as usize
    }

    pub fn qend(&self) -> usize {
        self.qpos as usize + self.length as usize
    }

    pub fn rbegin(&self) -> usize {
        self.rpos as usize
    }
    
    pub fn rend(&self) -> usize {
        self.rpos as usize + self.length as usize
    }

    pub fn qrange(&self) -> Range<usize> {
        (self.qpos as usize)..(self.qpos + self.length) as usize
    }

    pub fn rrange(&self) -> Range<usize> {
        (self.rpos as usize)..(self.rpos + self.length as u64) as usize
    }

    pub fn extend_left(&mut self, by: usize) {
        self.qpos -= by as u32;
        self.rpos -= by as u64;
        self.length += by as u32;
    }

    pub fn extend_right(&mut self, by: usize) {
        self.length += by as u32;
    }


    pub fn set(&mut self, other: &Self) {
        self.qpos = other.qpos;
        self.rpos = other.rpos;
        self.length = other.length
    }

    pub fn merge_into(&mut self, other: &Self) -> bool {
        let self_start = self.qpos;
        let self_end = self.qpos + self.length;
        let other_start = other.qpos;
        let other_end = other.qpos + other.length;
        let mut has_overlap = false;

        //                       overlap                      containment  
        //              a                      b             c             d
        //  self:   .........              ..........   .........       ...
        //  other:     ............     .....              ....       ...........
        //  Ignore c, as other is fully contained in self)  
        if other_start >= self_start && other_start <= self_end && other_end > self_end { // a + d
            let overlap = other_end - self_end;
            self.length += overlap;
            has_overlap = true;
        }
        if other_start < self_start && other_end >= self_start && other_end <= self_end { // b + d
            let overlap = other_start - self_start;
            self.qpos = other.qpos;
            self.rpos = other.rpos;
            self.length += overlap;
            has_overlap = true;
        }
        if other_start >= self_start && other_end <= self_end { // c
            has_overlap = true;
        }

        return has_overlap
    }

    pub fn contains(&self, other: &Self) -> bool {
        self.qbegin() <= other.qbegin() && self.qend() >= other.qend()
    }

    pub fn rpos_sorted_merge_into(&mut self, other: &Self) -> SeedOverlap {
        //             overlap               containment  
        //              a                        c      
        //  self:   .........                .........    
        //  other:     ............            ....     
        //  Ignore c, as other is fully contained in self)

        if other.qpos < self.qpos {
            eprintln!("Weird!!!");
        }
        
        assert!(other.qpos >= self.qpos || other.length > self.length);

        let self_start = self.qpos;
        let self_end = self.qpos + self.length;
        let other_start = other.qpos;
        let other_end = other.qpos + other.length;

    
        if other_start >= self_start && other_end <= self_end {
            // eprintln!("{}    {} {} {} {}", "SeedOverlap::ContainedOther", self_start, self_end, other_start, other_end);
            return SeedOverlap::ContainedOther;
        }
        
        if other_start >= self_start && other_start <= self_end && other_end > self_end { 
            let overlap = other_end - self_end;
            self.length += overlap;
            // eprintln!("{}    {} {} {} {}", "SeedOverlap::OffsetFwdOther", self_start, self_end, other_start, other_end);
            return SeedOverlap::OffsetFwdOther;
        }

        if self_start >= other_start && self_end <= other_end {
            self.qpos = other.qpos;
            self.length = other.length;
            // eprintln!("{}    {} {} {} {}", "SeedOverlap::ContainedSelf", self_start, self_end, other_start, other_end);
            return SeedOverlap::ContainedSelf;
        }

        // eprintln!("{}    {} {} {} {}", "SeedOverlap::NoOverlap", self_start, self_end, other_start, other_end);
        SeedOverlap::NoOverlap
    }

    pub fn reverse(&mut self, read_length: usize) {
        // eprintln!("Reverse: {} {} {} ", read_length, self.length, self.qpos);
        self.qpos = read_length as u32 - self.length - self.qpos;
        // eprintln!("Reverse: -> {}", self.qpos);
    }
    
    pub fn offset(&self) -> i32 {
        (self.qbegin() as i64 - self.rbegin() as i64) as i32
    }

    pub fn offsets(&self,  read_length: usize) -> (i64, i64) {
        (self.rbegin() as i64 - self.qbegin() as i64, self.rbegin() as i64 - (read_length as i64 - self.length as i64 - self.qbegin() as i64))
    }
}


#[derive(Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Anchor {
    pub reference: u64, // 8
    pub seed_count: u32,
    pub mismatches: u32, // 16
    pub forward: bool,
    pub orientation_set: bool,
    pub flagged_for_indel: bool,
    pub flag: u8, //
    pub counter1: u16,
    pub counter2: u16, // 24
    pub seeds: Vec<AnchorSeed>, // 40
    pub score: i32, // 44
    pub cigar: Option<Cigar>,
    pub reference_cigar_range: Range<usize>,
}




pub trait ToString {
    fn ts(&self) -> String;
}

impl ToString for &[u8] {
    fn ts(&self) -> String {
        String::from_utf8_lossy(self).to_string()
    }
}

pub fn hamming(query: &[u8], reference: &[u8]) -> u64 {
    zip(query, reference).fold(0, |acc, (a,b)| acc + (a != b) as u64)
}

impl Anchor {
    #[inline(always)]
    pub fn from_seed(seed: &Seed) -> Self {
        Self {
            reference: seed.rval,
            seed_count: 1,
            mismatches: seed.mismatch as u32,
            forward: true,
            orientation_set: false,
            flagged_for_indel: false,
            flag: 0,
            counter1: 0,
            counter2: 0,
            score: 0,
            seeds: vec! [AnchorSeed{ qpos: seed.qpos, rpos: seed.rpos, length: seed.length as u32 }],
            cigar: None,
            reference_cigar_range: 0..0,
        }
    }

    pub fn whole_align(&mut self, aligner: &mut (impl Align + Heuristic), query: &[u8], reference: &[u8], free_ends: usize, mut max_score: i32) -> Status {
        let (mut qr, mut rr) = self.whole(query.len(), reference.len());
        
        self.cigar = Some(Cigar::new());

        // Left dove 
        let ql_dove = min(free_ends, qr.start);
        let rl_dove = min(free_ends, rr.start);

        // Right dove new ends
        let q_end = min(qr.end + free_ends, query.len());
        let r_end = min(rr.end + free_ends, reference.len());

        // Right dove sizes
        let qr_dove = q_end - qr.end;
        let rr_dove = r_end - rr.end;

        qr.start -= ql_dove;
        rr.start -= rl_dove;
        qr.end = q_end;
        rr.end = r_end;

        let q = &query[qr.clone()];
        let r = &reference[rr.clone()];


        aligner.set_ends_free(ql_dove as i32, qr_dove as i32, rl_dove as i32, rr_dove as i32);

        let (score, cigar, status) = aligner.align(q, r);

        if matches!(status, Status::OK) {
            self.score = score;
            self.cigar.as_mut().unwrap().0.extend_from_slice(&cigar.0);
        }
        
        // let q_inserts = cigar.count_trailing_chars(b'I');
        // let r_deletions = cigar.count_trailing_chars(b'D');
        // let q_softclip: usize = min(r_deletions, qr_dove);
        // let r_offset: usize = min(q_inserts, r_dove);

        // // Update cigar and set reference starting point
        // let lcigar = self.cigar.as_mut().unwrap();
        // lcigar.0.extend_from_slice(&cigar.0[0..cigar.0.len() - (r_offset + q_softclip)]);
        // lcigar.add_softclip(q_softclip);
        // self.reference_cigar_range.end = rr.1.end - r_offset;


        status
    }

    pub fn smart_align(&mut self, aligner: &mut (impl Align + Heuristic), query: &[u8], reference: &[u8], free_ends: usize, mut max_score: i32) -> Status {
        // Accurate alignment of flanks first.
        // Add threshold later and do hamming first, and if the score can possibly improve with perfect alignment, do that

        self.cigar = Some(Cigar::new());

        let mut alignment_score = 0;

        // eprintln!("Left {}", max_score);
        aligner.set_max_alignment_score(max_score + 1);

        // eprintln!("Align with max score {}", max_score + 1);
        let (score,  status, qs, rs) = self.align_left_flank(aligner, query, reference, free_ends);
        match status { 
            Status::OK => {
                assert!(score != std::i32::MIN);
            }, 
            _ => { 
                // eprintln!("Drop after left {}", score);
                return Status::Dropped
            },
        }
        // eprintln!("Max score after left: {} .... {}", max_score, score);
        max_score += score;
        alignment_score += score;

        // eprintln!("Max score before middle: {}", max_score);
        let (score, status) = match self.align_middle(query, reference, &mut max_score) {
            Ok(res) => res,
            Err(ar) => {
                match ar {
                    AlignmentError::QueryRangeError(s) => {
                        println!("Error: {}", s);
                        println!("Q: {}", String::from_utf8_lossy(query));
                        println!("Self: {}", self);
                    },
                    AlignmentError::ReferenceRangeError(s) => {
                        println!("Error: {}", s);
                        println!("Q: {}", String::from_utf8_lossy(query));
                        println!("Self: {}", self);
                    },
                    AlignmentError::InvalidAlignmentError(s) => {
                        println!("Error: {}", s);
                        println!("Q: {}", String::from_utf8_lossy(query));
                        println!("Self: {}", self);
                    },
                    AlignmentError::InvalidRangeError(s) => {
                        println!("Error: {}", s);
                        println!("Q: {}", String::from_utf8_lossy(query));
                        println!("Self: {}", self);
                    },
                };
                panic!("Non recoverable error")
            },
        };
        alignment_score += score;

        match status { 
            Status::OK => {
            }, 
            _ => { 
                // eprintln!("Drop after middle {}", score);
                return Status::Dropped
            },
        }

        aligner.set_max_alignment_score(max_score + 1);
        let (score, status) = self.align_right_flank(aligner, query, reference, free_ends);
        
        match status { 
            Status::OK => {
                assert!(self.reference_cigar_range.start < self.reference_cigar_range.end);
            }, 
            _ => { 
                // eprintln!("Drop after right {}", score);
                return Status::Dropped
            },
        }

        alignment_score += score;

        // eprintln!("Glorious Test {}", alignment_score);
        // print_alignment(&query, &reference[self.reference_cigar_range.clone()], &self.cigar().0);
        self.score = alignment_score;
        
        Status::OK
    }


    pub fn align_left_flank(&mut self, aligner: &mut impl Align, query: &[u8], reference: &[u8], free_ends: usize) -> (i32, Status, usize, usize) {
        

        // Accurate alignment of flanks first.
        // Add threshold later and do hamming first, and if the score can possibly improve with perfect alignment, do that

        let mut lr = self.left_flank();

        if lr.0.len() == 0 || lr.1.len() == 0 {
            let lcigar = self.cigar.as_mut().unwrap();
            lcigar.add_softclip(max(lr.0.len(), lr.1.len()));
            self.reference_cigar_range.start = lr.1.start - max(lr.0.len(), lr.1.len());
            // eprintln!("{:?} {:?}", lr.0, lr.1);
            // Assumes match penalty score is 0, and other scores are negative.
            return (0, Status::OK, 0, 0)
        }

        let q_dove = min(free_ends, lr.0.start);
        let r_dove = min(free_ends, lr.1.start);

        // eprintln!("---LEFT---\nqdove = {}, rdove = {} .... {} {}", q_dove, r_dove, lr.0.start, lr.1.start);
        lr.0.start -= q_dove;
        lr.1.start -= r_dove;


        let q = &query[lr.clone().0];
        let r = &reference[lr.clone().1];

        
        // eprintln!("{}{}\n{}{}",
        //     " ".repeat(max(0, r.len() as i32 - q.len() as i32) as usize), String::from_utf8_lossy(q), 
        //     " ".repeat(max(0, q.len() as i32 - r.len() as i32) as usize), String::from_utf8_lossy(r));

        aligner.set_ends_free(q_dove as i32, 0, r_dove as i32, 0);
        // aligner.set_ends_free(100,100,100,100);

        let (score, cigar, status) = aligner.align(q, r);
        
        if !matches!(status, Status::OK) {
            return (std::i32::MIN, Status::Dropped, 0, 0)
        }

        assert!(score != std::i32::MIN);

        
        // eprintln!("S-----------------------------");
        // eprintln!("Score: {} (Free Q {}, Free R {})", score, q_dove, r_dove);
        // print_alignment(q, r, &cigar.0);



        let q_inserts = cigar.count_leading_chars(b'I');
        let r_deletions = cigar.count_leading_chars(b'D');
        let q_softclip: usize = min(r_deletions, q_dove);
        let r_offset: usize = min(q_inserts, r_dove);
        
        // Update cigar and set reference starting point
        let lcigar = self.cigar.as_mut().unwrap();
        lcigar.add_softclip(q_softclip);
        lcigar.0.extend_from_slice(&cigar.0[(r_offset + q_softclip)..]);
        self.reference_cigar_range.start = lr.1.start + r_offset;

        // eprintln!("q_dove {}, r_dove {}, qinsert {}, rdel {}\n -> {} {}", q_dove, r_dove, q_inserts, r_deletions, 0, r_offset);
        // eprintln!("Q {:?}  R {:?}", lr.0, lr.1);
        // eprintln!("{}\n{}\n{}",
        //     String::from_utf8_lossy(q),
        //     String::from_utf8_lossy(&r[r_offset..]),
        //     String::from_utf8_lossy(&cigar.0[r_offset..]));
        // print_alignment(q, &r[r_offset..], &lcigar.0);
        // eprintln!("E-----------------------------");

        (score, status, 0, &lr.1.start + r_offset)
    }

    pub fn cigar(&mut self) -> &mut Cigar {
        assert!(self.cigar.is_some());
        self.cigar.as_mut().unwrap()
    }

    pub fn align_right_flank(&mut self, aligner: &mut impl Align, query: &[u8], reference: &[u8], free_ends: usize) -> (i32, Status) {
        // Accurate alignment of flanks first.
        // Add threshold later and do hamming first, and if the score can possibly improve with perfect alignment, do that

        let mut rr = self.right_flank(query.len(), reference.len());

        if rr.0.len() == 0 || rr.1.len() == 0 {
            let q_softclip = rr.0.len();
            let lcigar = self.cigar.as_mut().unwrap();
            lcigar.add_softclip(q_softclip);
            self.reference_cigar_range.end = rr.1.end;
            return (0, Status::OK)
        }

        let q_end = min(rr.0.end + free_ends, query.len());
        let r_end = min(rr.1.end + free_ends, reference.len());
        let q_dove = q_end - rr.0.end;
        let r_dove = r_end - rr.1.end;

        // eprintln!("---RIGHT---\nqend = {}, rend = {} ..from.. q: {} r: {}", q_end, r_end, rr.0.end, rr.1.end);
        rr.0.end = q_end;
        rr.1.end = r_end;

        let q = &query[rr.0.clone()];
        let r = &reference[rr.1.clone()];

        // eprintln!("{}\n{}",
        //     String::from_utf8_lossy(q), 
        //     String::from_utf8_lossy(r));
        
        aligner.set_ends_free(0, q_dove as i32, 0,  r_dove as i32);
        // aligner.set_ends_free(100,100,100,100);

        let (score, cigar, status) = aligner.align(q, r);

        if !matches!(status, Status::OK) { return (std::i32::MIN, status) };

        // eprintln!("Score: {} (Free Q {}, Free R {})", score, q_dove, r_dove);
        // print_alignment(q, r, &cigar.0);

        let q_inserts = cigar.count_trailing_chars(b'I');
        let r_deletions = cigar.count_trailing_chars(b'D');
        let q_softclip: usize = min(r_deletions, q_dove);
        let r_offset: usize = min(q_inserts, r_dove);

        // Update cigar and set reference starting point
        let lcigar = self.cigar.as_mut().unwrap();
        lcigar.0.extend_from_slice(&cigar.0[0..cigar.0.len() - (r_offset + q_softclip)]);
        lcigar.add_softclip(q_softclip);
        self.reference_cigar_range.end = rr.1.end - r_offset;


        (score, status)
    }


    pub fn align_middle(&mut self, query: &[u8], reference: &[u8], max_score: &mut i32) -> AlignmentResult {
        let mut current_i = 0;
        let mut next_i = 1;

        let matches = (&self.seeds.first().unwrap()).length as usize;
        self.cigar().add_matches(matches);

        let mut score = 0;
        // eprintln!("Score align middle begin: {}", *max_score);
        while next_i < self.seeds.len() {
            let middle_range = self.between(&self.seeds[current_i], &self.seeds[next_i]);

            if middle_range.0.start > middle_range.0.end {
                return Err(AlignmentError::InvalidRangeError(format!("Invalid Range {:?}", middle_range.0)));
            }
            if middle_range.1.start > middle_range.1.end {
                return Err(AlignmentError::InvalidRangeError(format!("Invalid Range {:?}", middle_range.0)));
            }
            if middle_range.0.end >= query.len() {
                return Err(AlignmentError::QueryRangeError(format!("MRange {:?}.. Q len {} R len {}", middle_range, query.len(), reference.len())));
            }        
            if middle_range.1.end >= reference.len() {
                return Err(AlignmentError::ReferenceRangeError(format!("MRange {:?}.. Q len {} R len {}", middle_range, query.len(), reference.len())));
            }

            let middle_q = &query[middle_range.0.clone()];
            let middle_r = &reference[middle_range.1.clone()];

            let mut mismatches = 0;
            zip(middle_q, middle_r).for_each(|(q,r)| {
                self.cigar().0.push(if *q == *r { b'M' } else { mismatches += 1; b'X' });
            });
            *max_score -= mismatches * 4;
            score -= mismatches * 4;

            // eprintln!("Score align middle iter:  {} ... mismatches {}", *max_score, mismatches);

            if *max_score < 0 { 
                // eprintln!("Max score drop middle: {}", *max_score);
                return Ok((std::i32::MIN, Status::Dropped))
            };

            let matches = (&self.seeds[next_i]).length as usize;
            self.cigar().add_matches(matches);

            current_i += 1;
            next_i += 1;
        }

        return Ok((score, Status::OK));
    }

    pub fn extend_seeds(&mut self, query: &[u8], reference: &[u8]) {

        // Check orientation before this !
        if !self.orientation_set { 
            return
        }
        // To left
        let left_range = self.left_flank();

        if left_range.0.start + left_range.0.len() > query.len() {
            eprintln!("{}", self);
            self.visualize_alignment(query, reference);
            
            let valid_seeds = self.validate_seeds(query, reference);
            eprintln!("Seeds valid? {}", valid_seeds);

            panic!("Issues here {:?} {}", left_range, query.len())
        }

        /////////////////////////////////////////////////////////////////////////////////////
        // Left Flank
        /////////////////////////////////////////////////////////////////////////////////////
        let left_q = &query[left_range.0.clone()];
        let left_r = &reference[left_range.1.clone()];

        let by = zip(left_q, left_r)
            .rev()
            .enumerate()
            .find(|(i, (a, b))| a != b);
        
        match by {
            Some((by, _)) => self.seeds.first_mut().unwrap().extend_left(by),
            None =>  {
                self.seeds.first_mut().unwrap().extend_left(left_q.len())
            },
        }

        /////////////////////////////////////////////////////////////////////////////////////
        // Middle
        /////////////////////////////////////////////////////////////////////////////////////
        let mut current_i = 0;
        let mut next_i = 1;
        while next_i < self.seeds.len() {
            let middle_range = self.between(&self.seeds[current_i], &self.seeds[next_i]);
            
            let middle_q = &query[middle_range.0.clone()];
            let middle_r = &reference[middle_range.1.clone()];

            // eprintln!("Middle:\n{}\n{}", String::from_utf8_lossy(middle_q), String::from_utf8_lossy(middle_r));

            // First extend from right to left to see if we can merge
            let by = zip(middle_q, middle_r)
                .rev()
                .enumerate()
                .find(|(i, (a, b))| a != b);

            match by {
                Some((by, _)) => self.seeds[next_i].extend_left(by),
                None =>  {
                    // merge with righthand neighbor 
                    let right_len = self.seeds[next_i].length as usize;
                    self.seeds[current_i].extend_right(middle_q.len() + right_len);
                    self.seeds.remove(next_i);
                    continue
                },
            }

            // Cannot merge seeds.
            if by.is_some() {
                let by = zip(middle_q, middle_r)
                    .enumerate()
                    .find(|(i, (a, b))| a != b);

                match by {
                    Some((by, _)) => self.seeds[current_i].extend_right(by),
                    None => panic!("This should not happen."),
                }
            }

            current_i += 1;
            next_i += 1;
        }

        /////////////////////////////////////////////////////////////////////////////////////
        // Right Flank
        /////////////////////////////////////////////////////////////////////////////////////
        let right_range = self.right_flank(query.len(), reference.len());
        let right_q = &query[right_range.0.clone()];
        let right_r = &reference[right_range.1.clone()];

        let by = zip(right_q, right_r)
            .enumerate()
            .find(|(i, (a, b))| a != b);
        
        match by {
            Some((by, _)) => {
                // eprintln!("{}", String::from_utf8_lossy(right_q));
                // eprintln!("{}", String::from_utf8_lossy(right_r));
                // eprintln!("Before {:?}", self.seeds.last());
                self.seeds.last_mut().unwrap().extend_right(by);
                // eprintln!("After  {:?}", self.seeds.last());

                // let last = self.seeds.last_mut().unwrap();
                // if last.qend() > query.len() || last.rend() > reference.len() {
                //     panic!("By {:?} -- qe{} ql{} re{} rl{}", by, last.qend(), query.len(), last.rend(), reference.len())
                // }
            },
            None =>  {
                self.seeds.last_mut().unwrap().extend_right(right_q.len())
            },
        }


        // if by.is_some_and(|x| { x.0 > 10 }) || by.is_none() {
        //     eprintln!("-----\n{:?}, length of range {}, first seed: {:?}", right_range.0, right_q.len(), self.seeds.last().unwrap());
        //     eprintln!("{}", String::from_utf8_lossy(right_q));
        //     eprintln!("{}", String::from_utf8_lossy(right_r));
        //     eprintln!("Extend right by {:?}", by);
        // }
    }

    pub fn hamming(&self, query: &[u8], reference: &[u8]) -> u64 {
        let (qr, rr) = self.whole(query.len(), reference.len());
        triple_accel::hamming(&query[qr], &reference[rr]) as u64
    }

    pub fn gap_iter(&self) -> impl Iterator<Item = (Range<usize>, Range<usize>)> + '_ {
        zip(&self.seeds[1..], &self.seeds[0..self.seeds.len() - 1])
            .map(|(curr,next)| {
                ((curr.qend()..next.qbegin()), (curr.rend()..next.rbegin()))
        })
    }
    
    pub fn get_indel(&self, other: &Self, read_length: usize) -> i32 {
        let seed_self = &self.seeds[0];
        let seed_other = &other.seeds[0];
        let offset = if !self.orientation_set || !other.orientation_set {
            let offsets1 = seed_self.offsets(read_length);
            let offsets2 = seed_other.offsets(read_length);
            min(offsets1.0 - offsets2.0, offsets1.1 - offsets2.1)
        } else {
            let offset1 = seed_self.offset();
            let offset2 = seed_other.offset();
            (offset1 - offset2) as i64
        };
        
        // if offset.abs() < 10 {
        //     // eprintln!("Indel?! {}", offset);
        //     // eprintln!("{}\n{}", self, other);
        //     println!("Indel {} {} {}", offset, !self.orientation_set, !other.orientation_set);
        // }
        offset as i32
    }

    pub fn whole(&self, read_length: usize, ref_length: usize) -> (Range<usize>, Range<usize>) {
        //requires seeds sorted in ascending order
        let s: &AnchorSeed = self.seeds.first().unwrap();

        let q_overhang_length = s.qbegin();
        let r_overhang_length = s.rbegin();
        let left_overhang_length = min(q_overhang_length, r_overhang_length);

        let q_overhang_length = read_length - s.qend();
        let r_overhang_length = ref_length - s.rend();
        let right_overhang_length = min(q_overhang_length, r_overhang_length);

        (((s.qbegin() - left_overhang_length)..s.qend() + right_overhang_length),((s.rbegin() - left_overhang_length)..s.rend() + right_overhang_length))
    }

    pub fn left_flank(&self) -> (Range<usize>, Range<usize>) {
        //requires seeds sorted in ascending order
        let s: &AnchorSeed = self.seeds.first().unwrap();

        let q_overhang_length = s.qbegin();
        let r_overhang_length = s.rbegin();
        let overhang_length = min(q_overhang_length, r_overhang_length);

        (((s.qbegin() - overhang_length)..s.qbegin()),((s.rbegin() - overhang_length)..s.rbegin()))
    }

    pub fn right_flank(&self, read_length: usize, ref_length: usize) -> (Range<usize>, Range<usize>) {
        //requires seeds sorted in ascending order
        let s: &AnchorSeed = self.seeds.last().unwrap();

        let q_overhang_length = read_length - s.qend();
        let r_overhang_length = ref_length - s.rend();
        let overhang_length = min(q_overhang_length, r_overhang_length);
        
        ((s.qend()..s.qend() + overhang_length),(s.rend()..s.rend() + overhang_length))
    }

    pub fn between(&self, s1: &AnchorSeed, s2: &AnchorSeed) -> (Range<usize>, Range<usize>) {
        ((s1.qend()..s2.qbegin()), (s1.rend()..s2.rbegin()))
    }

    pub fn all_seeds_valid(&self, query: &[u8], reference: &[u8]) -> bool {
        self.seeds.iter().all(|s| {
            hamming(&query[s.qrange()], &reference[s.rrange()]) == 0
        })
    }

    pub fn validate_seeds(&self, query: &[u8], reference: &[u8]) -> bool {
        // println!("Validate: {}\nLENGTH: {}", self, self.seeds.len());

        self.seeds.iter().all(|s| {
            // println!("Seed: {}", s);
            // println!("Q: {}", String::from_utf8_lossy(&query[s.qrange()]));
            // println!("R: {}", String::from_utf8_lossy(&reference[s.rrange()]));

            hamming(&query[s.qrange()], &reference[s.rrange()]) == 0
        })
    }

    pub fn valid_seed_check(&self, query: &[u8], reference: &[u8]) {
        self.seeds.iter().for_each(|s| {
            println!("Seed: {}", s);
            println!("\tQ: {}", String::from_utf8_lossy(&query[s.qrange()]));
            println!("\tR: {}", String::from_utf8_lossy(&reference[s.rrange()]));
        });
    }

    pub fn valid_seed_count(&self, query: &[u8], reference: &[u8]) -> usize {
        self.seeds.iter().filter(|&s| {
            hamming(&query[s.qrange()], &reference[s.rrange()]) == 0
        }).count()
    }

    pub fn are_all_seeds_valid_any_config(&self, query: &[u8], query_rc: &[u8], reference: &[u8]) -> bool {
        self.seeds.iter().all(|fwd| {
            let mut rev = fwd.clone();
            rev.reverse(query.len());
            let rseed: &[u8] = &reference[fwd.rrange()];

            if self.forward {
                hamming(&query[fwd.qrange()], rseed) == 0 ||
                hamming(&query_rc[rev.qrange()], rseed) == 0 
            } else {
                hamming(&query[fwd.qrange()], rseed) == 0 ||
                hamming(&query_rc[rev.qrange()], rseed) == 0 
            }
        })
    }

    pub fn are_all_seeds_valid(&self, query: &[u8], reference: &[u8]) -> bool {
        self.seeds.iter().all(|seed| {
            hamming(&query[seed.qrange()], &reference[seed.rrange()]) == 0 
        })
    }

    pub fn any_orientation_valid(&mut self, rec: &RefFastqRecord, rec_rc: &OwnedFastqRecord, reference: &[u8]) -> bool {
        // println!("Any orientation valid");
        // println!("F REC   {}/{}", self.valid_seed_count(rec.seq(), reference), self.seeds.len());
        // println!("F RECRC {}/{}", self.valid_seed_count(rec_rc.seq(), reference), self.seeds.len());
        // self.valid_seed_check(rec.seq(), reference);
        // self.valid_seed_check(rec_rc.seq(), reference);
        if self.validate_seeds(rec.seq(), reference) {
            // eprintln!("YES 1");
            return true 
        };
        if self.validate_seeds(rec_rc.seq(), reference) {
            // eprintln!("YES 2");
            return true 
        };


        self.seeds.first_mut().unwrap().reverse(rec.seq().len());

        // println!("R REC   {}/{}", self.valid_seed_count(rec.seq(), reference), self.seeds.len());
        // println!("R RECRC {}/{}", self.valid_seed_count(rec_rc.seq(), reference), self.seeds.len());
        // self.valid_seed_check(rec.seq(), reference);
        // self.valid_seed_check(rec_rc.seq(), reference);
        if self.validate_seeds(rec.seq(), reference) {
            // eprintln!("YES 3");
            return true 
        };
        if self.validate_seeds(rec_rc.seq(), reference) {
            // eprintln!("YES 4");
            return true 
        };

        // println!("________________________________________{}", self.reference);

        return false
    }
    
    // pub fn resolve_orientation(&mut self, rec: &RefFastqRecord, rec_rc: &OwnedFastqRecord, reference: &[u8]) -> Result<()> {
    //     if self.validate_seeds(rec.seq(), reference) {
    //         //TODO: Code here is janky.
    //         return true 
    //     };
    //     if self.validate_seeds(rec_rc.seq(), reference) {
    //         eprintln!("Validate 2");
    //         return true 
    //     };
    //     self.seeds.first_mut().unwrap().reverse(rec.seq().len());
    //     if self.validate_seeds(rec.seq(), reference) {
    //         eprintln!("Validate 3");
    //         return true 
    //     };
    //     if self.validate_seeds(rec_rc.seq(), reference) {
    //         eprintln!("Validate 4");
    //         return true 
    //     };
        
    //     Error
    // }


    pub fn set_config(&mut self, config: &AnchorSeedConfig, read_length: usize) {
        type ASC = AnchorSeedConfig;
        match config {
            ASC::QueryRCSeedRC => {
                self.set_forward(false, read_length);
            },
            ASC::QuerySeedRC => {
                self.reverse_seeds(read_length);
                self.forward = true;
            },
            ASC::QueryRCSeed => {
                self.forward = false;
            },
            ASC::QuerySeed => {
                self.forward = true;
            },
            _ => {}
        }
    }

    pub fn set_forward(&mut self, forward: bool, read_length: usize) -> &Self {
        self.orientation_set = true;
        self.forward = forward;
        assert!(self.seeds.len() == 1);
        if !forward {
            self.seeds.first_mut().unwrap().reverse(read_length);
        }
        self
    }

    pub fn reverse_seeds(&mut self, read_length: usize) -> &Self {
        self.seeds.iter_mut().for_each(|s| {
            s.reverse(read_length);
        });
        self
    }

    pub fn visualize_alignment(&self, query: &[u8], reference: &[u8]) -> bool {
        let mut old: &AnchorSeed = self.seeds.first().unwrap();
        eprintln!("Q {} {}", query.len(), String::from_utf8_lossy(query));
        eprintln!("{:?}", old);
        let mut qspace = &query[0..old.qbegin()];
        let mut qseed = &query[old.qrange()];

        eprint!("Alignment Visualization:\n{}{}", qspace.ts().color(Color::Red), qseed.ts().color(Color::Green));

        self.seeds.iter().skip(1).enumerate().for_each(|(i, s)| {
            qspace = &query[old.qend()..s.qbegin()];
            qseed = &query[s.qrange()];
            eprint!("{}{}", qspace.ts().color(Color::Red), qseed.ts().color(Color::Green));

            old = s;
        });

        qspace = &query[old.qend()..];
        eprintln!("{}", qspace.ts().color(Color::Red));

        let mut old: &AnchorSeed = self.seeds.first().unwrap();

        let mut rstart = old.rbegin() - old.qbegin();
        if old.rbegin() < old.qbegin() {
            rstart = 0;
            eprint!("{}", String::from_utf8(vec![b' '; old.qbegin() - old.rbegin()]).unwrap());
        }
        let mut rspace = &reference[rstart..old.rbegin()];
        let mut rseed = &reference[old.rrange()];
        eprint!("{}{}", rspace.ts().color(Color::Red), rseed.ts().color(Color::Green));

        self.seeds.iter().skip(1).enumerate().for_each(|(i, s)| {
            rspace = &reference[old.rend()..s.rbegin()];
            rseed = &reference[s.rrange()];
            eprint!("{}{}", rspace.ts().color(Color::Red), rseed.ts().color(Color::Green));

            old = s;
        });
        
        rspace = &reference[old.rend()..min((old.rend() + qspace.len()), reference.len())];
        eprintln!("{}", rspace.ts().color(Color::Red));

        true
    }

    pub fn reference_pos(&self, read_length: usize) -> (u64, u64) {
        let seed = self.seeds.first().unwrap();
        let start = seed.rpos - seed.qpos as u64;
        (start, start + read_length as u64)
    }

    pub fn add_seed(&mut self, seed: &Seed, read_length: u32) {
        self.seed_count += 1;

        let s: &mut AnchorSeed = self.seeds.first_mut().unwrap();

        if s.qpos == seed.qpos && s.rpos == seed.rpos {
            let _ = s.clone();
            if s.length > seed.length as u32 {
                // eprintln!("----Replace\nFirst: qpos {}, rpos {}, len {}\nToAdd: {}", sc.qpos, sc.rpos, sc.length, seed.to_string());
                self.mismatches = seed.mismatch as u32;
                s.length = seed.length as u32;
            }
            if self.seeds.len() > 1 {
                panic!("Expected only one seed");
            }
            return
        }

        let mut aseed = AnchorSeed {
            qpos: seed.qpos,
            rpos: seed.rpos,
            length: seed.length as u32,
        };


        if !self.orientation_set {
            if aseed.contains(s) {
                eprintln!("Return1");
                s.set(&mut aseed);
                return
            }

            self.forward = seed.qpos > s.qpos && seed.rpos > s.rpos;
            eprintln!("Set direction ->> qpos {} rpos {} len {}\n{}\n--->  Forward? {}", s.qpos, s.rpos, s.length, seed.to_string(), self.forward);
            self.orientation_set = true;
            if !self.forward {
                s.reverse(read_length as usize);
            }
        }
        
        if !self.forward {
            aseed.reverse(read_length as usize);
        }

        if aseed.qpos < s.qpos {
            eprintln!("Return {}  {}", s, aseed);
            // eprintln!("\n\n\n-----\nAnchor: {} {} Size: {} ... {}", self.forward, self.forward_set, self.seed_count, self.seeds.len());
            // eprintln!("Anchor: {}", self.to_string());
            // eprintln!("Seed: {}", seed.to_string());
            // eprintln!("qpos {}, rpos {}, length {} ", aseed.qpos, aseed.rpos, aseed.length);
            // eprintln!("{} {}", self.seeds.first().unwrap().qpos, aseed.qpos);
            if self.orientation_set && self.seed_count == 1 {
                if !self.forward { s.reverse(read_length as usize); }
                self.orientation_set = false;
            }
            return
        }

        // Assume seeds come sorted by rpos. This makes the logic for merging seeds a lot easier.
        // After adding the second seed, orientation is clear. 
        assert!(aseed.qpos >= self.seeds.first().unwrap().qpos);
        match self.seeds.last_mut().unwrap().rpos_sorted_merge_into(&aseed) {
            SeedOverlap::NoOverlap => self.seeds.push(aseed),
            SeedOverlap::ContainedSelf => {},
            _ => {},
            // SeedOverlap::OffsetFwdOther => (),
            // SeedOverlap::OffsetFwdSelf => todo!(),
            // SeedOverlap::ContainedOther => todo!(),
        }
    }


    pub fn core_matches(&self) -> usize {
        self.seeds.iter().fold(0, |acc, seed| acc + seed.length as usize)
    }

    pub fn indels(&self) -> usize {
        if self.seeds.len() <= 1 { return 0 };
        self.seeds.iter()
            .zip(self.seeds.iter().skip(1))
            .fold(0, |acc, (seed1, seed2)| {
                acc + (seed2.qpos as usize - seed1.qpos as usize).abs_diff(seed2.rpos as usize - seed1.rpos as usize)
            })
    }
}

impl Display for Anchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let first = self.seeds.first().unwrap();
        let mut seeds_vstr = String::new();
        let mut seeds_str = String::new();

        for seed in &self.seeds {
            let seed_char = if self.orientation_set { 
                if self.forward { b'>' } else { b'<' }
            } else { b'X' };
            let spaces = String::from_utf8(vec![b' '; (seed.qpos) as usize]).unwrap();
            let xes = String::from_utf8(vec![seed_char; seed.length as usize]).unwrap();
            seeds_vstr += &spaces;
            seeds_vstr += &xes;
            seeds_vstr += "\n";

            seeds_str += &format!(" (qpos {} rpos {} len {})", seed.qpos, seed.rpos, seed.length);
        }

        write!(f, "{} --- Ref: {}, qpos: {}, rpos: {}, seed_count {}, mismatches: {}, core_matches: {} -- offset: {}\n{}\n{}",
            if self.orientation_set { 
                if self.forward { ">>>>>" } else { "<<<<<" }
            } else { "XXXXX" },
            self.reference,
            first.qpos,
            first.rpos,
            self.seed_count,
            self.mismatches,
            self.seeds.iter().fold(0, |acc, seed| acc + seed.length),
            first.rpos as i64 - first.qpos as i64,
            seeds_vstr,
            seeds_str)
    }
}

impl Default for Anchor {
    fn default() -> Self {
        Self { 
            reference: Default::default(),
            seed_count: Default::default(), 
            mismatches: Default::default(),
            forward: true,
            orientation_set: false,
            flagged_for_indel: false,
            flag: 0u8,
            counter1: 0,
            counter2: 0,
            score: 0,
            seeds: Vec::new(),
            cigar: None,
            reference_cigar_range: 0..0,
        }
    }
}

#[derive(Debug, Clone)]
pub enum AnchorSeedConfig {
    QuerySeed,
    QuerySeedRC,
    QueryRCSeed,
    QueryRCSeedRC,
    None,
}

pub fn get_seed_config(seed: &AnchorSeed, query: &[u8], query_rc: &[u8], reference: &[u8]) -> AnchorSeedConfig {

    fn seed_match(query_seed: &[u8], reference_seed: &[u8]) -> bool {
        hamming(query_seed, reference_seed) == 0
    }
    
    let mut seed_rc = seed.clone();
    seed_rc.reverse(query.len());

    let reference_seed = &reference[seed.rrange()];

    type ASC = AnchorSeedConfig;
    if seed_match(&query_rc[seed_rc.qrange()], reference_seed) {
        return ASC::QueryRCSeedRC;
    }
    if seed_match(&query[seed.qrange()], reference_seed) {
        return ASC::QuerySeed;
    }
    if seed_match(&query[seed_rc.qrange()], reference_seed) {
        return ASC::QuerySeedRC;
    }
    if seed_match(&query_rc[seed.qrange()], reference_seed) {
        return ASC::QueryRCSeed;
    }

    eprintln!("{}", hamming(&query_rc[seed_rc.qrange()], reference_seed));
    eprintln!("{}", hamming(&query[seed.qrange()], reference_seed));
    eprintln!("{}", hamming(&query[seed_rc.qrange()], reference_seed));
    eprintln!("{}", hamming(&query_rc[seed.qrange()], reference_seed));
    
    ASC::None
}

pub fn seed_match(seed: &AnchorSeed, query: &[u8], reference: &[u8]) -> bool {
    fn seed_match(query_seed: &[u8], reference_seed: &[u8]) -> bool {
        hamming(query_seed, reference_seed) == 0
    }
    seed_match(&query[seed.qrange()], &reference[seed.rrange()])
}

#[derive(Clone)]
pub struct Alignment {
    pub reference_id: u64,
    pub position: u32,
    pub forward: bool,
    pub cigar: Cigar,
}

pub type Alignments<'a> = &'a [Alignment];

impl Alignment {
    fn valid(&self) -> bool {
        true
    }
}