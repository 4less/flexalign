use std::{cmp::{max, min}, fmt::Display, iter::{zip}, ops::Range};

use bioreader::sequence::{fastq_record::{OwnedFastqRecord, RefFastqRecord}};
use thiserror::Error;

use crate::align::{common::{Align, Heuristic, Status}, errors::{AlignmentError, AlignmentResult}, sam::Cigar};

use super::anchor::AnchorSeed;


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

    pub fn offset_with_other(&self, other: &Self, read_length: usize) -> ((i64, bool, u64), (i64, bool, u64)) {
        let (self_fwd, self_rev) = self.offsets(read_length);
        let (other_fwd, other_rev) = other.offsets(read_length);

        let diff_fwd = self_fwd.abs_diff(other_fwd);
        let diff_rev = self_rev.abs_diff(other_rev);

        (
            (self_fwd, true, diff_fwd),
            (self_rev, true, diff_rev),
        )
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

    pub fn query_range(&self) -> Range<usize> {
        (self.qpos as usize)..(self.qpos as usize + self.length as usize)
    }

    pub fn reference_range(&self) -> Range<usize> {
        (self.rpos as usize)..(self.rpos as usize + self.length as usize)
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