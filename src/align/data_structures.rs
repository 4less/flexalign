use std::{cmp::min, fmt::Display};


#[derive(Clone, Debug)]
pub struct Seed {
    pub qpos: u32,
    pub rpos: u64,
    pub rval: u64,
    pub mismatch: u8,
    pub length: u8,
}

impl Seed {
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

}


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Anchor {
    pub reference: u64, //
    pub seed_count: u32,
    pub mismatches: u32, //
    pub forward: bool,
    pub forward_set: bool,
    pub seeds: Vec<AnchorSeed>,
}

impl Anchor {
    pub fn set_forward(&mut self, forward: bool, read_length: usize) {
        self.forward_set = true;
        self.forward = forward;
        if !forward {
            self.seeds.first_mut().unwrap().reverse(read_length);
        }
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


        if !self.forward_set {
            if aseed.contains(s) {
                eprintln!("Return1");
                s.set(&mut aseed);
                return
            }

            self.forward = seed.qpos > s.qpos && seed.rpos > s.rpos;
            eprintln!("Set direction ->> qpos {} rpos {} len {}\n{}\n--->  Forward? {}", s.qpos, s.rpos, s.length, seed.to_string(), self.forward);
            self.forward_set = true;
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
            if self.forward_set && self.seed_count == 1 {
                if !self.forward { s.reverse(read_length as usize); }
                self.forward_set = false;
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
            let seed_char = if self.forward_set { 
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
            if self.forward_set { 
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
            forward_set: false,
            seeds: Vec::new(),
        }
    }
}