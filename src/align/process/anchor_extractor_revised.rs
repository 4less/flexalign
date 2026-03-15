use std::{
    cmp::{max, min},
    mem::swap,
};

use crate::{
    align::{
        common::{
            AnchorExtractor, AnchorPair, AnchorScore, PairedAnchorExtractor, SeedGroupPairList,
            SeedGroupPairedList, StdAnchorScore,
        },
        data_structures::{anchor::{Anchor, SeedGroupPaired}, common::Seed},
        stats::{Stats},
    },
    flexalign::time,
};


#[derive(Clone)]
pub struct StdAnchorExtractor {
    pub anchors: Vec<Anchor>,
    pub indices: Vec<usize>,
    pub other_indices: Vec<usize>,
    pub anchor_map: micromap::Map<u64, u32, 64>,
    pub groups: Vec<(u32, u32)>,
}

#[derive(Clone)]
pub struct StdPairedAnchorExtractor {
    pub groups: SeedGroupPairedList,
    pub groups_paired: SeedGroupPairList,

    pub indices: Vec<usize>,
    pub other_indices: Vec<usize>,

    pub anchors: Vec<AnchorPair>,
    pub anchors_fwd: Vec<Anchor>,
    pub anchors_rev: Vec<Anchor>,
}


#[inline(always)]
pub fn seed_group_indices_module(seeds: &[Seed], groups: &mut Vec<(u32, u32)>) -> usize {
    let mut last_idx = 0;
    // let mut groups = Vec::new();
    let mut max_size = 0;
    groups.clear();

    for i in 1..seeds.len() {
        let prev = &seeds[i - 1];
        let next = &seeds[i];
        if prev.rval != next.rval {
            groups.push((last_idx as u32, i as u32));
            max_size = max(max_size, i - last_idx);
            last_idx = i;
        }
    }

    groups.push((last_idx as u32, seeds.len() as u32));
    max_size = max(max_size, seeds.len() - last_idx);

    // let acc = self.groups.iter().fold(0, |acc, (start, end)| { acc + (end-start) });
    // if seeds.len() != acc as usize {
    //     panic!("{} {}", acc, seeds.len());
    // }

    max_size
}


#[inline(always)]
pub fn seed_group_indices_paired_module(
    seeds: &[Seed],
    groups: &mut SeedGroupPairedList,
    forward: bool,
) -> usize {
    if seeds.len() == 0 {
        return 0;
    };

    let mut last_idx = 0;
    // let mut groups = Vec::new();
    let mut max_size = 0;

    for i in 1..seeds.len() {
        let prev = &seeds[i - 1];
        let next = &seeds[i];
        if prev.rval != next.rval {
            let size = i - last_idx;

            groups.push(SeedGroupPaired {
                reference: prev.rval,
                _dummy: false,
                start: last_idx as u32,
                size: size as u16,
                forward: forward,
            });
            max_size = max(max_size, size);
            last_idx = i;
        }
    }

    groups.push(SeedGroupPaired {
        reference: seeds.last().unwrap().rval,
        _dummy: false,
        start: last_idx as u32,
        size: (seeds.len() - last_idx) as u16,
        forward: forward,
    });

    max(max_size, seeds.len() - last_idx)
}



pub fn group_into_anchor<'a>(
    seeds: &[Seed],
    read_length: usize,

) -> () 
{
    if seeds.is_empty() {
        return;
    }

    let mut anchor_tmp = Vec::new();

    let first = &seeds[0];    
    anchor_tmp.push(Anchor::from_seed(&first));

    for seed in &seeds[1..] {
        for a in &anchor_tmp {
            let (offset_first, fwd, indel_first) = a.closest_offset_seed(&seed, read_length);
            // let (offset_first, fwd, indel_first) = first.closest_offset(&seed, read_length);
            //TODO
        }
    }
}


#[inline(always)]
pub fn group_into_anchor_module<'a>(
    seeds_extern: &[Seed],
    start: usize,
    end: usize,
    read_length: usize,
    indices: &'a mut Vec<usize>,
    other_indices: &'a mut Vec<usize>,
    anchors: &mut Vec<Anchor>,
) {
    let seeds = &seeds_extern[start..end];

    // eprintln!("-- Seed group -- {} ... {} - {}", seeds.first().unwrap().rval, start, end);
    // for seed in seeds.iter() {
    //     eprintln!("{} --- {}", seed, seed.reverse(read_length));
    // }

    let _anchor_group_index = anchors.len();
    let mut _added_anchors = 0;
    let mut indel_flag = false;

    if !seeds.is_empty() {
        // Group by exact offset. If there are seeds left, distribute them onto

        indices.clear();
        indices.extend(0..seeds.len());
        other_indices.clear();

        while !indices.is_empty() {
            let first = seeds.get(*indices.first().unwrap()).unwrap();

            // Set first seed as an anchor.
            let mut a: Anchor = Anchor::from_seed(first);

            let mut offset = None;
            let mut forward = None;

            // Find seeds that have 0 indels with respect to the first anchor
            for index in indices.iter().skip(1) {
                let next: &Seed = &seeds[*index];
                let (offset_first, fwd, indel_first) = first.closest_offset(&next, read_length);

                if indel_first == 0 && forward.is_none_or(|val| val == fwd) {
                    if offset.is_none() {
                        offset = Some(offset_first);
                        a.set_forward(fwd, read_length);
                        forward = Some(fwd);
                    }
                    a.add_seed(&next, read_length as u32);
                } else {
                    if indel_first > 0 && indel_first < 10 {
                        indel_flag = true;
                        a.flagged_for_indel = true;
                        a.add_seed(&next, read_length as u32);
                        // eprintln!("INDICATION FOR INDEL: {}", indel_first);
                    } else {
                        other_indices.push(*index);
                    }
                    // eprintln!("No Samey: \n{} \nOffset {}  Fed? {}  IndelFirst {}", next, offset_first, fwd, indel_first);
                }
            }

            anchors.push(a);
            _added_anchors += 1;

            indices.clear();
            swap(other_indices, indices);
        }

        // let mut local_anchors = &anchors[anchor_group_index..];
        // for i in 0..local_anchors.len() {
        //     for j in i+1..local_anchors.len() {
        //         let a = &local_anchors[i];
        //         let b = &local_anchors[j];

        //         let offset_dist = a.get_indel(b, read_length);
        //         if offset_dist.abs() < 10 {
        //             indel_flag = true;
        //         }
        //     }
        // }
    }
    // eprintln!("END -- Seed group -- {}", seeds.first().unwrap().rval );

    if indel_flag {
        println!("Indel")
    }
}

impl StdAnchorExtractor {
    pub fn new() -> Self {
        Self {
            anchors: Vec::new(),
            indices: Vec::new(),
            other_indices: Vec::new(),
            anchor_map: micromap::Map::default(),
            groups: Vec::new(),
        }
    }

    // pub fn seed_group_indices(&mut self, seeds: &[Seed]) -> usize {
    //     let mut last_idx = 0;
    //     // let mut groups = Vec::new();
    //     let mut max_size = 0;
    //     self.groups.clear();

    //     for i in 1..seeds.len() {
    //         let prev = &seeds[i-1];
    //         let next = &seeds[i];
    //         if prev.rval != next.rval {
    //             self.groups.push((last_idx as u32, i as u32));
    //             max_size = max(max_size, i - last_idx);
    //             last_idx = i;
    //         }
    //     }

    //     self.groups.push((last_idx as u32, seeds.len() as u32));
    //     max_size = max(max_size, seeds.len() - last_idx);

    //     // let acc = self.groups.iter().fold(0, |acc, (start, end)| { acc + (end-start) });
    //     // if seeds.len() != acc as usize {
    //     //     panic!("{} {}", acc, seeds.len());
    //     // }

    //     max_size
    // }

    pub fn group_into_anchor(
        &mut self,
        seeds_extern: &[Seed],
        start: usize,
        end: usize,
        read_length: usize,
    ) {
        let seeds = &seeds_extern[start..end];

        let _anchor_group_index = self.anchors.len();

        if seeds.is_empty() {
            return;
        };

        // Group by exact offset. If there are seeds left, distribute them onto

        self.indices.clear();
        self.indices.extend(0..seeds.len());
        self.other_indices.clear();

        while !self.indices.is_empty() {
            let first = seeds.get(self.indices[0]).unwrap();

            // Set first seed as an anchor.
            let mut a: Anchor = Anchor::from_seed(first);

            // eprintln!("Starting with anchor: {}", first.to_string());

            let mut offset = None;
            let mut forward = None;

            // Find seeds that have 0 indels with respect to the first anchor
            for index in self.indices.iter().skip(1) {
                let next: &Seed = &seeds[*index];

                let (offset_first, fwd, indel_first) = first.closest_offset(&next, read_length);

                if indel_first == 0 && (forward.is_none() || forward.unwrap() == fwd) {
                    if offset.is_none() {
                        // eprintln!("Establish offset: {} {} {}", offset_first, fwd, indel_first);
                        // eprintln!("Between: {} {}", );
                        offset = Some(offset_first);
                        a.set_forward(fwd, read_length);
                        forward = Some(fwd);
                    }
                    a.add_seed(&next, read_length as u32);
                } else {
                    self.other_indices.push(*index);
                }
            }
            self.anchors.push(a);
            self.indices.clear();
            swap(&mut self.other_indices, &mut self.indices);
        }
    }
}

impl AnchorExtractor for StdAnchorExtractor {
    fn generate(&mut self, seeds: &[Seed], read_length: usize, stats: &mut Stats) -> &mut [Anchor] {
        self.anchors.clear();
        self.anchor_map.clear();

        // let max_size = self.seed_group_indices(seeds);
        //TODO: Revisit and check function
        let max_size = seed_group_indices_module(seeds, &mut self.groups);

        let (_duration, _) = time(|| {
            glidesort::sort_by_key(&mut self.groups, |(start, end)| {
                -1i32 * (end - start) as i32
            });
        });

        let skip_threshold = min(max_size as i32, 32i32) - 10;

        stats.anchors += self.groups.len();
        for i in 0..min(8, self.groups.len()) {
            let (start, end) = self.groups[i];
            let group_size = end - start;
            // eprintln!("{} < {} == {} ({})", end-start, max_size - 5, (end - start) < (max_size as u32 - 5), max_size);
            if (group_size as i32) < skip_threshold {
                // eprintln!("Skip {} {}, {}, {},  {}", start, end, end-start, self.options.args.ranges, skip_threshold);
                continue;
            };

            // self.group_into_anchor(seeds, start as usize, end as usize, read_length);
            group_into_anchor_module(
                seeds,
                start as usize,
                end as usize,
                read_length,
                &mut self.indices,
                &mut self.other_indices,
                &mut self.anchors,
            );
        }

        &mut self.anchors
    }

    fn retrieve(&self) -> &[Anchor] {
        &self.anchors
    }

    fn retrieve_mut(&mut self) -> &mut [Anchor] {
        &mut self.anchors
    }
}

impl StdPairedAnchorExtractor {
    pub fn new() -> Self {
        Self {
            anchors_fwd: Vec::new(),
            anchors_rev: Vec::new(),
            anchors: Vec::new(),
            indices: Vec::new(),
            other_indices: Vec::new(),
            groups: Vec::new(),
            groups_paired: Vec::new(),
        }
    }
}

pub fn insert_size(
    a_fwd: Option<&Anchor>,
    a_rev: Option<&Anchor>,
    read_length_fwd: usize,
    read_length_rev: usize,
) -> Option<i64> {
    if a_fwd.is_none() || a_rev.is_none() {
        return None;
    };

    let span_fwd = a_fwd.as_ref().unwrap().reference_pos(read_length_fwd);
    let span_rev = a_rev.as_ref().unwrap().reference_pos(read_length_rev);

    Some(if span_fwd.0 < span_rev.0 {
        span_rev.0 as i64 - span_fwd.1 as i64
    } else {
        span_fwd.0 as i64 - span_rev.1 as i64
    })
}

// pub fn pair_anchors(anchors_fwd: &Vec<Anchor>, anchors_rev: &Vec<Anchor>, anchor_pairs: &Vec<AnchorPair>, max_dist: usize, read_length: usize) {
//     let _ = max_dist;

//     for a_fwd in anchors_fwd {
//         let span_fwd = a_fwd.reference_pos(read_length);
//         for a_rev in anchors_rev {
//             let span_rev = a_rev.reference_pos(read_length);

//             let insert_size = if span_fwd.0 < span_rev.0 { span_rev.0 as i64 - span_fwd.1 as i64 } else { span_fwd.0 as i64 - span_rev.1 as i64 };
//         }
//     }
// }

impl PairedAnchorExtractor for StdPairedAnchorExtractor {
    fn generate(
        &mut self,
        seeds_fwd: &[Seed],
        seeds_rev: &[Seed],
        read_length_fwd: usize,
        read_length_rev: usize,
        stats: &mut Stats,
    ) -> &mut [AnchorPair] {
        let _ = stats;

        self.groups.clear();
        self.groups_paired.clear();
        self.anchors_fwd.clear();
        self.anchors_rev.clear();
        self.anchors.clear();

        seed_group_indices_paired_module(seeds_fwd, &mut self.groups, true);
        let _fwd_size = self.groups.len();
        seed_group_indices_paired_module(seeds_rev, &mut self.groups, false);

        // eprintln!("Anchors (Fwd -> Rev): {} -> {}", fwd_size, self.groups.len());

        glidesort::sort_by_key(&mut self.groups, |e| (e.reference, e.forward));

        let mut current_idx = 0;
        let mut next_idx;

        if self.groups.is_empty() {
            return &mut self.anchors;
        }

        // eprintln!("ANCHOR PROCEDURE\n{:?}", self.groups);

        while current_idx < self.groups.len() - 1 {
            next_idx = current_idx + 1;
            let current = &self.groups[current_idx];
            let next = &self.groups[next_idx];

            // eprintln!("Iterate .. Current {}, Next {}", current, next);

            self.anchors_fwd.clear();
            self.anchors_rev.clear();

            if current.reference == next.reference {
                // eprintln!("DEAL WITH REFERENCE:  {}", current.reference);
                // process single(

                assert!(!current.forward);
                assert!(next.forward);

                group_into_anchor_module(
                    seeds_rev,
                    current.start as usize,
                    current.start as usize + current.size as usize,
                    read_length_rev,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_rev,
                );
                group_into_anchor_module(
                    seeds_fwd,
                    next.start as usize,
                    next.start as usize + next.size as usize,
                    read_length_fwd,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_fwd,
                );


                if self.anchors_fwd.len() <= 1 && self.anchors_rev.len() <= 1 {
                    self.anchors
                        .push(AnchorPair(self.anchors_fwd.pop(), self.anchors_rev.pop()));

                    let a = *self.anchors.last().as_ref().unwrap();
                    let _a_fwd = a.0.as_ref().unwrap();
                    let _a_rev = a.1.as_ref().unwrap();

                    // eprintln!("Insert size: {:?}   {:?}, {:?}", insert_size(Some(a_fwd), Some(a_rev), read_length), a_fwd.reference_pos(read_length), a_rev.reference_pos(read_length));
                } else {
                    for a_fwd in &self.anchors_fwd {
                        for a_rev in &self.anchors_rev {
                            match insert_size(
                                Some(a_fwd),
                                Some(a_rev),
                                read_length_fwd,
                                read_length_rev,
                            ) {
                                Some(is) => {
                                    if is < 1000 {
                                        self.anchors.push(AnchorPair(
                                            Some(a_fwd.clone()),
                                            Some(a_rev.clone()),
                                        ));
                                        self.anchors
                                            .last_mut()
                                            .unwrap()
                                            .resolve_orientation(read_length_fwd, read_length_rev);
                                    }
                                }
                                None => panic!(
                                    "This if branch is only entered if anchor has both reads"
                                ),
                            };
                        }
                    }
                }
                current_idx += 2;
            } else if current.forward {
                // eprintln!("SINGLE FWD DEAL WITH REFERENCE:  {}", current.reference);
                group_into_anchor_module(
                    seeds_fwd,
                    current.start as usize,
                    current.start as usize + current.size as usize,
                    read_length_fwd,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_fwd,
                );

                while !self.anchors_fwd.is_empty() {
                    self.anchors.push(AnchorPair(self.anchors_fwd.pop(), None));
                }
                current_idx += 1;
            } else {
                // eprintln!("SINGLE REV DEAL WITH REFERENCE:  {}", current.reference);
                group_into_anchor_module(
                    seeds_rev,
                    current.start as usize,
                    current.start as usize + current.size as usize,
                    read_length_rev,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_rev,
                );
                while !self.anchors_rev.is_empty() {
                    self.anchors.push(AnchorPair(None, self.anchors_rev.pop()));
                }
                current_idx += 1;
            }
        }

        if current_idx < self.groups.len() {
            let current = &self.groups[current_idx];
            if current.forward {
                // eprintln!("Extra SINGLE FWD DEAL WITH REFERENCE:  {}", current.reference);
                group_into_anchor_module(
                    seeds_fwd,
                    current.start as usize,
                    current.start as usize + current.size as usize,
                    read_length_fwd,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_fwd,
                );
    
                while !self.anchors_fwd.is_empty() {
                    self.anchors.push(AnchorPair(self.anchors_fwd.pop(), None));
                }
                current_idx += 1;
            } else {
                // eprintln!("Extra SINGLE REV DEAL WITH REFERENCE:  {}", current.reference);
                group_into_anchor_module(
                    seeds_rev,
                    current.start as usize,
                    current.start as usize + current.size as usize,
                    read_length_rev,
                    &mut self.indices,
                    &mut self.other_indices,
                    &mut self.anchors_rev,
                );
                while !self.anchors_rev.is_empty() {
                    self.anchors.push(AnchorPair(None, self.anchors_rev.pop()));
                }
                current_idx += 1;
            }
        }

        glidesort::sort_by_key(&mut self.anchors, |AnchorPair(a_fwd, a_rev)| {
            let s1 = a_fwd.as_ref().map_or(0, |a| StdAnchorScore::score(a));
            let s2 = a_rev.as_ref().map_or(0, |a| StdAnchorScore::score(a));

            -(s1 + s2)
        });

        &mut self.anchors
    }

    fn retrieve(&self) -> &[AnchorPair] {
        &self.anchors
    }

    fn retrieve_mut(&mut self) -> &mut [AnchorPair] {
        &mut self.anchors
    }
}
