use std::process::exit;

use itertools::Itertools;

use crate::align::data_structures::{anchor::Anchor, common::Seed};

pub trait SeedsToAnchors {
    fn group_into_anchor<'a>(seeds: &[Seed], read_length: usize, anchors: &mut Vec<Anchor>) -> ();
}

#[derive(Clone)]
pub struct GroupSeeds;
impl SeedsToAnchors for GroupSeeds {
    fn group_into_anchor<'a>(seeds: &[Seed], read_length: usize, anchors: &mut Vec<Anchor>) -> () {
        if seeds.is_empty() {
            return;
        }

        let mut anchor_tmp = Vec::new();

        let first = &seeds[0];

        
        anchor_tmp.push(Anchor::from_seed(&first));

        // eprintln!("--------------------------------------\nStarting seed: {}\n", first);
        for seed in &seeds[1..] {
            // eprintln!("####### Place Seed: {}", seed);
            let has_match = anchor_tmp.iter_mut().any(|a| {
                if let Some(seed_match) = a.match_seed(seed, read_length) {
                    // eprintln!("---------\nCorrect:\nSeed offsets: {:?}\nAnchor offsets: {:?}", seed.offsets(read_length), a.query_reference_alignment());
                    a.add_seed(seed, read_length as u32);
                    return true;
                }
                // eprintln!("---------\nFalse:\nSeed offsets: {:?}\nAnchor offsets: {:?}", seed.offsets(read_length), a.query_reference_alignment());
                // eprintln!("Match with first seed: \n{:?}", seed.offset_with_other(first, read_length));
                false
            });

            if !has_match {
                anchor_tmp.push(Anchor::from_seed(seed));
            }
        }

        let max_indel = 10;
        if anchor_tmp.len() > 1 {
            // eprintln!("Seeds group into multiple anchors: ");
            // for (idx, a) in anchor_tmp.iter().enumerate() {
            //     eprintln!("{} {}", idx, a);
            // }

            // eprintln!("Try merge....");
            anchor_tmp.iter().tuple_combinations().for_each(|(a, b)| {
                let mut new_anchor = a.try_merge(b, max_indel, read_length);

                if let Some(mut a ) = new_anchor {
                    a.0.flagged_for_indel = true;
                    anchors.push(a.0);
                }
            });
            // eprintln!("Done trying merge....");
            
        } else {
            anchors.push(anchor_tmp[0].clone());
        }

        if anchors.is_empty() {
            anchors.extend_from_slice(&anchor_tmp.as_slice());
        }
    }
}