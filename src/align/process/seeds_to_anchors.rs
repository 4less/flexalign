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
                if let Some(_seed_match) = a.match_seed(seed, read_length) {
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
            // Clusters that merge across an indel become one anchor; every cluster that merges
            // with none is kept as it is. Keeping only the merges dropped all the other loci on
            // the reference -- on a chromosome, every other copy of a repeat, the true one among
            // them, whenever any two clusters anywhere on it merged.
            let mut merged = vec![false; anchor_tmp.len()];
            for (i, j) in (0..anchor_tmp.len()).tuple_combinations() {
                if let Some(mut a) = anchor_tmp[i].try_merge(&anchor_tmp[j], max_indel, read_length) {
                    a.0.flagged_for_indel = true;
                    anchors.push(a.0);
                    merged[i] = true;
                    merged[j] = true;
                }
            }
            if merged.iter().any(|&m| m) {
                anchors.extend(anchor_tmp.iter().zip(&merged).filter(|(_, &m)| !m).map(|(a, _)| a.clone()));
            }
            // eprintln!("Done trying merge....");
            
        } else {
            anchors.push(anchor_tmp[0].clone());
        }

        if anchors.is_empty() {
            anchors.extend_from_slice(&anchor_tmp.as_slice());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(qpos: u32, rpos: u64) -> Seed {
        Seed { rpos, rval: 0, qpos, mismatch: 0, length: 31, flag: 0 }
    }

    /// Two clusters that merge across an indel used to take every other locus on the reference
    /// with them: only the merge was kept.
    #[test]
    fn a_merge_keeps_the_other_loci_on_the_reference() {
        // Offsets 1000 and 1004 merge (a 4-base indel); 50_000 is another copy.
        let seeds = [seed(0, 1000), seed(60, 1064), seed(20, 50_020)];
        let mut anchors = Vec::new();
        GroupSeeds::group_into_anchor(&seeds, 150, &mut anchors);
        let starts: Vec<u64> = anchors.iter().flat_map(|a| a.seeds.iter().map(|s| s.rpos)).collect();
        assert!(anchors.iter().any(|a| a.flagged_for_indel), "the near clusters merge: {anchors:?}");
        assert!(starts.contains(&50_020), "the far locus survives the merge: {anchors:?}");
    }
}
