use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{align::{common::{AnchorPair, PairedAnchorSorter}, data_structures::anchor::{get_seed_config, Anchor, AnchorSeed, AnchorSeedConfig}, stats::Stats}, database::common::FlexalignDatabase};

pub type FixAnchorResult = Result<(), FixAnchorError>;

#[derive(thiserror::Error, Debug)] 
pub enum FixAnchorError {
    #[error("{0}")]
    UnresolvableOrientation(String),
}

#[derive(Clone)]
pub struct PairedAnchorHeuristicSorter<'a, D: FlexalignDatabase> {
    pub db: &'a D,
}

impl<'a, D: FlexalignDatabase> PairedAnchorHeuristicSorter<'a, D> {
    pub fn new(db: &'a D) -> Self {
        Self { db }
    }

    pub fn fix_anchor(a: &mut Anchor, query: &[u8], query_rc: &[u8], reference: &[u8]) -> FixAnchorResult {
        // println!("\n##################################################################################################\n___QueryLength: {}\n____ BEFORE {}", query.len(), a);
        let v = a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference);

        let first: &AnchorSeed = &a.seeds[0];
        assert!(first.qend() <= query.len());

        if !v {// initial configuration is incorrect
            assert!(a.orientation_set || a.seeds.len() <= 1); 

            let first_seed_config = get_seed_config(a.seeds.first().unwrap(), query, query_rc, reference);
            // let v = a.are_all_seeds_valid(if a.forward { rec_fwd } else { rec_fwd_revc }, reference);
            
            let first: &AnchorSeed = &a.seeds[0];
            assert!(first.qend() <= query.len());
            
            type ASC = AnchorSeedConfig;
            match &first_seed_config {
                ASC::None => {
                    // This means during the anchor building phase, two seeds must have been merged that actually do not work together.
                    // This can happend for k-mers that appear both as their regular and their reverse complement in a single query.
                    // let any = a.seeds.iter().any(|s| matches!(get_seed_config(s, query, query_rc, reference), ASC::None));
                    let index = a.seeds.iter().position(|s| matches!(get_seed_config(s, query, query_rc, reference), ASC::None));
                    match index {
                        Some(index) => {
                            let new_seed = a.seeds[index].clone();
                            a.seeds.clear();
                            a.seeds.push(new_seed.clone());

                            let config = get_seed_config(&new_seed, query, query_rc, reference);
                            a.set_config(&config, query.len());

                            let first: &AnchorSeed = &a.seeds[0];
                            assert!(first.qend() <= query.len());
                        },
                        None => panic!("Nothing correct?"),
                    }
                    
                },
                config => a.set_config(config, query.len()),
            }
            let first: &AnchorSeed = &a.seeds[0];
            assert!(first.qend() <= query.len());
            
            let v = a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference);
                        
            if !v {
                let _ = a.seeds.split_off(1);
                assert!(a.seeds.len() == 1);
                
                let first: &AnchorSeed = &a.seeds[0];
                assert!(first.qend() <= query.len());
                // assert!(a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference));
            } else {
                let first: &AnchorSeed = &a.seeds[0];
                assert!(first.qend() <= query.len());
            }

            let first: &AnchorSeed = &a.seeds[0];
            assert!(first.qend() <= query.len());
        } else {
            let first: &AnchorSeed = &a.seeds[0];
            assert!(first.qend() <= query.len());
        }


        if !a.seeds.as_slice().windows(2).all(|w: &[AnchorSeed]| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()) {
            return Err(FixAnchorError::UnresolvableOrientation(format!("Cannot resolve orientation of reads")))
        } else {
            a.seed_status.set_valid();
        }

        let first: &AnchorSeed = &a.seeds[0];
        assert!(first.qend() <= query.len());

        Ok(())
    }
}

impl<'a, D: FlexalignDatabase> PairedAnchorSorter for PairedAnchorHeuristicSorter<'a, D> {
    fn sort(&self, mut anchors: &mut [AnchorPair], 
            rec_fwd: &RefFastqRecord, rec_fwd_revc: &OwnedFastqRecord,
            rec_rev: &RefFastqRecord, rec_rev_revc: &OwnedFastqRecord, stats: &mut Stats) {
        let _ = stats;

        anchors.iter_mut().for_each(|AnchorPair(a1, a2)| {
            let reference: &&[u8] = match a1 {
                Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
                None => &self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
            };

            // Treat each anchor in three stages.
            // 1. Is initial configuration correct?
            // 2. Is any configuration correct for all seeds?
            // 3. Troubleshooting - there are mixed seeds for this anchor.
            match a1 {
                Some(a) => {
                    let query = if a.forward { rec_fwd.seq() } else { rec_fwd_revc.seq() };
                    let first: &AnchorSeed = &a.seeds[0];
                    assert!(first.qend() <= query.len());

                    let fix_result = Self::fix_anchor(a, rec_fwd.seq(), rec_fwd_revc.seq(), reference);

                    let query = if a.forward { rec_fwd.seq() } else { rec_fwd_revc.seq() };

                    let first: &AnchorSeed = &a.seeds[0];
                    assert!(first.qend() <= query.len());

                    a.score = (query.len() as u64 - a.hamming(query, reference)) as i32;

                    if fix_result.is_ok() && !a.seeds.as_slice().windows(2).all(|w: &[AnchorSeed]| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()) {
                        panic!("{}", a);
                    }

                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    let query = if a.forward { rec_rev.seq() } else { rec_rev_revc.seq() };
                    // println!("\n___QueryLength: {}\n____ AFTER {}", query.len(), a);
                    let first: &AnchorSeed = &a.seeds[0];
                    // println!("\n --- {}, {} -- {}", query.len(), first.qbegin(), first.qend());
                    assert!(first.qend() <= query.len());

                    let fix_result = Self::fix_anchor(a, rec_rev.seq(), rec_rev_revc.seq(), reference);
                    let query = if a.forward { rec_rev.seq() } else { rec_rev_revc.seq() };

                    // println!("\n___QueryLength: {}\n____ AFTER {}", query.len(), a);
                    let first: &AnchorSeed = &a.seeds[0];
                    assert!(first.qend() <= query.len());

                    a.score = (query.len() as u64 - a.hamming(query, reference)) as i32;

                    if fix_result.is_ok() && !a.seeds.as_slice().windows(2).all(|w: &[AnchorSeed]| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()) {
                        panic!("{}", a);
                    }

                }, _ => {},
            }                        

        });


        glidesort::sort_by_key(&mut anchors,|AnchorPair(a1, a2)| {
            let s1 = match a1 {
                Some(a) => a.score,
                None => 0,
            };
            let s2 = match a2 {
                Some(a) => a.score,
                None => 0,
            };

            - ((s1 + s2) as i64)
        });


        anchors.iter().for_each(|AnchorPair(a1, a2)| {
            match a1 {
                Some(a) => {
                    if a.seed_status.is_valid() && a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("U 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seed_status.is_valid() && a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("U 2  {}", a);
                    }
                }, _ => {},
            }
        });
    }
}

