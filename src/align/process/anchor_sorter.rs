use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{align::{common::{AnchorPair, PairedAnchorSorter}, data_structures::{get_seed_config, Anchor, AnchorSeedConfig}, stats::Stats}, database::common::FlexalignDatabase};

#[derive(Clone)]
pub struct PairedAnchorHeuristicSorter<'a, D: FlexalignDatabase> {
    pub db: &'a D,
}

impl<'a, D: FlexalignDatabase> PairedAnchorHeuristicSorter<'a, D> {
    pub fn new(db: &'a D) -> Self {
        Self { db }
    }

    pub fn fix_anchor(a: &mut Anchor, query: &[u8], query_rc: &[u8], reference: &[u8]) -> () {
        let v = a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference);

        if !v {// initial configuration is incorrect
            assert!(a.orientation_set || a.seeds.len() <= 1); 

            let first_seed_config = get_seed_config(a.seeds.first().unwrap(), query, query_rc, reference);
            // let v = a.are_all_seeds_valid(if a.forward { rec_fwd } else { rec_fwd_revc }, reference);
            
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
                        },
                        None => panic!("Nothing correct?"),
                    }
                    
                },
                config => a.set_config(&config, query.len()),
            }
            
            let v = a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference);
                        
            if !v {
                println!("\n_fix anchor_Initial {:?} ... Orientation Forward? {}", first_seed_config, a.forward);
                println!("\nAnchor {:?}", a);
                for s in a.seeds.iter() {
                    println!("{:?} <- {}", get_seed_config(s, query, query_rc, reference), s);
                }
                let _ = a.seeds.split_off(1);
                assert!(a.seeds.len() == 1);
                assert!(a.are_all_seeds_valid(if a.forward { query } else { query_rc }, reference));
            }
        }
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

            match a1 {
                // Treat each anchor in three stages.
                // 1. Is initial configuration correct?
                // 2. Is any configuration correct for all seeds?
                // 3. Troubleshooting - there are mixed seeds for this anchor.
                Some(a) => Self::fix_anchor(a, rec_fwd.seq(), rec_fwd_revc.seq(), reference)
                , _ => {},
            }

            match a2 {
                // Treat each anchor in three stages.
                // 1. Is initial configuration correct?
                // 2. Is any configuration correct for all seeds?
                // 3. Troubleshooting - there are mixed seeds for this anchor.
                Some(a) => Self::fix_anchor(a, rec_rev.seq(), rec_rev_revc.seq(), reference)
                , _ => {},
            }
        });

    
        anchors.iter_mut().for_each(|AnchorPair(a1, a2)| {


            // println!("----------\nBEGIN--/1 {:?}", a1);
            // println!("BEGIN--/2 {:?}", a2);
            // match a1 {
            //     Some(a) => {
            //         if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
            //             eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
            //             panic!("X 1  {}", a);
            //         }
            //     }, _ => {},
            // }

            // match a2 {            //     Some(a) => {
            //         if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
            //             eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
            //             panic!("X 2  {}", a);
            //         }
            //     }, _ => {},
            // }

            let reference: &&[u8] = match a1 {
                Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
                None => &self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
            };

            if !a1.as_ref().is_some_and(|s| s.orientation_set) || !a2.as_ref().is_some_and(|s| s.orientation_set) {
                let _a1_valid = match a1 {
                    Some(a) => {
                        // println!("_______________________________________________FWD");
                        a.any_orientation_valid(rec_fwd, rec_fwd_revc, reference)
                    },
                    None => true,
                };
                let _a2_valid = match a2 {
                    Some(a) => {
                        // println!("_______________________________________________REV");
                        a.any_orientation_valid(rec_rev, rec_rev_revc, reference)
                    },
                    None => true,
                };
                // eprintln!("Orientation not set. {} {}", a1_valid, a2_valid);
                return ()
            }

            match a1 {
                Some(a) => {
                    let query = if a.forward { rec_fwd.seq() } else { rec_fwd_revc.seq() };
                    if query.len() == 0 { 
                        a.score = 0i32;
                    } else {

                        a.extend_seeds(query, reference);
                        a.score = a.core_matches() as i32;                      
                        a.score = (query.len() as u64 - a.hamming(query, reference)) as i32;
                        // eprintln!("Set score {}", a.score);

                        if !a.seeds.as_slice().windows(2).all(|w: &[crate::align::data_structures::AnchorSeed]| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()) {
                            panic!("{}", a);
                        }

                        if a.flagged_for_indel {
                            eprintln!("Heyu");
                        }
                    }
                    // eprintln!("{}", query.len());


                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        panic!("B 1  {}", a);
                    }
                },
                None => (),
            };
            match a2 {
                Some(a) => {
                    // println!("1 Extended {}", a);
                    let query = if a.forward { rec_rev.seq() } else { rec_rev_revc.seq() };
                    if query.len() == 0 { 
                        a.score = 0i32;
                    } else {
                        a.extend_seeds(query, reference);
                        a.score = a.core_matches() as i32;  
                        a.score = (query.len() as u64 - a.hamming(query, reference)) as i32;

                        // println!("2 Extended {}", a);

                        if !a.seeds.as_slice().windows(2).all(|w: &[crate::align::data_structures::AnchorSeed]| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()) {
                            panic!("{}", a);
                        }

                        // eprintln!("Set score {}", a.score);
                        if a.flagged_for_indel {
                            eprintln!("Heyu");
                        }
                    }
                    // eprintln!("{}", query.len());
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        panic!("B 2  {}", a);
                    }
                    // println!("3 Extended {}", a);
                },
                None => (),
            };

            match a1 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("V 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("V 2  {}", a);
                    }
                }, _ => {},
            }
            // println!("END--/1 {:?}", a1);
            // println!("END--/2 {:?}\n-------------------------------", a2);
        });

        // println!("AGAIN");
        // anchors.iter().for_each(|AnchorPair(a1, a2)| {
        //     println!(" /1 {:?}", a1);
        //     println!(" /2 {:?}", a2);
        // });
        anchors.iter().for_each(|AnchorPair(a1, a2)| {
            match a1 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("T 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("T 2  {}", a);
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
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("U 1  {}", a);
                    }
                }, _ => {},
            }

            match a2 {
                Some(a) => {
                    if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                        eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
                        panic!("U 2  {}", a);
                    }
                }, _ => {},
            }
        });
    }
}

