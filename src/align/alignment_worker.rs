#[derive(Clone)]
pub struct AlignmentWorker<'a, const C: usize, const S: usize, const L: usize>{
    pub seedlist: Vec<VRange<'a>>,
    pub stats: Stats,
    pub cs: ClosedSyncmer::<C,S,L>,
    
}



impl<'a, const C: usize, const S: usize, const L: usize> FnMut<(&'_ RefFastqRecord<'_>, &'_ RefFastqRecord<'_>)> for AlignmentWorker<'a, C, S, L> {
    extern "rust-call" fn call_mut(&mut self, args: (&'_ RefFastqRecord<'_>, &'_ RefFastqRecord<'_>)) -> Self::Output {
        // Your implementation here
        
        let mut kmer_count = 0;

        let (rec1, rec2) = args;

        let mut iter1 = KmerIter::<K, true>::new(rec1.seq());
        let mut iter2 = KmerIter::<K, true>::new(rec2.seq());

        self.seedlist.clear();
        
        // iter1.set_assume_perfect_data(rec1.perfect());
        // iter2.set_assume_perfect_data(rec2.perfect());

        iter1.set_assume_perfect_data(false);
        iter2.set_assume_perfect_data(false);


        for (_pos, kmer_fwd, kmer_rev) in iter1 {
            self.stats.kmers_processed_fwd += 1;

            let kmer = min(kmer_fwd, kmer_rev);
            let cmer = kmer.middle::<C>();
            if !self.cs.is_minimizer(cmer.0) { continue };

            self.stats.minimizer_fwd += 1;

            let fmer = kmer.flanks::<F>();


            kmer_count += 1;

            // println!("{} {} {}", kmer.to_string().expect("Correct k-mer"), cmer.to_string().expect("Correct k-mer"), kmer.0);

            let range = match db.flexmap.get(cmer.0) {
                Some(range) => range,
                None => continue,
            };
            // println!("Range {}", range.positions.len());
            
            // println!("Range___________________________________\n{}", range.to_verbose_string::<20,40>());

            // let mut str = String::new();


            // match &range.header {
            //     Some(header) => {
            //         assert_eq!(header.len(), range.positions.len());
            //         for idx in 0..header.len() {
            //             let (val, pos) = VData::<20,40>::get(range.positions[idx].0);
            //             kmer_count += val;

            //             // let reference = &db.id2reference[val as usize];
            //             // str.push_str(&format!("{}: {} {}\n", header[idx].to_string(), reference, pos));
            //         }
            //     },
            //     None => {
            //         for idx in 0..range.positions.len() {
            //             let (val, pos) = VData::<20,40>::get(range.positions[idx].0);
            //             kmer_count += val;
            //             // let reference: &String = &db.id2reference[val as usize];
            //             // str.push_str(&format!("................: {} {}\n", reference, pos));
            //         }
            //     }
            // }
        
            kmer_count += range.positions.len() as u64;
            seedlist.push(range);
            kmer_count += seedlist.first().map_or(0, |r| r.positions.len()) as u64;
        }
        

        let (duration, _) = time(|| seedlist.sort_unstable());
        
        stats.seeds_fwd += seedlist.len();
        stats.time_seed_sorting_fwd += duration;


        if seedlist.len() < 3 {
            println!("Seeds {}", seedlist.len());
            println!("{}", rec1.to_string());
        }

        seedlist.clear();

        for (_pos, kmer_fwd, kmer_rev) in iter2 {
            let kmer = min(kmer_fwd, kmer_rev);
            let cmer = kmer.middle::<C>();
            if !cs.is_minimizer(cmer.0) { continue };

            let fmer = kmer.flanks::<F>();

            kmer_count += 1;

            let range = match db.flexmap.get(cmer.0) {
                Some(range) => range,
                None => continue,
            };

                            // println!("Range {}", range.positions.len());
                
                // println!("Range___________________________________\n{}", range.to_verbose_string::<20,40>());

                // let mut str = String::new();


                // match &range.header {
                //     Some(header) => {
                //         assert_eq!(header.len(), range.positions.len());
                //         for idx in 0..header.len() {
                //             let (val, pos) = VData::<20,40>::get(range.positions[idx].0);
                //             kmer_count += val;

                //             // let reference = &db.id2reference[val as usize];
                //             // str.push_str(&format!("{}: {} {}\n", header[idx].to_string(), reference, pos));
                //         }
                //     },
                //     None => {
                //         for idx in 0..range.positions.len() {
                //             let (val, pos) = VData::<20,40>::get(range.positions[idx].0);
                //             kmer_count += val;
                //             // let reference: &String = &db.id2reference[val as usize];
                //             // str.push_str(&format!("................: {} {}\n", reference, pos));
                //         }
                //     }
                // }

            kmer_count += range.positions.len() as u64;
            seedlist.push(range);
            seedlist.sort_unstable();
            kmer_count += seedlist.first().map_or(0, |r| r.positions.len()) as u64;
        }

        ()
    }
}

impl<'a, const C: usize, const S: usize, const L: usize> FnOnce<(&'_ RefFastqRecord<'_>, &'_ RefFastqRecord<'_>)> for AlignmentWorker<'a, C, S, L> {
    type Output = (); // Set this to the correct output type
    extern "rust-call" fn call_once(mut self, args: (&'_ RefFastqRecord<'_>, &'_ RefFastqRecord<'_>)) -> Self::Output {
        self.call_mut(args)
    }
}
