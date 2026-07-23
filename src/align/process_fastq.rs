use std::{fs::File, io::{self}, path::PathBuf, str::FromStr, sync::{Arc, Mutex}};

use bioreader::{parallel::fastq::{read_fastq_paired_end_state_par, read_fastq_single_end_state_par}, sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord}, utils::is_gzip};
use flate2::read::GzDecoder;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;
use log::info;

use crate::{
    GOLDSTD_EVAL, align::{
        common::{NoSAMOutput, Or},
        modular_workflow::{Modular, ModularPE}, 
        process::{
            alignment::LIBWFA2Alignment, anchor_extender::StdPairedAnchorExtender, anchor_extractor::{StdAnchorExtractor, StdPairedAnchorExtractor}, kmer_extractor::StdKmerExtractor, output::{StdPAFOutput, StdSAMOutput}, range_extractor::StdRangeExtractor, seed_extractor::StdSeedExtractor
        }, 
        stats::Stats, 
        workflow
    }, database::common::FlexalignDatabase, io::output_buffer::{OutputBuffer, OutputTarget}, options::Options};


pub fn process_fastq_wrapper<
        const K: usize, 
        const C: usize, 
        const F: usize, 
        const S: usize, 
        const L: usize,
        const HEADER_THRESHOLD: usize,
        FM: FlexalignDatabase + Clone + Sync + Send,
    >(options: &Options, db: &FM) {

    for (fwd, rev_option) in options.fwd.iter().zip(options.rev.iter()) {

        let file_fwd = match File::open(fwd) {
            Err(why) => panic!("couldn't open {}: {}", &fwd.to_str().unwrap(), why),
            Ok(file) => file,
        };

        
        let stdout_writer = Arc::new(Mutex::new(OutputTarget::Stdout(io::stdout())));
        let stdout_buffer_fwd = OutputBuffer::new(Arc::clone(&stdout_writer), 2usize.pow(24));

        let mut handler_fwd: workflow::Standard<K, C, F, S, L, HEADER_THRESHOLD, ClosedSyncmer<C, S, L>, FM> = 
            workflow::Standard::new(&db, ClosedSyncmer::<C,S,L>::new(), &options, stdout_buffer_fwd);

        let fwd_gzip = is_gzip(fwd).expect(format!("Cannot check if file is gzipped. Check file: {}", fwd.to_str().unwrap()).as_str());

        let stats;

        // Distinguish between single- and paired-end reads
        match rev_option {
            // Paired-end reads
            Some(rev) => { 
                log::info!("Paired-end processing");
                info!("Iterate {} {}", &fwd.to_str().unwrap(), &rev.to_str().unwrap());
                let file_rev = match File::open(rev) {
                    Err(why) => panic!("couldn't open {}: {}", &rev.to_str().unwrap(), why),
                    Ok(file) => file,
                };

                let rev_gzip = is_gzip(rev).expect(format!("Cannot check if file is gzipped. Check file: {}", rev.to_str().unwrap()).as_str());

                if fwd_gzip != rev_gzip {
                    panic!("Reads must either both be compressed (.gz) or uncompressed.")
                };

                let stdout_buffer_rev = OutputBuffer::new(Arc::clone(&stdout_writer), 2usize.pow(24));
                let mut handler_rev: workflow::Standard<K, C, F, S, L, HEADER_THRESHOLD, ClosedSyncmer<C, S, L>, FM> = 
                    workflow::Standard::new(&db, ClosedSyncmer::<C,S,L>::new(), &options, stdout_buffer_rev);


                let worker = move |rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats| {
                    handler_fwd.run(rec_fwd, stats);
                    handler_rev.run(rec_rev, stats);
                };
                
                let mut state = Stats::default();
                
                eprintln!("GOLD STANDARD {}", GOLDSTD_EVAL);
                if GOLDSTD_EVAL {
                    let mapqeval = state.gold_std_evaluation.as_mut().unwrap();

                    let fp_read1_outpath = PathBuf::from_str(&format!("{}.fp", fwd.to_str().unwrap())).unwrap();
                    let fp_read1_file = File::create_new(fp_read1_outpath.clone()).unwrap();
                    let fp_read1_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read1_file)));
                    
                    // let fp_read2_outpath = PathBuf::from_str(&format!("{}.fp", rev.to_str().unwrap())).unwrap();
                    // let fp_read2_file = File::create_new(fp_read2_outpath.clone()).unwrap();
                    // let fp_read2_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read2_file)));
                    
                    // eprintln!("Set false reads outpt path {:?} {:?}", fp_read1_outpath, fp_read2_outpath);
                    // mapqeval.set_fp_read_output(Some(fp_read1_writer), Some(fp_read2_writer));

                    // exit(9);
                    mapqeval.set_fp_read_output(Some(fp_read1_writer), None); 
                }

                if fwd_gzip {
                    stats = read_fastq_paired_end_state_par(
                        GzDecoder::new(file_fwd),
                        GzDecoder::new(file_rev),
                        usize::pow(2, 24),
                        options.args.threads,
                        state,
                        worker,
                    );
                } else {
                    stats = read_fastq_paired_end_state_par(
                        file_fwd,
                        file_rev,
                        usize::pow(2, 24),
                        options.args.threads,
                        state,
                        worker,
                    );
                }
            },
            // Single-end read
            None => {
                log::info!("Single-end processing");

                let worker = move |rec: &RefFastqRecord, stats: &mut Stats| {
                    handler_fwd.run(rec, stats);
                };

                if fwd_gzip {
                    stats = read_fastq_single_end_state_par(
                        GzDecoder::new(file_fwd),
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                } else {
                    stats = read_fastq_single_end_state_par(
                        file_fwd,
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                }
            },
        }

        eprintln!("{}", stats.as_ref().unwrap());
        if options.args.debug {
            let s = stats.as_ref().unwrap();
            let pc = |n: usize| if s.dbg_checked > 0 { 100.0 * n as f64 / s.dbg_checked as f64 } else { 0.0 };
            eprintln!(
                "[debug] reads checked (true ref known): {}\n\
                 [debug] true anchor present among all anchors: {} ({:.2}%)\n\
                 [debug] true anchor is best after extension:   {} ({:.2}%)",
                s.dbg_checked,
                s.dbg_true_anchor_present, pc(s.dbg_true_anchor_present),
                s.dbg_true_anchor_best, pc(s.dbg_true_anchor_best),
            );
        }
        // stats.as_ref().unwrap().plot_mapq();
        // dbg!(stats);
    };

}



pub fn process_fastq_wrapper_modular<
        'a,
        const K: usize, 
        const C: usize, 
        const F: usize, 
        const S: usize, 
        const L: usize,
        const HEADER_THRESHOLD: usize,
        FM: FlexalignDatabase + Clone + Sync + Send,
    >(options: &Options, db: &FM) {

    for (index, (fwd, rev_option)) in options.fwd.iter().zip(options.rev.iter()).enumerate() {

        let file_fwd = match File::open(fwd) {
            Err(why) => panic!("couldn't open {}: {}", &fwd.to_str().unwrap(), why),
            Ok(file) => file,
        };

        eprintln!("Process: {:?} {:?}", fwd, rev_option);


        let out_buffer = if options.output_prefix.is_some() {
            let path: &std::path::PathBuf = options.output_prefix.as_ref().unwrap().get(index).expect(&format!("There is no output for input {:?}", fwd));
            let file_writer = Arc::new(Mutex::new(OutputTarget::File(File::create(path).expect(&format!("Cannot open output file {:?}", path)))));
            OutputBuffer::new(Arc::clone(&file_writer), 2usize.pow(24))
        } else {
            let stdout_writer = Arc::new(Mutex::new(OutputTarget::Stdout(io::stdout())));
            OutputBuffer::new(Arc::clone(&stdout_writer), 2usize.pow(24))
        };

        // let output: StdPAFOutput = StdPAFOutput::new(stdout_buffer);
        let output: Or<StdPAFOutput, NoSAMOutput> = Or::<StdPAFOutput, NoSAMOutput> {
            a: Some(StdPAFOutput::new(out_buffer)),
            b: None,
        };


        let mut modular_fwd = Modular {
            options,
            db,
            kmer_extractor: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
            range_extractor: StdRangeExtractor::<K, C, F, FM>::new(db),
            seed_extractor: StdSeedExtractor::<K, C, F>::new(
                options.args.max_best_flex,
                options.args.max_range_size,
                options.args.min_ranges,
                options.args.ranges,
                options.args.mask_flank_mult,
            ),
            anchor_extractor: StdAnchorExtractor::new(),
            rec_rev: OwnedFastqRecord::new(),
            output: output.clone(),
        };        


        let fwd_gzip = is_gzip(fwd).expect(format!("Cannot check if file is gzipped. Check file: {}", fwd.to_str().unwrap()).as_str());

        let stats;


        // Distinguish between single- and paired-end reads
        match rev_option {
            // Paired-end reads
            Some(rev) => { 
                info!("Iterate {} {}", &fwd.to_str().unwrap(), &rev.to_str().unwrap());
                let file_rev = match File::open(rev) {
                    Err(why) => panic!("couldn't open {}: {}", &rev.to_str().unwrap(), why),
                    Ok(file) => file,
                };

                let rev_gzip = is_gzip(rev).expect(format!("Cannot check if file is gzipped. Check file: {}", rev.to_str().unwrap()).as_str());

                if fwd_gzip != rev_gzip {
                    panic!("Reads must either both be compressed (.gz) or uncompressed.")
                };


                let mut modular_rev = Modular {
                    options,
                    db,
                    kmer_extractor: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
                    range_extractor: StdRangeExtractor::<K, C, F, FM>::new(db),
                    seed_extractor: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult,
                    ),
                    anchor_extractor: StdAnchorExtractor::new(),
                    rec_rev: OwnedFastqRecord::new(),
                    // output_paf: Some(output),
                    // output_sam: None::<NoSAMOutput>,
                    output: output.clone(),
                };  


                let mut modular_pe = ModularPE {
                    options,
                    db,
                    kmer_extractor_fwd: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
                    kmer_extractor_rev: StdKmerExtractor::<K, C, ClosedSyncmer<C, S, L>>::default(),
                    range_extractor_fwd: StdRangeExtractor::<K, C, F, FM>::new(db),
                    range_extractor_rev: StdRangeExtractor::<K, C, F, FM>::new(db),
                    seed_extractor_fwd: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult,
                    ),
                    seed_extractor_rev: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult,
                    ),
                    anchor_extractor: StdPairedAnchorExtractor::new(),
                    anchor_extender: StdPairedAnchorExtender::new(db),
                    align: LIBWFA2Alignment::default(),
                    output: output,
                    rec_fwd_revc: OwnedFastqRecord::new(),
                    rec_rev_revc: OwnedFastqRecord::new(),
                };  


                let _worker = move |rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats| {
                    modular_fwd.run(rec_fwd, stats);
                    modular_rev.run(rec_rev, stats);
                };

                let worker_pe = move |rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats| {
                    
                    modular_pe.run(rec_fwd, rec_rev, stats);
                };

                
                let mut state = Stats::default();
                
                eprintln!("GOLD STANDARD {}", GOLDSTD_EVAL);
                if GOLDSTD_EVAL {
                    let mapqeval = state.gold_std_evaluation.as_mut().unwrap();

                    let fp_read1_outpath = PathBuf::from_str(&format!("{}.fp", fwd.to_str().unwrap())).unwrap();
                    let fp_read1_file = File::create(fp_read1_outpath.clone()).unwrap();
                    let fp_read1_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read1_file)));
                    
                    let fp_read2_outpath = PathBuf::from_str(&format!("{}.fp", rev.to_str().unwrap())).unwrap();
                    let fp_read2_file = File::create(fp_read2_outpath.clone()).unwrap();
                    let fp_read2_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read2_file)));
                    
                    // eprintln!("Set false reads outpt path {:?} {:?}", fp_read1_outpath, fp_read2_outpath);
                    // mapqeval.set_fp_read_output(None, Some(fp_read2_writer)); 
                    mapqeval.set_fp_read_output(Some(fp_read1_writer), Some(fp_read2_writer)); 
                }


                if fwd_gzip {
                    stats = read_fastq_paired_end_state_par(
                        GzDecoder::new(file_fwd),
                        GzDecoder::new(file_rev),
                        usize::pow(2, 24),
                        options.args.threads,
                        state,
                        worker_pe,//worker,
                    );
                } else {
                    stats = read_fastq_paired_end_state_par(
                        file_fwd,
                        file_rev,
                        usize::pow(2, 24),
                        options.args.threads,
                        state,
                        worker_pe,//worker,
                    );
                }
            },
            // Single-end read
            None => {

                let worker = move |rec: &RefFastqRecord, stats: &mut Stats| {
                    modular_fwd.run(rec, stats);
                };

                if fwd_gzip {
                    stats = read_fastq_single_end_state_par(
                        GzDecoder::new(file_fwd),
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                } else {
                    stats = read_fastq_single_end_state_par(
                        file_fwd,
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                }
            },
        }

        eprintln!("{}", stats.as_ref().unwrap());
        if options.args.debug {
            let s = stats.as_ref().unwrap();
            let pc = |n: usize| if s.dbg_checked > 0 { 100.0 * n as f64 / s.dbg_checked as f64 } else { 0.0 };
            eprintln!(
                "[debug] reads checked (true ref known): {}\n\
                 [debug] true anchor present among all anchors: {} ({:.2}%)\n\
                 [debug] true anchor is best after extension:   {} ({:.2}%)",
                s.dbg_checked,
                s.dbg_true_anchor_present, pc(s.dbg_true_anchor_present),
                s.dbg_true_anchor_best, pc(s.dbg_true_anchor_best),
            );
        }
        // stats.as_ref().unwrap().plot_mapq();
        // dbg!(stats);
    };

}

