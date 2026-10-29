use std::{fs::File, io::{self}, path::PathBuf, str::FromStr, sync::{Arc, Mutex}};

use bioreader::{fastq_byte_reader::FastqPairedByteReader, fastq_reader::PairedFastqReader, parallel::fastq::{read_fastq_paired_end_state_par, read_fastq_single_end_state_par}, sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord}, utils::is_gzip};
use flate2::read::{GzDecoder, MultiGzDecoder};
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;
use log::info;

use crate::{
    GOLDSTD_EVAL, align::{
        modular_workflow::{self, FastHitScratch, Modular, ModularPE}, 
        process::{
            anchor_extractor::StdAnchorExtractor, output::make_output, range_extractor::StdRangeExtractor, seed_extractor::StdSeedExtractor, syncmer_extractor::SyncmerKmerExtractor
        }, 
        stats::Stats, 
        workflow
    }, database::common::FlexalignDatabase, io::{output_buffer::{OutputBuffer, OutputTarget}, threaded_reader::ThreadedReader}, options::Options};


/// Where the false-positive reads of `reads_file` are written (gold-standard evaluation only).
///
/// `--output-prefix-fp-reads PREFIX` puts them at `PREFIX<read file name>.fp`; without it they land
/// beside the reads as `<reads_file>.fp`, which is what the code did unconditionally before -- the
/// flag existed but was never read, so passing it silently changed nothing.
fn fp_reads_path(options: &Options, reads_file: &std::path::Path) -> String {
    match &options.args.output_prefix_fp_reads {
        Some(prefix) => {
            let name = reads_file
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "reads".to_string());
            format!("{prefix}{name}.fp")
        }
        None => format!("{}.fp", reads_file.to_string_lossy()),
    }
}


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
                
                if GOLDSTD_EVAL {
                    let mapqeval = state.gold_std_evaluation.as_mut().unwrap();

                    // --output-prefix-fp-reads, when given, decides where these go; otherwise they
                    // land beside the reads as <reads>.fp (the previous hardcoded behaviour). The
                    // flag was declared and never read, so setting it silently did nothing.
                    let fp_read1_outpath = PathBuf::from_str(&fp_reads_path(options, fwd)).unwrap();
                    // `create_new` fails if the file exists, so a SECOND run panicked. Truncate
                    // instead, matching the equivalent path in the modular wrapper below.
                    let fp_read1_file = File::create(fp_read1_outpath.clone()).unwrap();
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
        if let Some(path) = &options.args.time_log {
            let input = fwd.to_str().unwrap_or("input");
            if let Err(e) = stats.as_ref().unwrap().append_time_log(path, input) {
                eprintln!("warning: could not write --time-log {}: {}", path, e);
            }
        }
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



/// Bytes of FASTQ each worker takes from the shared reader per batch (per mate file).
///
/// 16 MiB batches left 16 threads ~26% idle on a 650 MB sample: 41 batches over 16 workers is
/// 2.5 rounds, so the last round ran on 9 threads while 7 waited -- a tail that costs a full
/// batch of wall time. Smaller batches balance the tail at the price of more lock acquisitions,
/// which are cheap next to the work a batch carries: 4 MiB took the wall from 2.47 s to 2.09 s,
/// and once the fast-hit path made each pair cheaper, 1 MiB measured 1.13 s against 1.20 s.
const READ_BATCH_BYTES: usize = 1 << 20;

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
        // Phases of one input, outside the per-pair stages: the SAM header is written on one
        // thread before any read is processed, and on a large database it rivals the whole pass.
        let input_started = std::time::Instant::now();


        let out_buffer = if options.output_prefix.is_some() {
            let path: &std::path::PathBuf = options.output_prefix.as_ref().unwrap().get(index).expect(&format!("There is no output for input {:?}", fwd));
            let file_writer = Arc::new(Mutex::new(OutputTarget::File(File::create(path).expect(&format!("Cannot open output file {:?}", path)))));
            OutputBuffer::new(Arc::clone(&file_writer), 2usize.pow(24))
        } else {
            let stdout_writer = Arc::new(Mutex::new(OutputTarget::Stdout(io::stdout())));
            OutputBuffer::new(Arc::clone(&stdout_writer), 2usize.pow(24))
        };

        // let output: StdPAFOutput = StdPAFOutput::new(stdout_buffer);
        let output = make_output(options, out_buffer, db);
        let header_written = input_started.elapsed();
        // Reset when the parallel pass starts, after the worker setup and the fast-hit probe.
        let mut pass_started = input_started;


        let mut modular_fwd = Modular {
            options,
            db,
            kmer_extractor: SyncmerKmerExtractor::<K, C, S, L>::default(),
            range_extractor: StdRangeExtractor::<K, C, F, FM>::new(db),
            seed_extractor: StdSeedExtractor::<K, C, F>::new(
                options.args.max_best_flex,
                options.args.max_range_size,
                options.args.min_ranges,
                options.args.ranges,
                options.args.mask_flank_mult, options.args.flank_slack,
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
                    kmer_extractor: SyncmerKmerExtractor::<K, C, S, L>::default(),
                    range_extractor: StdRangeExtractor::<K, C, F, FM>::new(db),
                    seed_extractor: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult, options.args.flank_slack,
                    ),
                    anchor_extractor: StdAnchorExtractor::new(),
                    rec_rev: OwnedFastqRecord::new(),
                    // output_paf: Some(output),
                    // output_sam: None::<NoSAMOutput>,
                    output: output.clone(),
                };  


                let (anchor_extractor, anchor_extender, align) = modular_workflow::paired_align_half::<K, _>(options, db);
                let mut modular_pe = ModularPE {
                    options,
                    db,
                    kmer_extractor_fwd: SyncmerKmerExtractor::<K, C, S, L>::default(),
                    kmer_extractor_rev: SyncmerKmerExtractor::<K, C, S, L>::default(),
                    range_extractor_fwd: StdRangeExtractor::<K, C, F, FM>::new(db),
                    range_extractor_rev: StdRangeExtractor::<K, C, F, FM>::new(db),
                    seed_extractor_fwd: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult, options.args.flank_slack,
                    ),
                    seed_extractor_rev: StdSeedExtractor::<K, C, F>::new(
                        options.args.max_best_flex,
                        options.args.max_range_size,
                        options.args.min_ranges,
                        options.args.ranges,
                        options.args.mask_flank_mult, options.args.flank_slack,
                    ),
                    anchor_extractor,
                    anchor_extender,
                    align,
                    output: output,
                    rec_fwd_revc: OwnedFastqRecord::new(),
                    rec_rev_revc: OwnedFastqRecord::new(),
                    fast: FastHitScratch::default(),
                };  


                let _worker = move |rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats| {
                    modular_fwd.run(rec_fwd, stats);
                    modular_rev.run(rec_rev, stats);
                };

                // The fast-hit probe: try the path on the input's first pairs, once, before the
                // workers start, and switch it off for the run when too few succeed. One decision
                // for the whole run, from the same pairs every time, so the output is the same
                // run to run (a per-worker probe depended on thread scheduling).
                if options.args.fast_hit_on() {
                    let probe_started = std::time::Instant::now();
                    let mut probe_stats = Stats::default();
                    let (mut tried, mut hit) = (0usize, 0usize);
                    {
                        let f1 = File::open(fwd).expect("reads");
                        let f2 = File::open(rev).expect("reads");
                        let (b1, b2): (Box<dyn std::io::Read + Send>, Box<dyn std::io::Read + Send>) = if fwd_gzip {
                            (Box::new(GzDecoder::new(f1)), Box::new(GzDecoder::new(f2)))
                        } else {
                            (Box::new(f1), Box::new(f2))
                        };
                        let byte_reader = Arc::new(Mutex::new(FastqPairedByteReader::new(b1, b2, READ_BATCH_BYTES)));
                        let mut reader = PairedFastqReader::new(byte_reader, READ_BATCH_BYTES);
                        while let Some((r1, r2)) = reader.next() {
                            tried += 1;
                            if modular_pe.probe_fast_hit(&r1, &r2, &mut probe_stats) {
                                hit += 1;
                            }
                            if tried >= modular_workflow::FAST_HIT_PROBE {
                                break;
                            }
                        }
                    }
                    let off = tried > 0 && hit * 5 < tried * modular_workflow::FAST_HIT_MIN_FIFTHS;
                    modular_workflow::set_fast_hit_enabled(!off);
                    eprintln!("Fast-hit probe: {} of the first {} pairs placed -> {} ({:.1?})", hit, tried, if off { "off for this run" } else { "on" }, probe_started.elapsed());
                }

                let worker_pe = move |rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats| {
                    
                    modular_pe.run(rec_fwd, rec_rev, stats);
                };

                
                let mut state = Stats::default();
                
                if GOLDSTD_EVAL {
                    let mapqeval = state.gold_std_evaluation.as_mut().unwrap();

                    let fp_read1_outpath = PathBuf::from_str(&fp_reads_path(options, fwd)).unwrap();
                    let fp_read1_file = File::create(fp_read1_outpath.clone()).unwrap();
                    let fp_read1_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read1_file)));
                    
                    let fp_read2_outpath = PathBuf::from_str(&fp_reads_path(options, rev)).unwrap();
                    let fp_read2_file = File::create(fp_read2_outpath.clone()).unwrap();
                    let fp_read2_writer = Arc::new(Mutex::new(OutputTarget::File(fp_read2_file)));
                    
                    // eprintln!("Set false reads outpt path {:?} {:?}", fp_read1_outpath, fp_read2_outpath);
                    // mapqeval.set_fp_read_output(None, Some(fp_read2_writer)); 
                    mapqeval.set_fp_read_output(Some(fp_read1_writer), Some(fp_read2_writer)); 
                }


                pass_started = std::time::Instant::now();
                if fwd_gzip {
                    // Each mate file is inflated on a thread of its own, ahead of the workers
                    // (io::threaded_reader): under the shared reader lock the inflate starved them.
                    stats = read_fastq_paired_end_state_par(
                        ThreadedReader::new(MultiGzDecoder::new(file_fwd)),
                        ThreadedReader::new(MultiGzDecoder::new(file_rev)),
                        READ_BATCH_BYTES,
                        options.args.threads,
                        state,
                        worker_pe,//worker,
                    );
                } else {
                    stats = read_fastq_paired_end_state_par(
                        file_fwd,
                        file_rev,
                        READ_BATCH_BYTES,
                        options.args.threads,
                        state,
                        worker_pe,//worker,
                    );
                }
            },
            // Single-end read
            None => {
                // The single-end path runs `Modular`, which stops at anchors and never aligns, so
                // it has no CIGAR to put in a SAM record. SAM is the default now, so this is the
                // common case rather than a user error: say what is being written and carry on.
                // `Args::emit_sam` already resolved the format to PAF here.
                if !options.args.paf {
                    log::info!("single-end input: writing PAF (SAM needs the gapped alignment that \
                                only the paired-end pipeline runs)");
                }

                let worker = move |rec: &RefFastqRecord, stats: &mut Stats| {
                    modular_fwd.run(rec, stats);
                };

                if fwd_gzip {
                    stats = read_fastq_single_end_state_par(
                        GzDecoder::new(file_fwd),
                        READ_BATCH_BYTES,
                        options.args.threads,
                        worker,
                    );
                } else {
                    stats = read_fastq_single_end_state_par(
                        file_fwd,
                        READ_BATCH_BYTES,
                        options.args.threads,
                        worker,
                    );
                }
            },
        }

        if crate::flexalign::STAGE_TIMING.load(std::sync::atomic::Ordering::Relaxed) {
            eprintln!(
                "Input phases: SAM header {:.3?}, setup {:.3?}, read pass {:.3?}",
                header_written,
                pass_started.duration_since(input_started).saturating_sub(header_written),
                pass_started.elapsed(),
            );
        }
        let mut stats = stats;
        if let Ok(s) = stats.as_mut() {
            s.fast_hit_workers_off = if options.args.fast_hit_on() && !modular_workflow::fast_hit_enabled() { 1 } else { 0 };
        }
        eprintln!("{}", stats.as_ref().unwrap());
        if let Some(path) = &options.args.time_log {
            let input = fwd.to_str().unwrap_or("input");
            if let Err(e) = stats.as_ref().unwrap().append_time_log(path, input) {
                eprintln!("warning: could not write --time-log {}: {}", path, e);
            }
        }
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

