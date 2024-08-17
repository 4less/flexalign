use std::{fs::File, io::{self}, sync::{Arc, Mutex}};

use bioreader::{parallel::fastq::{read_fastq_paired_end_state_par, read_fastq_single_end_state_par}, sequence::fastq_record::RefFastqRecord, utils::is_gzip};
use flate2::read::GzDecoder;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;
use log::info;

use crate::{align::{workflow, stats::Stats}, database::common::FlexalignDatabase, io::output_buffer::{OutputBuffer, OutputTarget}, options::Options};

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
                
                if fwd_gzip {
                    stats = read_fastq_paired_end_state_par(
                        GzDecoder::new(file_fwd),
                        GzDecoder::new(file_rev),
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                } else {
                    stats = read_fastq_paired_end_state_par(
                        file_fwd,
                        file_rev,
                        usize::pow(2, 24),
                        options.args.threads,
                        worker,
                    );
                }
            },
            // Single-end read
            None => {

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

        eprintln!("{}", stats.unwrap());
        // dbg!(stats);
    };

}

