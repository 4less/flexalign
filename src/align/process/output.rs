use crate::{align::common::PAFOutput, io::output_buffer::OutputBuffer};


#[derive(Clone)]
pub struct StdPAFOutput {
    pub buffer: OutputBuffer,
}

impl StdPAFOutput {
    pub fn new(buffer: OutputBuffer) -> Self {
        Self {
            buffer
        }
    }
}

impl PAFOutput for StdPAFOutput {
    fn write(
        &mut self,
        query_name: &str,
        query_length: usize,
        query_start: i32,
        query_end: i32,
        fwd: bool,
        reference_name: &str,
        reference_length: usize,
        reference_start: i32,
        reference_end: i32,
        residue_matches: u32,
        alignment_block_length: usize,
        mapping_quality: u8,
    ) {
        self.buffer.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            query_name, 
            query_length,
            query_start,
            query_end,
            if fwd { '+' } else { '-' },
            reference_name,
            reference_length,
            reference_start,
            reference_end,
            residue_matches,
            alignment_block_length,
            mapping_quality));
    }
}


#[derive(Clone)]
pub struct StdSAMOutput {
    pub buffer: OutputBuffer,
}

impl StdSAMOutput {
    pub fn new(buffer: OutputBuffer) -> Self {
        Self {
            buffer
        }
    }
}

