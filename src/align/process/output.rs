use crate::{align::{common::{Or, PAFOutput, SAMOutput}, sam::Cigar}, database::common::FlexalignDatabase, io::output_buffer::OutputBuffer, options::Options};


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

impl SAMOutput for StdSAMOutput {
    fn write_header(&mut self, references: &[(&str, usize)]) {
        // SO:unsorted -- workers write records as reads stream past, in no reference order.
        let mut header = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for (name, length) in references {
            header.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
        }
        header.push_str(&format!("@PG\tID:flexalign\tPN:flexalign\tVN:{}\n", env!("CARGO_PKG_VERSION")));
        self.buffer.write(header);
    }

    fn write(
        &mut self,
        query_name: &str,
        flag: u16,
        reference_name: &str,
        reference_start: usize,
        mapping_quality: u8,
        cigar: &Cigar,
        mate_reference_name: Option<&str>,
        mate_reference_start: usize,
        template_length: i64,
        seq: &[u8],
        qual: &[u8],
    ) {
        // A CIGAR that does not fit the read would make an unreadable record; skip rather than
        // emit something samtools will reject.
        let cigar_string = match cigar.to_sam(seq.len()) {
            Some(c) => c,
            None => return,
        };

        // RNEXT is "=" when the mate is on this same reference, which is the common case.
        let rnext = match mate_reference_name {
            Some(name) if name == reference_name => "=",
            Some(name) => name,
            None => "*",
        };

        self.buffer.write(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n",
            query_name,
            flag,
            reference_name,
            reference_start + 1, // SAM POS is 1-based
            mapping_quality,
            cigar_string,
            rnext,
            if mate_reference_name.is_some() { mate_reference_start + 1 } else { 0 },
            template_length,
            String::from_utf8_lossy(seq),
            if qual.is_empty() { "*".to_string() } else { String::from_utf8_lossy(qual).to_string() },
            cigar.edit_distance(),
        ));
    }
}


/// Builds a run's output sink: SAM when `--sam` is given, PAF otherwise.
///
/// For SAM this also emits the header and flushes it, which must happen before any worker writes a
/// record -- the returned sink is cloned per worker, so the header would otherwise race the records
/// (or be duplicated into every clone).
pub fn make_output<D: FlexalignDatabase>(
    options: &Options,
    buffer: OutputBuffer,
    db: &D,
) -> Or<StdPAFOutput, StdSAMOutput> {
    if !options.args.sam {
        return Or { a: Some(StdPAFOutput::new(buffer)), b: None };
    }

    // The database exposes references by id without a count; ids are dense from 0, so walk until
    // one is missing. Zero-length entries are skipped: id 0 is flexmap's "dummy" sentinel (see
    // `flexmap::build`, which pushes it so real ids start at 1), and SAM requires LN >= 1 anyway.
    let mut references: Vec<(&str, usize)> = Vec::new();
    let mut id = 0usize;
    while let Some(name) = db.get_rname(id) {
        let length = db.get_reference(id).map_or(0, |r| r.len());
        if length > 0 {
            references.push((name, length));
        }
        id += 1;
    }

    let mut sam = StdSAMOutput::new(buffer);
    sam.write_header(&references);
    sam.buffer.flush();

    Or { a: None, b: Some(sam) }
}
