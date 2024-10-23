use std::cmp::min;

use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{align::eval::{BinaryEvaluator, MapqEvaluation}, database::common::FlexalignDatabase};

pub fn get_id_from_header(header_str: &str, db: &impl FlexalignDatabase) -> usize {
    let first_part_a = header_str.split('-').next().unwrap_or("");
    let first_part_b = header_str.splitn(3, '_').take(2).collect::<Vec<&str>>().join("_");
    let mut true_id = *db.get_rid(first_part_a).unwrap_or(&0);
    if true_id == 0 {
        true_id = *db.get_rid(&first_part_b).unwrap_or(&0);
    }
    true_id
}

pub fn correct(header_str: &[u8], reference: u64, db: &impl FlexalignDatabase) -> bool {
    let ref_string = db.get_rname(reference as usize).unwrap();
    let correct = &ref_string.as_bytes()[..min(ref_string.len(), header_str.len())] == &header_str[..min(ref_string.len(), header_str.len())];
    correct
}

pub fn evaluate(eval: &mut MapqEvaluation, refstr: &str, pseudo_mapq: u64, rec: &RefFastqRecord, _db: &impl FlexalignDatabase) {
    // let header_str = String::from_utf8_lossy(rec.head());
    // let first_part_a = header_str.split('-').next().unwrap_or("");
    // let first_part_b = header_str.splitn(3, '_').take(2).collect::<Vec<&str>>().join("_");
    // let mut true_id = *db.get_rid(first_part_a).unwrap_or(&0);

    // if true_id == 0 {
    //     true_id = *db.get_rid(&first_part_b).unwrap_or(&0);
    // }

    // if true_id == 0 {
    //     panic!("True id is {}", true_id);
    // }


    let correct = &refstr.as_bytes()[..min(refstr.len(), rec.head().len())] == &rec.head()[..min(refstr.len(), rec.head().len())];
    // eprintln!("{}\t{}\t{}\t{}", ref_string, header_str, correct, pseudo_mapq);

    eval.add(correct, pseudo_mapq);
}