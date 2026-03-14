use std::{cmp::{min, max}, ops::Range};

use derive_builder::Builder;

fn calculate_ranges(
    seq_len: usize,
    seed_size: usize,
    mut seed_num: usize,
    min_gap_size: usize,
) -> Option<Vec<Range<usize>>> {
    if seed_size > seq_len {
        return None;
    }

    while ((seed_num - 1) * (seed_size + min_gap_size)) + seed_size > seq_len {
        seed_num -= 1;
    }

    let total_gap = seq_len - (seed_num * seed_size); // Total gap to distribute
    let gap = if seed_num <= 1 { 0 } else {
        total_gap / (seed_num - 1)
    };  // Gap size between ranges

    let mut ranges = Vec::new();
    for i in 0..seed_num {
        let start = i * (gap + seed_size); // Position for the start of the range
        ranges.push(start..start + seed_size); // Create a range with size s
    }

    Some(ranges)
}

#[derive(Builder, Debug)]
#[builder(default)]
pub struct HasIndelsParams {
    seed_size: usize,
    seed_num: usize,
    min_gap_size: usize,
    max_mismatches: u32,
    max_indels: usize,
}

impl Default for HasIndelsParams {
    fn default() -> Self {
        Self {
            seed_size: 8,
            seed_num: 3,
            min_gap_size: 2,
            max_mismatches: 1,
            max_indels: 5,
        }
    }
}

fn has_indels(
    query: &[u8],
    reference: &[u8],
    ref_alignment: usize,
    params: &HasIndelsParams,
) -> Option<i32> {
    let sims = offset_similarities(
        query,
        reference,
        ref_alignment,
        params.max_indels,
        params.seed_num,
        params.seed_size,
        params.min_gap_size,
    )?;
    let sim = sims.iter()
        .filter(|(offset, hamming)| *offset != 0 && hamming <= &params.max_mismatches).min_by_key(|x| x.1)?;

    Some(sim.0)
}

fn offset_similarities(
    query: &[u8],
    reference: &[u8],
    ref_alignment: usize,
    max_indel: usize,
    seed_num: usize,
    seed_size: usize,
    min_gap_size: usize,
) -> Option<Vec<(i32, u32)>> {
    assert!(ref_alignment + seed_size <= reference.len());

    let len = min(query.len(), reference.len() - ref_alignment);
    let qpositions = calculate_ranges(len, seed_size, seed_num, min_gap_size)?;

    let mut offsets = Vec::new();
    for q in qpositions.clone() {
        let rrange = (ref_alignment + q.start)..(ref_alignment + q.start + seed_size);
        let zero_offset_hamm = triple_accel::hamming(&query[q.clone()], &reference[rrange]);
        if zero_offset_hamm == 0 {
            offsets.push((0, 0));
            continue;
        }
        let res = best_offset_similarity(
            &query[q.clone()],
            reference,
            ref_alignment + q.start,
            max_indel,
        );
        if let Some(res) = res {
            offsets.push(res);
        }
    }

    Some(offsets)
}

fn best_offset_similarity(
    query: &[u8],
    reference: &[u8],
    ref_alignment: usize,
    max_indel: usize,
) -> Option<(i32, u32)> {
    let ql = query.len();
    let rl = reference.len();
    let offset_start = -(min(max_indel, ref_alignment) as i32);
    let offset_end = min(max_indel as i32, (rl - ref_alignment - ql) as i32);

    let offsets = offset_start..offset_end;
    // offsets.iter().map
    let best_sim = offsets
        .map(|offset| {
            let begin = (ref_alignment as i32 + offset) as usize;
            let ref_window = &reference[begin..begin + ql];
            (offset, triple_accel::hamming(query, ref_window))
        })
        .max_by_key(|p| -(p.1 as i64));

    best_sim
}


// Unit tests
#[cfg(test)] // Ensures the code is only compiled and run during testing
mod tests {
    
    use std::hint::black_box;

    use test::Bencher;

    use super::*; // Import the functions from the parent module

    #[test]
    fn test_best_offset_similarity() {
        //    ACGTCGTG
        // CCGACACGTCATGGCTAGCT
        let query = "ACGTCGTG";
        let reference = "CCGACACGTCATGGCTAGCT";
        let reference_offset = 3;

        assert_eq!(
            best_offset_similarity(query.as_bytes(), reference.as_bytes(), reference_offset, 5),
            Some((2, 1))
        )
    }

    #[test]
    fn test_offset_similarities() {
        //    ACGTCGTG
        // CCGACACGTCATGGCTAGCT
        let query =
            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference = "CCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTACGTGACTGCTGCGTCGTGAAAGCTCGTGTATCG";
        let reference_offset = 0;

        assert_eq!(
            offset_similarities(
                query.as_bytes(),
                reference.as_bytes(),
                reference_offset,
                5,
                3,
                8,
                5
            ),
            Some(vec![(0, 1), (0, 0), (-2, 1)])
        );

        //CCGACTCG               GATCTATC               CTATCTAC
        //CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACT
        //CCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC

        //    ACGTCGTG
        // CCGACACGTCATGGCTAGCT
        let query =
            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC";
        let reference_offset = 7;

        assert_eq!(
            offset_similarities(
                query.as_bytes(),
                reference.as_bytes(),
                reference_offset,
                5,
                3,
                8,
                5
            ),
            Some(vec![(0, 1), (0, 0), (-2, 0)])
        );

        //CCGACTCGTCATGGCT   CTACGATCTATCGTAG   ATCGTAACCTATCTAC
        //CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACT
        //CCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC
        assert_eq!(
            offset_similarities(
                query.as_bytes(),
                reference.as_bytes(),
                reference_offset,
                5,
                3,
                16,
                0
            ),
            Some(vec![(0, 1), (0, 0), (-2, 7)])
        );
    }



    #[test]
    fn test_has_indels() {
        let query =            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC";
        let reference_offset = 7;

        let has_indel_params = HasIndelsParams::default();

        assert_eq!(has_indels(query.as_bytes(), reference.as_bytes(), reference_offset, &has_indel_params), Some(-2));

        let query =            "CCGACTCGTCATGGCTCTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset = 7;

        let has_indel_params = HasIndelsParams::default();

        assert_eq!(has_indels(query.as_bytes(), reference.as_bytes(), reference_offset, &has_indel_params), None);


        let query =            "CCGACTCGTCATGTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset = 7;

        let has_indel_params = HasIndelsParams::default();

        assert_eq!(has_indels(query.as_bytes(), reference.as_bytes(), reference_offset, &has_indel_params), Some(4));
    }



    #[bench]
    fn bench_has_indels1(b: &mut Bencher) {
        let query1 =            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference1 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC";
        let reference_offset1 = 7;

        let has_indel_params = HasIndelsParamsBuilder::default()
            .seed_num(1usize).seed_size(4).build().unwrap();

        let query2 =            "CCGACTCGTCATGGCTCTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference2 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset2 = 7;


        let query3 =            "CCGACTCGTCATGTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference3 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset3 = 7;


        b.iter(|| {
            black_box(has_indels(query1.as_bytes(), reference1.as_bytes(), reference_offset1, &has_indel_params));
            black_box(has_indels(query2.as_bytes(), reference2.as_bytes(), reference_offset2, &has_indel_params));
            black_box(has_indels(query3.as_bytes(), reference3.as_bytes(), reference_offset3, &has_indel_params));
        })
    }


    #[bench]
    fn bench_has_indels2(b: &mut Bencher) {

        let query1 =            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference1 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC";
        let reference_offset1 = 7;

        // let has_indel_params = HasIndelsParamsBuilder::default()
        //     .seed_size(16usize)
        //     .build().unwrap();

        let has_indel_params = HasIndelsParamsBuilder::default()
            .seed_num(1usize).seed_size(16).build().unwrap();

        let query2 =            "CCGACTCGTCATGGCTCTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference2 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset2 = 7;


        let query3 =            "CCGACTCGTCATGTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference3 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset3 = 7;


        b.iter(|| {
            black_box(has_indels(query1.as_bytes(), reference1.as_bytes(), reference_offset1, &has_indel_params));
            black_box(has_indels(query2.as_bytes(), reference2.as_bytes(), reference_offset2, &has_indel_params));
            black_box(has_indels(query3.as_bytes(), reference3.as_bytes(), reference_offset3, &has_indel_params));
        })
    }



    #[bench]
    fn bench_has_indels3(b: &mut Bencher) {

        let query1 =            "CCGACTCGTCATGGCTCTACTACGATCTATCGTAGCTGATCGTAACCTATCTACTACGTGACTGCTTCGTCGTAAAAGCT";
        let reference1 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATCGTAGCTGATCGTACTATCTACTAC";
        let reference_offset1 = 7;

        let query2 =            "CCGACTCGTCATGGCTCTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference2 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset2: usize = 7;

        let query3 =            "CCGACTCGTCATGTACTACGATCTATCTTCGTCGTAAAAGCT";
        let reference3 = "CGATGCTCCGACACGTCATGGCTGTACTACGATCTATC";
        let reference_offset3 = 7;

        

        
        // let has_indel_params = HasIndelsParamsBuilder::default()
        //     .seed_size(16usize)
        //     .build().unwrap();

        let has_indel_params = HasIndelsParamsBuilder::default()
            .seed_num(1usize).seed_size(8).max_indels(10).build().unwrap();



        b.iter(|| {
            black_box(has_indels(query1.as_bytes(), reference1.as_bytes(), reference_offset1, &has_indel_params));
            black_box(has_indels(query2.as_bytes(), reference2.as_bytes(), reference_offset2, &has_indel_params));
            black_box(has_indels(query3.as_bytes(), reference3.as_bytes(), reference_offset3, &has_indel_params));
        })
    }
}

