
use std::cmp::min;

use bioreader::sequence::fastq_record::RefFastqRecord;
use flexmap::values::VRange;
use kmerrs::consecutive::kmer::Kmer;
use kmerrs::kmer_utils::mask_coremer;

use crate::{align::{common::RangeExtractor, stats::Stats}, database::common::FlexalignDatabase};

// pub type Range<'a, const F: usize> = (usize, Kmer<F>, VRange<'a>, usize);

#[derive(Clone)]
pub struct Range<'a, const F: usize> {
    pub qpos: usize,
    pub fmer: Kmer<F>,
    pub vrange: VRange<'a>,
    pub range_size: u32,
    pub core_min_phredq: u16,
    pub flex_min_phredq: u16,
}

#[derive(Clone)]
pub struct StdRangeExtractor<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> {
    pub ranges: Vec<Range<'a, F>>,
    pub db: &'a D,
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> RangeExtractor<K, F> for StdRangeExtractor<'a, K, C, F, D> {
    fn retrieve(&self) -> &[Range<'_, F>] {
        &self.ranges
    }
    
    fn generate(&mut self, kmers: &[(usize, Kmer<K>)], rec: &RefFastqRecord, stats: &mut Stats) -> &[Range<'_, F>] {
        self.ranges.clear();

        for (pos, kmer) in kmers {
            let cmer = kmer.middle::<C>();
            let fmer = kmer.flanks::<F>();

            // Look up under the same de-biasing mix used at index build and syncmer selection.
            let result = self.db.get_vrange(mask_coremer::<C>(cmer.0));

            // eprintln!("{} {} -> {}", pos, kmer.to_string().unwrap(), result.is_some());
            let range: VRange<'a> = match result {
                Some(range) => range,
                None => { continue },
            };
            let range_len = (&range).positions.len();
            self.ranges.push(Range { qpos: *pos, fmer: fmer, vrange: range, range_size: range_len as u32, core_min_phredq: 0, flex_min_phredq: 0 });
        }
        
        for range in &mut self.ranges {
            let qual_mer = &rec.qual()[range.qpos..(range.qpos + C)];
            range.core_min_phredq = min(qual_mer.iter().map(|q| q - 33).min().unwrap() as u16, 40);
            // eprintln!("{} {} ... sort by {}", range.core_min_phredq, range.vrange.positions.len(), range.vrange.positions.len() * (41 - range.core_min_phredq as usize));
        }

        // self.ranges.sort_unstable_by_key(|r| r.vrange.positions.len() * (41 - r.core_min_phredq as usize));
        self.ranges.sort_unstable_by_key(|r| r.vrange.positions.len());
        
        // let (duration, _) = time(|| self.ranges.sort_unstable_by_key(|r| r.2.positions.len()));
        // stats.time_range_sorting += duration;
        stats.ranges += self.ranges.len();

        &self.ranges
    }
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> StdRangeExtractor<'a, K, C, F, D> {
    pub fn new(db: &'a D) -> Self {
        Self { ranges: Vec::new(), db }
    }
}
