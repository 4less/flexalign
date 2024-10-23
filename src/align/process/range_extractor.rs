use flexmap::values::VRange;
use kmerrs::consecutive::kmer::Kmer;

use crate::{align::{common::RangeExtractor, stats::Stats}, database::common::FlexalignDatabase, flexalign::time};

pub type Range<'a, const F: usize> = (usize, Kmer<F>, VRange<'a>, usize);

#[derive(Clone)]
pub struct StdRangeExtractor<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> {
    pub ranges: Vec<Range<'a, F>>,
    pub db: &'a D,
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> RangeExtractor<K, F> for StdRangeExtractor<'a, K, C, F, D> {
    fn retrieve(&self) -> &[Range<F>] {
        &self.ranges
    }
    
    fn generate(&mut self, kmers: &[(usize, Kmer<K>)], stats: &mut Stats) -> &[Range<F>] {
        self.ranges.clear();
        for (pos, kmer) in kmers {
            let cmer = kmer.middle::<C>();
            let fmer = kmer.flanks::<F>();


            // let (duration, result) = time(|| self.db.get_vrange(cmer.0));
            // stats.time_get_vranges += duration;

            let result = self.db.get_vrange(cmer.0);

            // eprintln!("{} {} -> {}", pos, kmer.to_string().unwrap(), result.is_some());
            let range: VRange<'a> = match result {
                Some(range) => range,
                None => continue,
            };
            let range_len = (&range).positions.len();
            self.ranges.push((*pos, fmer, range, range_len));
        }
        self.ranges.sort_unstable_by_key(|r| r.2.positions.len());
        
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

