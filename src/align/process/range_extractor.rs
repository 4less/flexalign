
use std::cmp::min;

use bioreader::sequence::fastq_record::RefFastqRecord;
use flexmap::range::AnyRange;
use kmerrs::consecutive::kmer::Kmer;
use kmerrs::kmer_utils::mask_coremer;

use crate::{align::{common::{RangeExtractor, UniqueHit}, stats::Stats}, database::common::FlexalignDatabase};

// pub type Range<'a, const F: usize> = (usize, Kmer<F>, VRange<'a>, usize);

#[derive(Clone)]
pub struct Range<'a, const F: usize> {
    pub qpos: usize,
    pub fmer: Kmer<F>,
    pub vrange: AnyRange<'a>,
    pub range_size: u32,
    pub core_min_phredq: u16,
    pub flex_min_phredq: u16,
}

#[derive(Clone)]
pub struct StdRangeExtractor<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> {
    pub ranges: Vec<Range<'a, F>>,
    pub db: &'a D,
    /// The read's keys and their (query position, flank), reused across reads.
    keys: Vec<u64>,
    key_info: Vec<(usize, Kmer<F>)>,
    /// (range size, index) sort keys and the reordered ranges, reused across reads.
    order: Vec<(u32, u32)>,
    sorted: Vec<Range<'a, F>>,
    /// Set by `unique_hits`, which already fetched every range of these keys into `ranges`:
    /// `generate` on the same k-mers then skips the lookup (the fast path's fallback pairs).
    prefetched: bool,
    /// The strided round's k-mers, reused across reads.
    subset: Vec<(usize, Kmer<K>)>,
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> RangeExtractor<K, F> for StdRangeExtractor<'a, K, C, F, D> {
    fn retrieve(&self) -> &[Range<'_, F>] {
        &self.ranges
    }
    
    fn generate(&mut self, kmers: &[(usize, Kmer<K>)], rec: &RefFastqRecord, stats: &mut Stats) -> &[Range<'_, F>] {
        let reuse = self.prefetched && self.key_info.len() == kmers.len()
            && self.key_info.iter().zip(kmers).all(|(k, (pos, _))| k.0 == *pos);
        self.prefetched = false;
        if !reuse {
            self.lookup(kmers);
        }
        
        for range in &mut self.ranges {
            let qual_mer = &rec.qual()[range.qpos..(range.qpos + C)];
            range.core_min_phredq = min(qual_mer.iter().map(|q| q - 33).min().unwrap() as u16, 40);
            // eprintln!("{} {} ... sort by {}", range.core_min_phredq, range.vrange.positions.len(), range.vrange.positions.len() * (41 - range.core_min_phredq as usize));
        }

        // Rarest first. Sorted through (len, index) keys rather than the ~70-byte ranges
        // themselves (the element moves were a visible share of a run); ties keep query order.
        self.order.clear();
        self.order.extend(self.ranges.iter().enumerate().map(|(i, r)| (r.vrange.len() as u32, i as u32)));
        self.order.sort_unstable();
        self.sorted.clear();
        self.sorted.extend(self.order.iter().map(|&(_, i)| self.ranges[i as usize].clone()));
        std::mem::swap(&mut self.ranges, &mut self.sorted);
        
        // let (duration, _) = time(|| self.ranges.sort_unstable_by_key(|r| r.2.positions.len()));
        // stats.time_range_sorting += duration;
        stats.ranges += self.ranges.len();

        &self.ranges
    }

    fn unique_hits_strided(&mut self, kmers: &[(usize, Kmer<K>)], stride: usize, max_positions: usize, out: &mut Vec<UniqueHit>) -> bool {
        self.subset.clear();
        self.subset.extend(kmers.iter().step_by(stride.max(1)).cloned());
        let subset = std::mem::take(&mut self.subset);
        self.lookup(&subset);
        self.subset = subset;
        self.prefetched = false;
        out.clear();
        Self::hits_of(&self.ranges, max_positions, out);
        true
    }

    fn unique_hits(&mut self, kmers: &[(usize, Kmer<K>)], max_positions: usize, out: &mut Vec<UniqueHit>) -> bool {
        self.lookup(kmers);
        self.prefetched = true;
        out.clear();
        Self::hits_of(&self.ranges, max_positions, out);
        true
    }
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> StdRangeExtractor<'a, K, C, F, D> {
    fn hits_of(ranges: &[Range<'a, F>], max_positions: usize, out: &mut Vec<UniqueHit>) {
        for r in ranges {
            // Headerless keys (at most HEADER_THRESHOLD positions) hold the positions directly.
            let n = r.vrange.len();
            if n <= max_positions && !r.vrange.has_header() {
                for j in 0..n {
                    let (rval, rpos) = r.vrange.value(j);
                    out.push(UniqueHit { qpos: r.qpos as u32, rval, rpos });
                }
            }
        }
    }
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> StdRangeExtractor<'a, K, C, F, D> {
    /// Fetches every key's range into `ranges` (unsorted, qualities unset).
    fn lookup(&mut self, kmers: &[(usize, Kmer<K>)]) {
        self.ranges.clear();
        self.keys.clear();
        self.key_info.clear();

        for (pos, kmer) in kmers {
            // Look up under the same de-biasing mix used at index build and syncmer selection.
            self.keys.push(mask_coremer::<C>(kmer.middle::<C>().0));
            self.key_info.push((*pos, kmer.flanks::<F>()));
        }

        // All of the read's keys in one call, so that the index can fetch them ahead.
        let db: &'a D = self.db;
        let (ranges, key_info) = (&mut self.ranges, &self.key_info);
        db.for_each_vrange(&self.keys, &mut |i, range| {
            if let Some(range) = range {
                let (qpos, fmer) = (key_info[i].0, key_info[i].1.clone());
                let range_size = range.len() as u32;
                ranges.push(Range { qpos, fmer, vrange: range, range_size, core_min_phredq: 0, flex_min_phredq: 0 });
            }
        });
    }
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> StdRangeExtractor<'a, K, C, F, D> {
    pub fn new(db: &'a D) -> Self {
        Self { ranges: Vec::new(), db, keys: Vec::new(), key_info: Vec::new(), order: Vec::new(), sorted: Vec::new(), prefetched: false, subset: Vec::new() }
    }
}
