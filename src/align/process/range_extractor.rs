
use std::cmp::min;

use bioreader::sequence::fastq_record::RefFastqRecord;
use flexmap::range::AnyRange;
use kmerrs::consecutive::kmer::Kmer;
use kmerrs::kmer_utils::mask_coremer;

use crate::{align::{common::{RangeExtractor, UniqueHit}, process::seed_extractor::MAX_KEPT_TIES, stats::Stats}, database::common::FlexalignDatabase};

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
    /// Per-flank distance scratch of the fast path's header resolution.
    flex_scratch: Vec<u8>,
    /// The keys `prefetch` fetches ahead, reused across reads.
    prefetch_keys: Vec<u64>,
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> RangeExtractor<K, F> for StdRangeExtractor<'a, K, C, F, D> {
    fn retrieve(&self) -> &[Range<'_, F>] {
        &self.ranges
    }

    fn prefetch(&mut self, kmers: &[(usize, Kmer<K>)]) {
        self.prefetch_keys.clear();
        self.prefetch_keys.extend(kmers.iter().map(|(_, kmer)| mask_coremer::<C>(kmer.middle::<C>().0)));
        self.db.prefetch_keys(&self.prefetch_keys);
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
        
        // let (duration, _) = stage_time(|| self.ranges.sort_unstable_by_key(|r| r.2.positions.len()));
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
        self.hits_of(max_positions, out);
        true
    }

    fn unique_hits(&mut self, kmers: &[(usize, Kmer<K>)], max_positions: usize, out: &mut Vec<UniqueHit>) -> bool {
        self.lookup(kmers);
        self.prefetched = true;
        out.clear();
        self.hits_of(max_positions, out);
        true
    }
}

impl<'a, const K: usize, const C: usize, const F: usize, D: FlexalignDatabase> StdRangeExtractor<'a, K, C, F, D> {
    /// The hits the fast path votes with. A headerless key (at most HEADER_THRESHOLD positions)
    /// gives every position; a headered key with at most `FAST_HIT_MAX_CELLS` positions gives
    /// the cells whose flanks match the read's exactly (a full K-mer match), resolved the way the
    /// seed extractor resolves them. A 15-mer core occurs several times in a few Gbp of
    /// reference, so on anything but a small reference nearly every key is headered: without
    /// these the path placed 1% of pairs on a 5 GB reference, against 78% on 1 GB.
    fn hits_of(&mut self, max_positions: usize, out: &mut Vec<UniqueHit>) {
        const FAST_HIT_MAX_CELLS: usize = 64;
        let mut idx = [0u32; MAX_KEPT_TIES];
        for r in &self.ranges {
            let n = r.vrange.len();
            if !r.vrange.has_header() {
                if n <= max_positions {
                    for j in 0..n {
                        let (rval, rpos) = r.vrange.value(j);
                        out.push(UniqueHit { qpos: r.qpos as u32, rval, rpos });
                    }
                }
            } else if n <= FAST_HIT_MAX_CELLS {
                if let Some(kept) = r.exact_flank_cells(max_positions.max(1), &mut self.flex_scratch, &mut idx) {
                    for &i in &idx[..kept] {
                        let (rval, rpos) = r.vrange.value(i as usize);
                        out.push(UniqueHit { qpos: r.qpos as u32, rval, rpos });
                    }
                }
            }
        }
    }
}

impl<'a, const F: usize> Range<'a, F> {
    /// The cells whose flank equals the read's: what `resolve_flex(max, 0, ..)` keeps when the
    /// best flank distance is 0, found without the distances. Distance 0 means every flank base
    /// matches, so the 32-bit flank cell equals the read's flank exactly. None when no cell
    /// matches or more than `max` do. `out` holds at least `max` (<= `MAX_KEPT_TIES`) indices.
    #[inline]
    pub fn exact_flank_cells(&self, max: usize, dist_buf: &mut Vec<u8>, out: &mut [u32]) -> Option<usize> {
        let AnyRange::Packed(p) = &self.vrange else {
            return match self.resolve_flex(max, 0, dist_buf, out) {
                Some((0, kept)) if kept > 0 => Some(kept),
                _ => None,
            };
        };
        let flex = self.fmer.0 as u32;
        let mut count = 0usize;
        for (i, c) in p.flex_bytes().chunks_exact(4).enumerate() {
            if u32::from_le_bytes([c[0], c[1], c[2], c[3]]) == flex {
                if count == max {
                    return None;
                }
                out[count] = i as u32;
                count += 1;
            }
        }
        (count > 0).then_some(count)
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
        Self { ranges: Vec::new(), db, keys: Vec::new(), key_info: Vec::new(), order: Vec::new(), sorted: Vec::new(), prefetched: false, subset: Vec::new(), flex_scratch: Vec::new(), prefetch_keys: Vec::new() }
    }
}

