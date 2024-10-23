use std::cmp::min;

use bioreader::sequence::fastq_record::RefFastqRecord;
use kmerrs::{consecutive::kmer::{Kmer, KmerIter}, minimizer::context_free::Minimizer};

use crate::align::{common::KmerExtractor, stats::Stats};

#[derive(Clone)]
pub struct StdKmerExtractor<const K: usize, const C: usize, M: Minimizer + Default> {
    pub kmers: Vec<(usize, Kmer<K>)>,
    pub minimizer: M,
}

impl<const K: usize, const C: usize, M: Minimizer + Default> 
        Default for StdKmerExtractor<K, C, M> {
    fn default() -> Self {
        Self { kmers: Vec::new(), minimizer: Default::default() }
    }
}

impl<
        const K: usize, 
        const C: usize,
        M: Minimizer + Default
    > KmerExtractor<K> for StdKmerExtractor<K, C, M> {
    fn generate(&mut self, rec: &RefFastqRecord, stats: &mut Stats) -> &[(usize, Kmer<K>)] {
        let iter = KmerIter::<K, true>::new(rec.seq());
        self.kmers.clear();
        for (pos, kmer_fwd, kmer_rev) in iter {
            stats.kmers_processed += 1;

            let cmer_fwd = kmer_fwd.middle::<C>();
            let cmer_rev = kmer_rev.middle::<C>();
            let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
            let cmer = min(cmer_fwd, cmer_rev);


            // let (duration, is_minimizer) = time(|| );
            // stats.time_get_minimizer += duration;
            // timing the minimizer takes like 5 more seconds.
            if !self.minimizer.is_minimizer(cmer.0) {
                continue;
            };

            stats.minimizer += 1;

            self.kmers.push((pos, kmer));
        }

        &self.kmers

    }

    fn retrieve(&self) -> &[(usize, Kmer<K>)] {
        &self.kmers
    }
}
