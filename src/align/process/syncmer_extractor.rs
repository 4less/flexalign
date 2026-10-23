//! Syncmer selection in O(1) per k-mer position.
//!
//! [`StdKmerExtractor`](super::kmer_extractor::StdKmerExtractor) asks `ClosedSyncmer::is_minimizer`
//! about every k-mer: that unpacks the canonical core-mer into its `L = C - S + 1` s-mers and takes
//! their minimum, ~25 operations per position, and it was the single hottest function of a run
//! (12% of CPU, with the extractor around it another 10%).
//!
//! The selection rule is: the core-mer is kept iff its FIRST s-mer is the minimum of all its s-mers
//! (an open syncmer, offset 0, with ties kept -- see `kmerr`'s `closed_syncmer.rs`, where the
//! `last()` half of the closed rule is commented out). Every one of those s-mers is a 2-bit word
//! over `S` consecutive read bases, in either the forward or the reverse-complement strand, so a
//! read's s-mer values can be computed once per strand, and "minimum over the core's `L` s-mers"
//! is a sliding-window minimum over `L` consecutive entries of one of those two arrays. Which
//! strand, and which entry must equal the minimum, follows from which orientation of the core is
//! canonical:
//!
//! * forward canonical: the core's s-mer `i` is the forward s-mer at read position `p + i`, so
//!   selected iff `F[p] == min F[p..p+L]`;
//! * reverse canonical: the core's s-mer `i` is the reverse complement of the forward s-mer at
//!   read position `p + L - 1 - i`, so its first s-mer is `R[p + L - 1]`, and selected iff
//!   `R[p + L - 1] == min R[p..p+L]`,
//!
//! where `p` is the core's start in the read, `F[q]` the 2-bit value of `read[q..q+S]`, and
//! `R[q]` that of its reverse complement. The window minima come from the van Herk / Gil-Werman
//! block decomposition (two passes, three compares per entry, no branches on data).
//!
//! The result is identical to the standard extractor by construction -- the test at the bottom
//! checks that on random reads with N's, where the k-mer iterator restarts.

use bioreader::sequence::fastq_record::RefFastqRecord;
use kmerrs::consecutive::kmer::{Kmer, KmerIter};
use kmerrs::kmer_utils::NUC_U8_TO_BIT;

use crate::align::{common::KmerExtractor, stats::Stats};

#[derive(Clone)]
pub struct SyncmerKmerExtractor<const K: usize, const C: usize, const S: usize, const L: usize> {
    pub kmers: Vec<(usize, Kmer<K>)>,
    fwd: Vec<u16>,
    rev: Vec<u16>,
    fwd_min: Vec<u16>,
    rev_min: Vec<u16>,
    scratch: Vec<u16>,
}

impl<const K: usize, const C: usize, const S: usize, const L: usize> Default
    for SyncmerKmerExtractor<K, C, S, L>
{
    fn default() -> Self {
        Self { kmers: Vec::new(), fwd: Vec::new(), rev: Vec::new(), fwd_min: Vec::new(), rev_min: Vec::new(), scratch: Vec::new() }
    }
}

/// Sliding-window minimum, window `w`, over `v`: `out[i] = min(v[i..i+w])` for every full window.
/// Van Herk / Gil-Werman: a suffix-min pass inside each block of `w` and a prefix-min pass, then
/// every window is the min of one suffix entry and one prefix entry.
#[inline]
fn window_min(v: &[u16], w: usize, prefix: &mut Vec<u16>, out: &mut Vec<u16>) {
    let n = v.len();
    out.clear();
    if n < w {
        return;
    }
    prefix.clear();
    prefix.resize(n, 0);
    // suffix minima, written into `out` first (it is resized to n, then truncated to n - w + 1)
    out.resize(n, 0);
    let mut i = 0;
    while i < n {
        let end = (i + w).min(n);
        // prefix within the block
        let mut m = v[i];
        for j in i..end {
            m = m.min(v[j]);
            prefix[j] = m;
        }
        // suffix within the block
        let mut m = v[end - 1];
        for j in (i..end).rev() {
            m = m.min(v[j]);
            out[j] = m;
        }
        i = end;
    }
    for i in 0..=(n - w) {
        out[i] = out[i].min(prefix[i + w - 1]);
    }
    out.truncate(n - w + 1);
}

impl<const K: usize, const C: usize, const S: usize, const L: usize> SyncmerKmerExtractor<K, C, S, L> {
    /// Core start inside the k-mer: `middle::<C>()` takes the central `C` bases.
    const CORE_OFFSET: usize = (K - C) / 2;

    /// Fills `fwd`/`rev` with the s-mer value at every read position (0 where the window holds a
    /// non-ACGT base; such windows lie only inside k-mers the iterator skips, so they are never
    /// consulted) and the window minima over `L` entries.
    fn prepare(&mut self, seq: &[u8]) {
        debug_assert_eq!(L, C - S + 1, "L must be C - S + 1, as in ClosedSyncmer");
        let n = seq.len();
        self.fwd.clear();
        self.rev.clear();
        if n < S {
            self.fwd_min.clear();
            self.rev_min.clear();
            return;
        }
        let mask: u32 = (1u32 << (2 * S)) - 1;
        let top = 2 * (S - 1);
        let mut f: u32 = 0;
        let mut r: u32 = 0;
        let mut valid = 0usize; // consecutive ACGT bases ending at the current position
        for (i, &b) in seq.iter().enumerate() {
            let code = if b >= b'A' && b < b'A' + 20 { NUC_U8_TO_BIT[(b - b'A') as usize] } else { 4 };
            if code > 3 {
                valid = 0;
                f = 0;
                r = 0;
            } else {
                valid += 1;
                f = ((f << 2) | code as u32) & mask;
                r = (r >> 2) | ((3 - code as u32) << top);
            }
            if i + 1 >= S {
                let ok = valid >= S;
                self.fwd.push(if ok { f as u16 } else { 0 });
                self.rev.push(if ok { r as u16 } else { 0 });
            }
        }
        window_min(&self.fwd, L, &mut self.scratch, &mut self.fwd_min);
        window_min(&self.rev, L, &mut self.scratch, &mut self.rev_min);
    }
}

impl<const K: usize, const C: usize, const S: usize, const L: usize> KmerExtractor<K>
    for SyncmerKmerExtractor<K, C, S, L>
{
    fn generate(&mut self, rec: &RefFastqRecord, stats: &mut Stats) -> &[(usize, Kmer<K>)] {
        self.select(rec.seq(), stats)
    }

    fn retrieve(&self) -> &[(usize, Kmer<K>)] {
        &self.kmers
    }
}

impl<const K: usize, const C: usize, const S: usize, const L: usize> SyncmerKmerExtractor<K, C, S, L> {
    /// The selection proper, on a bare sequence.
    pub fn select(&mut self, seq: &[u8], stats: &mut Stats) -> &[(usize, Kmer<K>)] {
        self.kmers.clear();
        self.prepare(seq);
        let (fwd, rev, fwd_min, rev_min) = (&self.fwd, &self.rev, &self.fwd_min, &self.rev_min);
        let mut processed = 0usize;
        for (pos, kmer_fwd, kmer_rev) in KmerIter::<K, true>::new(seq) {
            processed += 1;
            let cmer_fwd = kmer_fwd.middle::<C>();
            let cmer_rev = kmer_rev.middle::<C>();
            let p = pos + Self::CORE_OFFSET;
            let selected = if cmer_fwd < cmer_rev {
                fwd[p] == fwd_min[p]
            } else {
                rev[p + L - 1] == rev_min[p]
            };
            if !selected {
                continue;
            }
            let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
            self.kmers.push((pos, kmer));
        }
        stats.kmers_processed += processed;
        stats.minimizer += self.kmers.len();
        &self.kmers
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kmerrs::minimizer::context_free::Minimizer;
    use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

    const K: usize = 31;
    const C: usize = 15;
    const S: usize = 4;
    const L: usize = C - S + 1;

    /// xorshift: no dependency, reproducible.
    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    fn random_read(seed: &mut u64, len: usize, n_rate: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                let x = rng(seed);
                if n_rate > 0 && x % n_rate == 0 { b'N' } else { b"ACGT"[(x >> 8) as usize % 4] }
            })
            .collect()
    }

    /// What `StdKmerExtractor` does, verbatim, minus the record plumbing.
    fn reference(seq: &[u8], out: &mut Vec<(usize, Kmer<K>)>, processed: &mut usize) {
        let mut minimizer = ClosedSyncmer::<C, S, L>::new();
        out.clear();
        for (pos, kmer_fwd, kmer_rev) in KmerIter::<K, true>::new(seq) {
            *processed += 1;
            let cmer_fwd = kmer_fwd.middle::<C>();
            let cmer_rev = kmer_rev.middle::<C>();
            let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
            let cmer = std::cmp::min(cmer_fwd, cmer_rev);
            if !minimizer.is_minimizer(cmer.0) {
                continue;
            }
            out.push((pos, kmer));
        }
    }

    #[test]
    fn matches_std_extractor() {
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut fast = SyncmerKmerExtractor::<K, C, S, L>::default();
        let mut stats = Stats::default();
        let mut expected = Vec::new();
        let mut processed = 0usize;
        let mut selected = 0;
        for i in 0..3000 {
            let len = 20 + (rng(&mut seed) % 300) as usize;
            let n_rate = match i % 3 { 0 => 0, 1 => 40, _ => 7 };
            let seq = random_read(&mut seed, len, n_rate);
            let a = fast.select(&seq, &mut stats).to_vec();
            reference(&seq, &mut expected, &mut processed);
            assert_eq!(a, expected, "read {i}: {}", String::from_utf8_lossy(&seq));
            selected += a.len();
        }
        assert!(selected > 0);
        assert_eq!(stats.kmers_processed, processed);
        assert_eq!(stats.minimizer, selected);
    }

    #[test]
    fn window_min_matches_naive() {
        let mut seed = 42u64;
        let (mut scratch, mut out) = (Vec::new(), Vec::new());
        for _ in 0..500 {
            let n = (rng(&mut seed) % 60) as usize;
            let w = 1 + (rng(&mut seed) % 15) as usize;
            let v: Vec<u16> = (0..n).map(|_| (rng(&mut seed) % 256) as u16).collect();
            window_min(&v, w, &mut scratch, &mut out);
            if n < w {
                assert!(out.is_empty());
                continue;
            }
            for i in 0..=(n - w) {
                assert_eq!(out[i], *v[i..i + w].iter().min().unwrap());
            }
        }
    }
}
