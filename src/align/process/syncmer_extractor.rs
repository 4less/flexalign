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
//! Pass 1 walks the bases once for the s-mers of both strands and their window minima; pass 2
//! rolls the k-mers in registers and decides each with one compare. All buffers are reused.
//!
//! The result is identical to the standard extractor by construction -- the test at the bottom
//! checks that on random reads with N's, where the k-mer iterator restarts.

use bioreader::sequence::fastq_record::RefFastqRecord;
use kmerrs::consecutive::kmer::Kmer;

use crate::align::{common::KmerExtractor, stats::Stats};

#[derive(Clone)]
pub struct SyncmerKmerExtractor<const K: usize, const C: usize, const S: usize, const L: usize> {
    pub kmers: Vec<(usize, Kmer<K>)>,
    /// Forward / reverse-complement s-mer starting at each read position (2S <= 16 bits).
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

/// ASCII -> 2-bit code, 4 for anything that is not upper-case A/C/G/T. Lower case counts as N,
/// exactly as `kmerr::kmer_utils::allowed` has it, so the two selectors see the same k-mers.
static CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
};

/// Sliding-window minimum, window `w`, over `v`: `out[i] = min(v[i..i+w])` for every full window.
///
/// With AVX2 (the release build's target), 16 windows are reduced at once: the minimum of the
/// `w` shifted loads. Otherwise, van Herk / Gil-Werman: a suffix-min pass inside each block of
/// `w` and a prefix-min pass, then every window is the min of one suffix entry and one prefix
/// entry.
#[inline]
fn window_min(v: &[u16], w: usize, prefix: &mut Vec<u16>, out: &mut Vec<u16>) {
    let n = v.len();
    out.clear();
    if n < w {
        return;
    }
    let m = n - w + 1;
    out.resize(m, 0);
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        let _ = prefix;
        // SAFETY: every load reads v[i + j .. i + j + 16] with i + 16 <= m and j < w, so
        // i + j + 16 <= m - 1 + w = n; every store writes out[i..i+16] with i + 16 <= m.
        unsafe {
            use std::arch::x86_64::{__m256i, _mm256_loadu_si256, _mm256_min_epu16, _mm256_storeu_si256};
            let mut i = 0;
            while i + 16 <= m {
                let base = v.as_ptr().add(i);
                let mut acc: __m256i = _mm256_loadu_si256(base as *const __m256i);
                for j in 1..w {
                    acc = _mm256_min_epu16(acc, _mm256_loadu_si256(base.add(j) as *const __m256i));
                }
                _mm256_storeu_si256(out.as_mut_ptr().add(i) as *mut __m256i, acc);
                i += 16;
            }
            while i < m {
                let mut mn = *v.get_unchecked(i);
                for j in 1..w {
                    mn = mn.min(*v.get_unchecked(i + j));
                }
                *out.get_unchecked_mut(i) = mn;
                i += 1;
            }
        }
        return;
    }
    #[allow(unreachable_code)]
    {
        prefix.clear();
        prefix.resize(n, 0);
        out.resize(n, 0);
        let mut i = 0;
        while i < n {
            let end = (i + w).min(n);
            let mut mn = v[i];
            for j in i..end {
                mn = mn.min(v[j]);
                prefix[j] = mn;
            }
            let mut mn = v[end - 1];
            for j in (i..end).rev() {
                mn = mn.min(v[j]);
                out[j] = mn;
            }
            i = end;
        }
        for i in 0..m {
            out[i] = out[i].min(prefix[i + w - 1]);
        }
        out.truncate(m);
    }
}

impl<const K: usize, const C: usize, const S: usize, const L: usize> SyncmerKmerExtractor<K, C, S, L> {
    /// Core start inside the k-mer: `middle::<C>()` takes the central `C` bases.
    const CORE_OFFSET: usize = (K - C) / 2;
    const KMASK: u64 = (1u64 << (2 * K)) - 1;
    const KTOP: u32 = 2 * (K as u32 - 1);
    const SMASK: u32 = (1u32 << (2 * S)) - 1;
    const STOP: u32 = 2 * (S as u32 - 1);
    const CMASK: u64 = (1u64 << (2 * C)) - 1;
    /// What `Kmer::middle::<C>` shifts by (bits).
    const CSHIFT: u32 = K as u32 - C as u32;

    /// The selection proper, on a bare sequence. Pass 1 walks the bases once for the s-mers of
    /// both strands and their window minima; pass 2 rolls the forward and reverse-complement
    /// k-mers in registers (the recurrence of `KmerIter<K, true>`, restarting after any non-ACGT
    /// base, so it yields the same k-mers) and decides each with one compare.
    pub fn select(&mut self, seq: &[u8], stats: &mut Stats) -> &[(usize, Kmer<K>)] {
        debug_assert_eq!(L, C - S + 1, "L must be C - S + 1, as in ClosedSyncmer");
        debug_assert!(2 * K < 64 && 2 * S <= 16);
        self.kmers.clear();
        let n = seq.len();
        if n < K {
            return &self.kmers;
        }

        // Pass 1: s-mers of both strands at every position, then their window minima. Written
        // through pre-sized slices: a push per base re-checked capacity and re-read the Vec.
        let n_smers = n - S + 1;
        self.fwd.resize(n_smers, 0);
        self.rev.resize(n_smers, 0);
        {
            let (fwd, rev) = (self.fwd.as_mut_slice(), self.rev.as_mut_slice());
            let (mut f, mut r) = (0u32, 0u32);
            let mut run = 0usize;
            for (i, &b) in seq.iter().enumerate() {
                let code = CODE[b as usize] as u32;
                if code > 3 {
                    run = 0;
                } else {
                    run += 1;
                    f = ((f << 2) | code) & Self::SMASK;
                    r = (r >> 2) | ((3 - code) << Self::STOP);
                }
                if i + 1 >= S {
                    // 0 where the window holds a non-ACGT base: such windows lie only inside
                    // k-mers pass 2 skips, so they are never consulted.
                    let ok = run >= S;
                    // SAFETY: i + 1 - S < n - S + 1 = n_smers.
                    unsafe {
                        *fwd.get_unchecked_mut(i + 1 - S) = if ok { f as u16 } else { 0 };
                        *rev.get_unchecked_mut(i + 1 - S) = if ok { r as u16 } else { 0 };
                    }
                }
            }
        }
        window_min(&self.fwd, L, &mut self.scratch, &mut self.fwd_min);
        window_min(&self.rev, L, &mut self.scratch, &mut self.rev_min);

        // Pass 2: roll the k-mers and decide each. Everything the loop reads is a local slice
        // and the output a local Vec, so nothing is re-read from `self` per base.
        debug_assert!(self.fwd_min.len() == n_smers - L + 1);
        let (fwd, rev) = (self.fwd.as_slice(), self.rev.as_slice());
        let (fwd_min, rev_min) = (self.fwd_min.as_slice(), self.rev_min.as_slice());
        let mut kmers = std::mem::take(&mut self.kmers);
        kmers.reserve(n / 8);
        let (mut kf, mut kr) = (0u64, 0u64);
        let mut run = 0usize;
        let mut processed = 0usize;
        for (i, &b) in seq.iter().enumerate() {
            let code = CODE[b as usize] as u64;
            if code > 3 {
                run = 0;
                continue;
            }
            run += 1;
            kf = ((kf << 2) | code) & Self::KMASK;
            kr = (kr >> 2) | ((3 - code) << Self::KTOP);
            if run < K {
                continue;
            }
            processed += 1;
            let pos = i + 1 - K;
            let cmer_fwd = (kf >> Self::CSHIFT) & Self::CMASK;
            let cmer_rev = (kr >> Self::CSHIFT) & Self::CMASK;
            let p = pos + Self::CORE_OFFSET;
            let forward = cmer_fwd < cmer_rev;
            // In range by construction: p + L - 1 <= i + 1 - S < n_smers, and p <= n_smers - L
            // (the last full window), since p + L - 1 <= i + 1 - S.
            debug_assert!(p + L - 1 < n_smers && p < fwd_min.len());
            let selected = unsafe {
                if forward {
                    *fwd.get_unchecked(p) == *fwd_min.get_unchecked(p)
                } else {
                    *rev.get_unchecked(p + L - 1) == *rev_min.get_unchecked(p)
                }
            };
            if selected {
                kmers.push((pos, Kmer::<K>(if forward { kf } else { kr })));
            }
        }
        self.kmers = kmers;
        stats.kmers_processed += processed;
        stats.minimizer += self.kmers.len();
        &self.kmers
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

#[cfg(test)]
mod tests {
    use super::*;
    use kmerrs::consecutive::kmer::KmerIter;
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
                if n_rate > 0 && x % n_rate == 0 { b"NRYacgt"[(x >> 16) as usize % 7] } else { b"ACGT"[(x >> 8) as usize % 4] }
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
