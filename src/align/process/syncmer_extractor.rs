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
        let mut done = 0usize;
        #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
        if S == 4 {
            // 16 s-mers a step: codes of 19 bases via a nibble shuffle, validated against the four
            // letters, widened to 16-bit lanes and combined from four shifted loads.
            done = unsafe { smers4_avx2(seq, self.fwd.as_mut_slice(), self.rev.as_mut_slice()) };
        }
        if done < n_smers {
            let (fwd, rev) = (self.fwd.as_mut_slice(), self.rev.as_mut_slice());
            let (mut f, mut r) = (0u32, 0u32);
            let mut run = 0usize;
            // Resume the recurrence at base `done`: the s-mer at `done` needs bases done..done+S.
            for (i, &b) in seq.iter().enumerate().skip(done) {
                let code = CODE[b as usize] as u32;
                if code > 3 {
                    run = 0;
                } else {
                    run += 1;
                    f = ((f << 2) | code) & Self::SMASK;
                    r = (r >> 2) | ((3 - code) << Self::STOP);
                }
                if i + 1 >= S + done {
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
        kmers.clear();
        // At most one k-mer a base. Every k-mer is written to the next free slot and kept only
        // if selected: the strand (canonical orientation) and the selection are coin flips to
        // the branch predictor, so both tests are computed and picked, not branched on.
        kmers.reserve(n);
        let out = kmers.spare_capacity_mut();
        let mut kept = 0usize;
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
            let (sel_fwd, sel_rev) = unsafe {
                (*fwd.get_unchecked(p) == *fwd_min.get_unchecked(p), *rev.get_unchecked(p + L - 1) == *rev_min.get_unchecked(p))
            };
            let selected = (forward & sel_fwd) | (!forward & sel_rev);
            let kmer = if forward { kf } else { kr };
            // SAFETY: kept <= processed - 1 < n = the reserved capacity.
            unsafe { out.get_unchecked_mut(kept).write((pos, Kmer::<K>(kmer))) };
            kept += selected as usize;
        }
        // SAFETY: the first `kept` slots were written above.
        unsafe { kmers.set_len(kept) };
        self.kmers = kmers;
        stats.kmers_processed += processed;
        stats.minimizer += self.kmers.len();
        &self.kmers
    }
}

/// The forward and reverse-complement 4-mer value at 16 positions a step, as long as 32 bases
/// remain, 0 where a window holds a non-ACGT base. Returns how many
/// positions were written; the caller finishes the tail with the scalar recurrence.
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[target_feature(enable = "avx2")]
unsafe fn smers4_avx2(seq: &[u8], fwd: &mut [u16], rev: &mut [u16]) -> usize {
    use std::arch::x86_64::*;
    let n = seq.len();
    if n < 32 {
        return 0;
    }
    // low nibble -> code: A(0x41)->0, C(0x43)->1, G(0x47)->2, T(0x54)->3; anything else is
    // caught by the exact-letter comparison below.
    let table = _mm_setr_epi8(0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
    let nib = _mm_set1_epi8(0x0f);
    let (la, lc, lg, lt) = (_mm_set1_epi8(b'A' as i8), _mm_set1_epi8(b'C' as i8), _mm_set1_epi8(b'G' as i8), _mm_set1_epi8(b'T' as i8));
    let three = _mm256_set1_epi16(3);
    let mut p = 0usize;
    while p + 32 <= n {
        // bases p..p+16 and p+16..p+32 as codes (u8) and validity (0xff where ACGT); the views
        // at offsets 1, 2, 3 come from the concatenation (alignr takes the LOW half second).
        let raw0 = _mm_loadu_si128(seq.as_ptr().add(p) as *const __m128i);
        let raw1 = _mm_loadu_si128(seq.as_ptr().add(p + 16) as *const __m128i);
        let code = |v: __m128i| _mm_shuffle_epi8(table, _mm_and_si128(v, nib));
        let valid = |v: __m128i| _mm_or_si128(_mm_or_si128(_mm_cmpeq_epi8(v, la), _mm_cmpeq_epi8(v, lc)), _mm_or_si128(_mm_cmpeq_epi8(v, lg), _mm_cmpeq_epi8(v, lt)));
        let (c0, v0) = (code(raw0), valid(raw0));
        let (ch, vh) = (code(raw1), valid(raw1));
        let c1 = _mm_alignr_epi8(ch, c0, 1);
        let c2 = _mm_alignr_epi8(ch, c0, 2);
        let c3 = _mm_alignr_epi8(ch, c0, 3);
        let v1 = _mm_alignr_epi8(vh, v0, 1);
        let v2 = _mm_alignr_epi8(vh, v0, 2);
        let v3 = _mm_alignr_epi8(vh, v0, 3);
        // widen to 16 lanes of u16
        let w = |v: __m128i| _mm256_cvtepu8_epi16(v);
        let (w0, w1, w2, w3) = (w(c0), w(c1), w(c2), w(c3));
        let f = _mm256_or_si256(_mm256_or_si256(_mm256_slli_epi16(w0, 6), _mm256_slli_epi16(w1, 4)), _mm256_or_si256(_mm256_slli_epi16(w2, 2), w3));
        let r = _mm256_or_si256(_mm256_or_si256(_mm256_slli_epi16(_mm256_sub_epi16(three, w3), 6), _mm256_slli_epi16(_mm256_sub_epi16(three, w2), 4)), _mm256_or_si256(_mm256_slli_epi16(_mm256_sub_epi16(three, w1), 2), _mm256_sub_epi16(three, w0)));
        let ok = _mm256_cvtepi8_epi16(_mm_and_si128(_mm_and_si128(v0, v1), _mm_and_si128(v2, v3)));
        _mm256_storeu_si256(fwd.as_mut_ptr().add(p) as *mut __m256i, _mm256_and_si256(f, ok));
        _mm256_storeu_si256(rev.as_mut_ptr().add(p) as *mut __m256i, _mm256_and_si256(r, ok));
        p += 16;
    }
    p
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
