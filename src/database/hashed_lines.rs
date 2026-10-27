//! A reader of flexmap's line-hashed key table (`HashedKeys`) built from its public sections and
//! parameters. It answers as `HashedKeys::locate` does, but compares a line's 8 entries at once
//! (AVX2) instead of one after the other with a branch each, and hands a read's keys out in chunks:
//! locate a chunk and fetch its values, then call back, so the values arrive while the chunk is
//! located. The table layout is flexmap's (see `HashedKeys`); only the reading is done here.

use flexmap::memory::{prefetch, Section};
use flexmap::packed::{HashedKeys, PackedFlexmap, PackedRange, SlotLocator};

const ENTRIES: u64 = 8;
const DISP_BITS: u32 = 5;
const COUNT_BITS: u32 = 13;
/// In a line's first entry: keys homed at or before this line lie past it.
const OVERFLOW: u64 = 1 << 63;
/// Keys located before their values are read.
const CHUNK: usize = 32;
/// How far ahead of the one being located a key's line is fetched.
const KEY_AHEAD: usize = 8;

#[derive(Clone)]
pub struct HashedLines {
    lines: Section,
    homes: u64,
    /// A key's home line is the high word of key * mult (flexmap's `HashedKeys`); the key space need
    /// not be a power of two (a shard's is its key range).
    mult: u64,
    tag_bits: u32,
    start_bits: u32,
}

impl HashedLines {
    pub fn new(keys: &HashedKeys) -> Self {
        let [n_keys, _used, homes, start_bits] = keys.params();
        // As HashedKeys derives them: mult = floor(homes * 2^64 / keys), and enough tag bits for the
        // ceil(2^64 / mult) keys a home can hold.
        let n_keys = n_keys.max(homes + 1);
        let mult = (((homes as u128) << 64) / n_keys as u128) as u64;
        let per_home = (((1u128 << 64) + mult as u128 - 1) / mult as u128) as u64;
        let tag_bits = (64 - per_home.leading_zeros()).max(1);
        let lines = keys.sections().swap_remove(0);
        Self { lines, homes, mult, tag_bits, start_bits: start_bits as u32 }
    }

    #[inline(always)]
    fn line(&self, index: u64) -> *const u64 {
        unsafe { (self.lines.ptr() as *const u64).add((index * ENTRIES) as usize) }
    }

    #[inline(always)]
    fn home(&self, key: u64) -> u64 {
        ((key as u128 * self.mult as u128) >> 64) as u64
    }

    #[inline(always)]
    pub fn prefetch_key(&self, key: u64) {
        prefetch(self.line(self.home(key)) as *const u8);
    }

    /// (first cell, cells) of `key`, as `HashedKeys::locate`.
    #[inline(always)]
    pub fn locate(&self, key: u64) -> Option<(u64, u64)> {
        let home = self.home(key);
        let tag = key & ((1 << self.tag_bits) - 1);
        let low = self.tag_bits + DISP_BITS;
        let low_mask = (1u64 << low) - 1;
        for disp in 0..(1u64 << DISP_BITS) {
            let p = self.line(home + disp);
            if let Some(i) = match_in_line(p, low_mask, tag | disp << self.tag_bits) {
                let e = unsafe { *p.add(i) } & !OVERFLOW;
                let start = (e >> low) & ((1u64 << self.start_bits) - 1);
                return Some((start, (e >> (low + self.start_bits)) & ((1 << COUNT_BITS) - 1)));
            }
            if unsafe { *p } & OVERFLOW == 0 {
                return None;
            }
        }
        None
    }

    /// `f(i, index.get(keys[i]))` for every key in order.
    #[inline]
    pub fn for_each<'s, const C: usize, const F: usize, const T: usize>(
        &self,
        index: &'s PackedFlexmap<HashedKeys, C, F, T>,
        keys: &[u64],
        mut f: impl FnMut(usize, Option<PackedRange<'s>>),
    ) {
        let mut loc = [None::<(u64, u64)>; CHUNK];
        for (c, chunk) in keys.chunks(CHUNK).enumerate() {
            let base = c * CHUNK;
            for (j, &key) in chunk.iter().enumerate() {
                if let Some(&ahead) = keys.get(base + j + KEY_AHEAD) {
                    self.prefetch_key(ahead);
                }
                loc[j] = self.locate(key);
                if let Some((start, slots)) = loc[j] {
                    index.values.prefetch(start, slots);
                }
            }
            for (j, l) in loc[..chunk.len()].iter().enumerate() {
                f(base + j, l.map(|(start, slots)| index.values.range(start, slots)));
            }
        }
    }
}

/// The entry of the line at `p` whose low bits are `want`, if any (at most one is: a line's
/// entries of one displacement have distinct tags). An empty entry (0, but for the overflow bit)
/// never matches.
#[inline(always)]
fn match_in_line(p: *const u64, low_mask: u64, want: u64) -> Option<usize> {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    unsafe {
        use core::arch::x86_64::*;
        let m = _mm256_set1_epi64x(low_mask as i64);
        let w = _mm256_set1_epi64x(want as i64);
        let body = _mm256_set1_epi64x(!OVERFLOW as i64);
        let zero = _mm256_setzero_si256();
        #[inline(always)]
        unsafe fn hits(v: __m256i, m: __m256i, w: __m256i, body: __m256i, zero: __m256i) -> u32 {
            let eq = _mm256_cmpeq_epi64(_mm256_and_si256(v, m), w);
            let empty = _mm256_cmpeq_epi64(_mm256_and_si256(v, body), zero);
            _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_andnot_si256(empty, eq))) as u32
        }
        let a = _mm256_loadu_si256(p as *const __m256i);
        let b = _mm256_loadu_si256(p.add(4) as *const __m256i);
        let bits = hits(a, m, w, body, zero) | hits(b, m, w, body, zero) << 4;
        if bits == 0 { None } else { Some(bits.trailing_zeros() as usize) }
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
    {
        (0..ENTRIES as usize).find(|&i| {
            let e = unsafe { *p.add(i) } & !OVERFLOW;
            e & low_mask == want && e != 0
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locates_as_flexmap() {
        // 2^14 possible keys (and a shard's key space, not a power of two), every fifth one used,
        // and all of a dense block, whose homes overflow into the following lines.
        for keys in [1u64 << 14, 12_345] {
        let used: Vec<u64> = (0..keys).filter(|k| k % 5 == 0 || (4000..4200).contains(k)).collect();
        let each_key = |f: &mut dyn FnMut(u64, u64, u64)| {
            let mut start = 0;
            for (i, &k) in used.iter().enumerate() {
                let slots = 1 + (i as u64 % 7);
                f(k, start, slots);
                start += slots;
            }
        };
        let total: u64 = (0..used.len() as u64).map(|i| 1 + i % 7).sum();
        let table = HashedKeys::build(keys, used.len() as u64, total, &each_key).expect("fits");
        let lines = HashedLines::new(&table);
        assert!((0..lines.homes).any(|l| unsafe { *lines.line(l) } & OVERFLOW != 0), "no line overflows");
        for k in 0..keys {
            assert_eq!(lines.locate(k), table.locate(k), "key {k}");
        }
        }
    }
}
