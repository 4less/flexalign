//! The repeat BWT: an FM-index over the reference regions dense in cores the main index drops.
//!
//! The main index drops a core with more than `--max-range-size` positions (young LINE-1, Alu,
//! satellites). A read from such a region has most of its cores dropped, and the few left often
//! point to other copies. The repeat index (`repeat.rs`) keeps every dropped position keyed by its
//! 31-mer -- one entry per position however long the run, and nothing longer than 31 bases to tell
//! copies apart.
//!
//! Here the regions themselves are indexed instead: an FM-index answers "where does this read's
//! sequence occur" for a match of any length. A read's match is extended base by base along the
//! read until it is unique (or the read ends, or it hits a difference), so copies are told apart by
//! as much of the read as there is; a match still found in many copies is counted without listing
//! them, and positions are only produced once few remain.
//!
//! Only where it pays (`--bwt-min-dropped`, m): a region enters the index around every dropped
//! position with at least m dropped positions within one read length, and a read is searched only
//! when it carries at least m dropped cores -- so every read that is searched lies in the index. A
//! read with fewer is placed by its other cores, as without the BWT.
//!
//! Layout, for memory: the text is the forward strand of those regions, split at bases other than
//! ACGT (a read is searched in both orientations); the BWT is 2 bits a row with rank counts inline
//! (64-byte blocks of 192 rows: one cache line per rank); the suffix array is kept for every
//! `SAMPLE`th text position (a row is located in at most `SAMPLE - 1` LF steps); a table of
//! `PRE`-mer intervals starts every search `PRE` bases in. About 0.6 bytes a text base.
//!
//! Kept next to the index as `<index>.rbwt`, stamped like the repeat index (the index file's size
//! and modification time) and with the parameters it was built with; a mismatch rebuilds it.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::Mutex;

use flexmap::values::{HeaderSeq, VCell};
use flexmap::VD;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

use super::repeat::{as_bytes, as_bytes_mut, each_core, read_reference, stamp_of, PIECE};

const MAGIC: [u8; 8] = *b"FARBWT01";
/// Rows per rank block (one 64-byte cache line: four u32 counts, 384 bits of BWT).
const BB: u64 = 192;
/// Blocks per superblock (whose u64 counts the blocks' u32 counts are relative to).
const SB_SHIFT: u64 = 22;
/// Every `SAMPLE`th text position has its suffix-array row kept.
const SAMPLE: u64 = 32;
/// Searches start from the interval of their last `PRE` bases, looked up.
const PRE: usize = 10;

#[repr(C, align(64))]
#[derive(Clone, Copy, Default)]
struct Block {
    cnt: [u32; 4],
    bits: [u64; 6],
}

/// Each base repeated over a word: `x ^ PAT[c]` has 00 where a 2-bit base is `c`.
const PAT: [u64; 4] = [0, 0x5555_5555_5555_5555, 0xAAAA_AAAA_AAAA_AAAA, u64::MAX];
const LOW: u64 = 0x5555_5555_5555_5555;

#[inline(always)]
fn eq_bits(w: u64, c: usize) -> u64 {
    let x = w ^ PAT[c];
    !(x | (x >> 1)) & LOW
}

#[inline(always)]
fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// The text positions of the sampled rows, in row order: 4 bytes each while the text fits.
enum Samples {
    U32(Vec<u32>),
    U64(Vec<u64>),
}

impl Samples {
    #[inline(always)]
    fn get(&self, i: usize) -> u64 {
        match self {
            Samples::U32(v) => v[i] as u64,
            Samples::U64(v) => v[i],
        }
    }
    fn bytes(&self) -> usize {
        match self {
            Samples::U32(v) => v.len() * 4,
            Samples::U64(v) => v.len() * 8,
        }
    }
}

pub struct RepeatBwt {
    stamp: (u64, u64),
    cap: u64,
    read_len: u32,
    min_dropped: u32,
    /// Text length; the BWT has `n + 1` rows (row 0 is the empty suffix).
    n: u64,
    /// The row whose suffix is the whole text (its BWT symbol is the sentinel, stored as A).
    primary: u64,
    /// Rows starting with a symbol smaller than each base (the sentinel's row counts).
    c: [u64; 5],
    blocks: Vec<Block>,
    sb: Vec<[u64; 4]>,
    /// Rows whose text position is a multiple of `SAMPLE`, and their rank per 512 bits.
    marks: Vec<u64>,
    mark_cum: Vec<u64>,
    samples: Samples,
    /// Per `PRE`-mer (first base most significant), its row interval.
    pre: Vec<[u64; 2]>,
    /// The dropped cores' keys, for telling a dropped core from one absent from the reference.
    keys: Vec<u32>,
    table: Vec<u64>,
    table_shift: u32,
    /// The text's pieces: where each starts in the text (and the end), and where in the reference.
    piece_start: Vec<u64>,
    piece_ref: Vec<u32>,
    piece_pos: Vec<u64>,
}

/// Per-thread buffers of `resolve`.
#[derive(Default, Clone)]
pub struct BwtScratch {
    fwd: Vec<u8>,
    rc: Vec<u8>,
    order: Vec<usize>,
    covered: Vec<bool>,
    cov: Vec<usize>,
}

impl RepeatBwt {
    pub fn min_dropped(&self) -> usize {
        self.min_dropped as usize
    }

    pub fn text_len(&self) -> u64 {
        self.n
    }

    pub fn memory_bytes(&self) -> usize {
        self.blocks.len() * 64 + self.sb.len() * 32 + self.marks.len() * 8 + self.mark_cum.len() * 8
            + self.samples.bytes() + self.pre.len() * 16 + self.keys.len() * 4 + self.table.len() * 8
            + self.piece_start.len() * 8 + self.piece_ref.len() * 4 + self.piece_pos.len() * 8
    }

    pub fn dropped_keys(&self) -> usize {
        self.keys.len()
    }

    /// Whether `key` is a core the main index dropped.
    #[inline]
    pub fn is_dropped(&self, key: u64) -> bool {
        let key = key as u32;
        let mask = self.table.len() - 1;
        let mut slot = home(key, self.table_shift);
        loop {
            let e = self.table[slot];
            if e == 0 {
                return false;
            }
            if (e >> 32) as u32 == key.wrapping_add(1) {
                return true;
            }
            slot = (slot + 1) & mask;
        }
    }

    /// Occurrences of base `c` in BWT rows `[0, i)`.
    #[inline(always)]
    fn rank(&self, c: usize, i: u64) -> u64 {
        let b = (i / BB) as usize;
        let r = (i % BB) as usize;
        let blk = &self.blocks[b];
        let mut n = self.sb[b >> SB_SHIFT][c] + blk.cnt[c] as u64;
        let full = r / 32;
        for w in &blk.bits[..full] {
            n += eq_bits(*w, c).count_ones() as u64;
        }
        let rem = r % 32;
        if rem > 0 {
            n += (eq_bits(blk.bits[full], c) & ((1u64 << (2 * rem)) - 1)).count_ones() as u64;
        }
        // The sentinel is stored as A: the block counts leave it out, the popcount of its own
        // block does not.
        if c == 0 && self.primary < i && self.primary >= b as u64 * BB {
            n -= 1;
        }
        n
    }

    #[inline(always)]
    fn bwt_at(&self, row: u64) -> usize {
        let blk = &self.blocks[(row / BB) as usize];
        let r = (row % BB) as usize;
        ((blk.bits[r / 32] >> (2 * (r % 32))) & 3) as usize
    }

    #[inline(always)]
    fn step(&self, c: usize, lo: u64, hi: u64) -> (u64, u64) {
        (self.c[c] + self.rank(c, lo), self.c[c] + self.rank(c, hi))
    }

    #[inline(always)]
    fn marked(&self, row: u64) -> bool {
        (self.marks[(row / 64) as usize] >> (row % 64)) & 1 == 1
    }

    #[inline(always)]
    fn mark_rank(&self, row: u64) -> u64 {
        let w = (row / 64) as usize;
        let g = w / 8;
        let mut n = self.mark_cum[g];
        for x in &self.marks[g * 8..w] {
            n += x.count_ones() as u64;
        }
        n + (self.marks[w] & ((1u64 << (row % 64)) - 1)).count_ones() as u64
    }

    /// The text position of `row`'s suffix: LF steps back to a sampled row.
    fn locate(&self, row: u64) -> u64 {
        let (mut r, mut steps) = (row, 0u64);
        // The primary row's text position is 0, a multiple of SAMPLE, so it is marked: LF is never
        // asked of it.
        while !self.marked(r) {
            let c = self.bwt_at(r);
            r = self.c[c] + self.rank(c, r);
            steps += 1;
        }
        self.samples.get(self.mark_rank(r) as usize) + steps
    }

    /// The reference place of a match of `len` bases at text position `t`, unless it crosses
    /// from one piece into the next.
    #[inline]
    fn place(&self, t: u64, len: u64) -> Option<(u64, u64)> {
        let j = self.piece_start.partition_point(|&s| s <= t) - 1;
        (t + len <= self.piece_start[j + 1]).then(|| (self.piece_ref[j] as u64, self.piece_pos[j] + (t - self.piece_start[j])))
    }

    /// The match of `codes[j..e)` found by extending left from `e` until it reaches `stop` with at
    /// most `few` occurrences (then also returned as the second value), the read starts, or the
    /// next base does not occur: `(lo, hi, j)` (j == e when even the last `PRE` bases do not occur).
    fn extend(&self, codes: &[u8], e: usize, stop: usize, few: u64) -> ((u64, u64, usize), Option<(u64, u64, usize)>) {
        if e < PRE || codes[e - PRE..e].iter().any(|&c| c > 3) {
            return ((0, 0, e), None);
        }
        let p = codes[e - PRE..e].iter().fold(0usize, |p, &c| p << 2 | c as usize);
        let [mut lo, mut hi] = self.pre[p];
        if lo >= hi {
            return ((0, 0, e), None);
        }
        let mut j = e - PRE;
        let mut first = (hi - lo <= few && j <= stop).then_some((lo, hi, j));
        while j > 0 && first.is_none() {
            let c = codes[j - 1];
            if c > 3 {
                break;
            }
            let (l, h) = self.step(c as usize, lo, hi);
            if l >= h {
                break;
            }
            (lo, hi) = (l, h);
            j -= 1;
            if hi - lo <= few && j <= stop {
                first = Some((lo, hi, j));
            }
        }
        ((lo, hi, j), first)
    }

    /// The reference positions of the K-mers at read positions `qs` (dropped cores' K-mers, the
    /// read as sequenced), as `(i, VD cell)` pairs in `out` for `qs[i]`: the read is searched in
    /// both orientations, from the end of the rightmost K-mer not yet covered, extending the match
    /// left only until it holds that whole K-mer and has at most `max_locate` occurrences -- no
    /// further: a longer match would run across a sequencing error into copies that happen to
    /// share it and lose the true one, and every K-mer it covered with it. Short matches keep an
    /// error to the K-mers around it (as one exact K-mer each would) and leave the copies a base
    /// or two apart as candidates, for the alignment to choose among and the MAPQ to see. A match
    /// covering several K-mers serves them all; a K-mer whose match still has more occurrences
    /// when it reaches the read's start or a difference, or that is not found whole, gets none.
    pub fn resolve(&self, seq: &[u8], k: usize, qs: &[usize], max_locate: usize, s: &mut BwtScratch, out: &mut Vec<(u32, u64)>) {
        let len = seq.len();
        s.fwd.clear();
        s.fwd.extend(seq.iter().map(|&b| code(b)));
        s.rc.clear();
        s.rc.extend(s.fwd.iter().rev().map(|&c| if c > 3 { c } else { 3 - c }));
        for rc in [false, true] {
            let at = |i: usize| if rc { len - qs[i] - k } else { qs[i] };
            s.order.clear();
            s.order.extend(0..qs.len());
            s.order.sort_unstable_by_key(|&i| std::cmp::Reverse(at(i)));
            s.covered.clear();
            s.covered.resize(qs.len(), false);
            let codes: &[u8] = if rc { &s.rc } else { &s.fwd };
            for oi in 0..s.order.len() {
                let i = s.order[oi];
                if s.covered[i] {
                    continue;
                }
                let e = at(i) + k;
                let ((lo, hi, j), first) = self.extend(codes, e, at(i), max_locate as u64);
                if j > at(i) {
                    s.covered[i] = true;
                    continue;
                }
                s.cov.clear();
                s.cov.extend(s.order.iter().copied().filter(|&c| !s.covered[c] && at(c) >= j && at(c) + k <= e));
                let _ = (lo, hi);
                if let Some((flo, fhi, fj)) = first {
                    for row in flo..fhi {
                        if let Some((rid, rpos)) = self.place(self.locate(row), (e - fj) as u64) {
                            for &c in &s.cov {
                                // A K-mer left of where that match began sits at its offset from it.
                                if let Some(p) = (rpos + at(c) as u64).checked_sub(fj as u64) {
                                    out.push((c as u32, VD::set(rid, p)));
                                }
                            }
                        }
                    }
                }
                for &c in &s.cov {
                    s.covered[c] = true;
                }
            }
        }
    }
}

/// A copy of a stretch `[j, e)` of a read, in the orientation searched: the read laid on reference
/// `rid` with its first base (of that orientation) at `start`, reverse-complemented when `rc`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReadCopy {
    pub rid: u64,
    pub start: i64,
    pub rc: bool,
    pub j: u32,
    pub e: u32,
}

impl RepeatBwt {
    /// The match of `codes[j..e)` extended left from `e` as far as it goes: `(lo, hi, j)`, j == e
    /// when not even the last `PRE` bases occur.
    fn extend_full(&self, codes: &[u8], e: usize) -> (u64, u64, usize) {
        if e < PRE || codes[e - PRE..e].iter().any(|&c| c > 3) {
            return (0, 0, e);
        }
        let p = codes[e - PRE..e].iter().fold(0usize, |p, &c| p << 2 | c as usize);
        let [mut lo, mut hi] = self.pre[p];
        if lo >= hi {
            return (0, 0, e);
        }
        let mut j = e - PRE;
        while j > 0 {
            let c = codes[j - 1];
            if c > 3 {
                break;
            }
            let (l, h) = self.step(c as usize, lo, hi);
            if l >= h {
                break;
            }
            (lo, hi) = (l, h);
            j -= 1;
        }
        (lo, hi, j)
    }

    /// Every copy of the read's long exact stretches, in both orientations: the read is cut into
    /// maximal matches from its end leftwards, each restarting left of the base that ended the last
    /// (a sequencing error or a copy difference splits a match instead of steering it into copies
    /// that happen to share the error); each of at least `min_len` bases has all its copies listed
    /// into `out` when they are at most `max`. Returns whether every such match was listed whole --
    /// then `out` holds every place the read's long matches occur.
    pub fn read_copies(&self, seq: &[u8], min_len: usize, max: usize, s: &mut BwtScratch, out: &mut Vec<ReadCopy>) -> bool {
        let len = seq.len();
        s.fwd.clear();
        s.fwd.extend(seq.iter().map(|&b| code(b)));
        s.rc.clear();
        s.rc.extend(s.fwd.iter().rev().map(|&c| if c > 3 { c } else { 3 - c }));
        let mut complete = true;
        for rc in [false, true] {
            let codes: &[u8] = if rc { &s.rc } else { &s.fwd };
            let mut e = len;
            while e >= min_len {
                let (lo, hi, j) = self.extend_full(codes, e);
                if e - j >= min_len {
                    if hi - lo <= max as u64 {
                        for row in lo..hi {
                            if let Some((rid, rpos)) = self.place(self.locate(row), (e - j) as u64) {
                                out.push(ReadCopy { rid, start: rpos as i64 - j as i64, rc, j: j as u32, e: e as u32 });
                            }
                        }
                    } else {
                        complete = false;
                    }
                }
                if j == 0 {
                    break;
                }
                e = if j == e { e - 1 } else { j - 1 };
            }
        }
        out.sort_unstable();
        out.dedup();
        complete
    }
}

fn home(key: u32, shift: u32) -> usize {
    ((key as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> shift) as usize
}

fn make_table(keys: &[u32]) -> (Vec<u64>, u32) {
    let size = (keys.len() * 2).next_power_of_two().max(16);
    let shift = 64 - size.trailing_zeros();
    let mut table = vec![0u64; size];
    for (i, &k) in keys.iter().enumerate() {
        let mut slot = home(k, shift);
        while table[slot] != 0 {
            slot = (slot + 1) & (size - 1);
        }
        table[slot] = ((k as u64 + 1) << 32) | i as u64;
    }
    (table, shift)
}

// --- Per-read storage of the ranges handed to the seed extractor ---

/// Cells and flanks of the current read's BWT ranges: chunks that are never moved while a read is
/// processed, so the slices handed out stay valid until [`reset_scratch`] (called when a read
/// starts, with `refseq2::reset_scratch`).
struct Arena {
    chunks: Vec<Box<[u64]>>,
    chunk: usize,
    used: usize,
}

const ARENA_CHUNK: usize = 1 << 16;

thread_local! {
    static ARENA: RefCell<Arena> = const { RefCell::new(Arena { chunks: Vec::new(), chunk: 0, used: 0 }) };
}

/// `n` cells and `n` flanks equal to `flank`, from this thread's arena, filled by `cells`.
pub fn range_slices<'a>(n: usize, flank: u32, cells: impl Iterator<Item = u64>) -> (&'a [HeaderSeq], &'a [VCell]) {
    ARENA.with(|a| {
        let mut a = a.borrow_mut();
        // n cells (u64) and n flanks (u32, packed two a word).
        let words = n + n.div_ceil(2);
        let p = loop {
            if a.chunk == a.chunks.len() {
                a.chunks.push(vec![0u64; words.max(ARENA_CHUNK)].into_boxed_slice());
                a.used = 0;
            }
            let (chunk, used) = (a.chunk, a.used);
            if used + words <= a.chunks[chunk].len() {
                a.used += words;
                break unsafe { a.chunks[chunk].as_mut_ptr().add(used) };
            }
            a.chunk += 1;
            a.used = 0;
        };
        unsafe {
            let cells_p = p as *mut VCell;
            let mut i = 0;
            for c in cells.take(n) {
                cells_p.add(i).write(VCell(c));
                i += 1;
            }
            debug_assert_eq!(i, n);
            let flanks_p = p.add(n) as *mut u32;
            for j in 0..n {
                flanks_p.add(j).write(flank);
            }
            // HeaderSeq is a transparent u32.
            (std::slice::from_raw_parts(flanks_p as *const HeaderSeq, n), std::slice::from_raw_parts(cells_p as *const VCell, n))
        }
    })
}

/// Forgets this thread's ranges: called when a read (pair) starts.
pub fn reset_scratch() {
    ARENA.with(|a| {
        let mut a = a.borrow_mut();
        a.chunks.truncate(1);
        a.chunk = 0;
        a.used = 0;
    });
}

// --- Build ---

impl RepeatBwt {
    /// Finds the dropped cores (as the main index drops them: over `cap` positions) and their
    /// positions, keeps the regions dense in them, and indexes those.
    pub fn build<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
        reference: &Path,
        rid_of: &HashMap<String, usize>,
        cap: u64,
        read_len: u32,
        min_dropped: u32,
        threads: usize,
    ) -> Self {
        let threads = threads.max(1);
        let t0 = std::time::Instant::now();
        let seqs = read_reference(reference, rid_of);
        let pieces: Vec<(usize, usize, usize)> = seqs
            .iter()
            .enumerate()
            .flat_map(|(r, (_, s))| (0..s.len()).step_by(PIECE).map(move |a| (r, a, (a + PIECE).min(s.len()))))
            .collect();

        // Pass 1: positions per key.
        let size = 1usize << (2 * C);
        let mut counts = vec![0u32; size];
        let atomic: &[AtomicU32] = unsafe { std::slice::from_raw_parts(counts.as_mut_ptr() as *const AtomicU32, size) };
        let next = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    let mut cs = ClosedSyncmer::<C, S, L>::new();
                    while let Some(&(r, a, b)) = pieces.get(next.fetch_add(1, Ordering::Relaxed)) {
                        each_core::<K, C, F, S, L>(&seqs[r].1, a, b, &mut cs, |_, key, _| {
                            atomic[key as usize].fetch_add(1, Ordering::Relaxed);
                        });
                    }
                });
            }
        });
        let keys: Vec<u32> = (0..size as u32).filter(|&k| counts[k as usize] as u64 > cap).collect();

        // Pass 2: the dropped cores' positions, (sequence index << 34 | position).
        let found: Mutex<Vec<u64>> = Mutex::new(Vec::new());
        let next = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    let mut cs = ClosedSyncmer::<C, S, L>::new();
                    let mut local = Vec::new();
                    while let Some(&(r, a, b)) = pieces.get(next.fetch_add(1, Ordering::Relaxed)) {
                        each_core::<K, C, F, S, L>(&seqs[r].1, a, b, &mut cs, |pos, key, _| {
                            if counts[key as usize] as u64 > cap {
                                local.push((r as u64) << 34 | pos as u64);
                            }
                        });
                    }
                    found.lock().unwrap().extend(local);
                });
            }
        });
        drop(counts);
        let mut found = found.into_inner().unwrap();
        found.sort_unstable();
        found.dedup();
        let n_dropped = found.len();

        // The dense positions' windows, per sequence, merged: enough on either side for a read of
        // `read_len` that holds the K-mer.
        let span = (read_len as u64).saturating_sub(K as u64);
        let (mut win_seq, mut win_a, mut win_b): (Vec<u64>, Vec<u64>, Vec<u64>) = (Vec::new(), Vec::new(), Vec::new());
        let mut i = 0;
        while i < found.len() {
            let r = found[i] >> 34;
            let mut j = i;
            while j < found.len() && found[j] >> 34 == r {
                j += 1;
            }
            let ps: Vec<u64> = found[i..j].iter().map(|&x| x & ((1 << 34) - 1)).collect();
            let slen = seqs[r as usize].1.len() as u64;
            let (mut lo, mut hi) = (0usize, 0usize);
            for &p in &ps {
                while ps[lo] + span < p {
                    lo += 1;
                }
                while hi < ps.len() && ps[hi] <= p + span {
                    hi += 1;
                }
                if (hi - lo) < min_dropped as usize {
                    continue;
                }
                let a = p.saturating_sub(span);
                let b = (p + K as u64 + span).min(slen);
                if win_seq.last() == Some(&r) && *win_b.last().unwrap() >= a {
                    let last = win_b.len() - 1;
                    win_b[last] = win_b[last].max(b);
                } else {
                    win_seq.push(r);
                    win_a.push(a);
                    win_b.push(b);
                }
            }
            i = j;
        }
        drop(found);

        // The text: the windows' bases, split where a base is not ACGT; a piece shorter than K
        // cannot hold a match and is left out.
        let mut text: Vec<u8> = Vec::new();
        let (mut piece_start, mut piece_ref, mut piece_pos) = (Vec::new(), Vec::new(), Vec::new());
        for w in 0..win_seq.len() {
            let (rid, s) = (&seqs[win_seq[w] as usize].0, &seqs[win_seq[w] as usize].1);
            let (a, b) = (win_a[w] as usize, win_b[w] as usize);
            let mut x = a;
            while x < b {
                while x < b && code(s[x]) > 3 {
                    x += 1;
                }
                let y = (x..b).find(|&y| code(s[y]) > 3).unwrap_or(b);
                if y - x >= K {
                    piece_start.push(text.len() as u64);
                    piece_ref.push(*rid as u32);
                    piece_pos.push(x as u64);
                    text.extend(s[x..y].iter().map(|&c| code(c)));
                }
                x = y;
            }
        }
        piece_start.push(text.len() as u64);
        drop(seqs);
        let n = text.len() as u64;
        eprintln!(
            "Repeat BWT: {} dropped cores at {} positions; {} windows (>= {} dropped cores within {} bases): {:.3} Gbp of text, {:.1}s",
            keys.len(), n_dropped, win_seq.len(), min_dropped, read_len, n as f64 / 1e9, t0.elapsed().as_secs_f64()
        );

        let mut bwt = Self::from_text(&text, threads);
        let (table, table_shift) = make_table(&keys);
        bwt.cap = cap;
        bwt.read_len = read_len;
        bwt.min_dropped = min_dropped;
        bwt.keys = keys;
        bwt.table = table;
        bwt.table_shift = table_shift;
        bwt.piece_start = piece_start;
        bwt.piece_ref = piece_ref;
        bwt.piece_pos = piece_pos;
        eprintln!("Repeat BWT built in {:.1}s: {:.2} GB", t0.elapsed().as_secs_f64(), bwt.memory_bytes() as f64 / 1e9);
        bwt
    }

    /// The FM-index of `text` (codes 0..3), without the reference bookkeeping.
    fn from_text(text: &[u8], threads: usize) -> Self {
        let n = text.len() as u64;
        // Suffix array, 32-bit while it fits (half the memory of 64).
        let sa_at: Box<dyn Fn(usize) -> u64> = if n < i32::MAX as u64 {
            let sa = libsais::SuffixArrayConstruction::for_text(text)
                .in_owned_buffer32()
                .multi_threaded(libsais::ThreadCount::fixed(threads.min(u16::MAX as usize) as u16))
                .run()
                .expect("suffix array")
                .into_vec();
            Box::new(move |i| sa[i] as u64)
        } else {
            let sa = libsais::SuffixArrayConstruction::for_text(text)
                .in_owned_buffer64()
                .multi_threaded(libsais::ThreadCount::fixed(threads.min(u16::MAX as usize) as u16))
                .run()
                .expect("suffix array")
                .into_vec();
            Box::new(move |i| sa[i] as u64)
        };

        // Row 0 is the empty suffix (text position n), row i + 1 is the suffix array's i.
        let rows = n + 1;
        let nb = (rows / BB + 1) as usize;
        let mut blocks = vec![Block::default(); nb];
        let mut sb = vec![[0u64; 4]; (nb >> SB_SHIFT) + 1];
        let mut marks = vec![0u64; (rows as usize).div_ceil(64) + 1];
        let mut samples: Vec<u64> = Vec::with_capacity((n / SAMPLE + 1) as usize);
        let mut run = [0u64; 4];
        let mut primary = 0;
        for row in 0..nb as u64 * BB {
            let b = (row / BB) as usize;
            let off = (row % BB) as usize;
            if off == 0 {
                if b % (1 << SB_SHIFT) == 0 {
                    sb[b >> SB_SHIFT] = run;
                }
                let base = sb[b >> SB_SHIFT];
                blocks[b].cnt = [0, 1, 2, 3].map(|c| (run[c] - base[c]) as u32);
            }
            if row >= rows {
                continue;
            }
            let t = if row == 0 { n } else { sa_at(row as usize - 1) };
            let c = if t == 0 {
                primary = row;
                0
            } else {
                text[t as usize - 1] as usize
            };
            blocks[b].bits[off / 32] |= (c as u64) << (2 * (off % 32));
            if t != 0 {
                run[c] += 1;
            }
            if t % SAMPLE == 0 {
                marks[(row / 64) as usize] |= 1 << (row % 64);
                samples.push(t);
            }
        }
        drop(sa_at);
        let mut mark_cum = Vec::with_capacity(marks.len() / 8 + 1);
        let mut acc = 0u64;
        for g in marks.chunks(8) {
            mark_cum.push(acc);
            acc += g.iter().map(|x| x.count_ones() as u64).sum::<u64>();
        }
        mark_cum.push(acc);
        let mut c = [0u64; 5];
        c[0] = 1;
        for i in 0..4 {
            c[i + 1] = c[i] + run[i];
        }
        let samples = if n < u32::MAX as u64 { Samples::U32(samples.into_iter().map(|x| x as u32).collect()) } else { Samples::U64(samples) };
        let mut bwt = Self {
            stamp: (0, 0), cap: 0, read_len: 0, min_dropped: 0, n, primary, c, blocks, sb, marks, mark_cum, samples,
            pre: Vec::new(), keys: Vec::new(), table: vec![0; 16], table_shift: 60,
            piece_start: vec![0, n], piece_ref: vec![0], piece_pos: vec![0],
        };
        bwt.pre = (0..1usize << (2 * PRE))
            .map(|p| {
                let (mut lo, mut hi) = (0u64, n + 1);
                for k in (0..PRE).rev() {
                    if lo >= hi {
                        break;
                    }
                    (lo, hi) = bwt.step((p >> (2 * (PRE - 1 - k))) & 3, lo, hi);
                }
                [lo, hi.max(lo)]
            })
            .collect();
        bwt
    }

    // --- Files ---

    pub fn path_for(index_file: &str) -> String {
        format!("{index_file}.rbwt")
    }

    /// Loads `<index_file>.rbwt` if it was built from this index with these parameters; otherwise
    /// builds and writes it. Returns it and what was done.
    pub fn load_or_build<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
        index_file: &str,
        reference: &Path,
        rid_of: impl FnOnce() -> HashMap<String, usize>,
        cap: u64,
        read_len: u32,
        min_dropped: u32,
        threads: usize,
    ) -> (Self, &'static str) {
        let stamp = stamp_of(index_file).unwrap_or((0, 0));
        let path = Self::path_for(index_file);
        if let Ok(b) = Self::load(&path) {
            if b.stamp == stamp && b.cap == cap && b.read_len == read_len && b.min_dropped == min_dropped {
                return (b, "read");
            }
        }
        let mut b = Self::build::<K, C, F, S, L>(reference, &rid_of(), cap, read_len, min_dropped, threads);
        b.stamp = stamp;
        if let Err(e) = b.save(&path) {
            log::warn!("could not write the repeat BWT {path}: {e}");
        }
        (b, "built")
    }

    fn save(&self, path: &str) -> std::io::Result<()> {
        let tmp = format!("{path}.tmp{}", std::process::id());
        {
            let mut w = BufWriter::with_capacity(1 << 24, File::create(&tmp)?);
            w.write_all(&MAGIC)?;
            let (s32, s64): (&[u32], &[u64]) = match &self.samples {
                Samples::U32(v) => (v, &[]),
                Samples::U64(v) => (&[], v),
            };
            let head = [
                self.stamp.0, self.stamp.1, self.cap, self.read_len as u64, self.min_dropped as u64, self.n, self.primary,
                self.c[0], self.c[1], self.c[2], self.c[3], self.c[4],
                self.blocks.len() as u64, self.sb.len() as u64, self.marks.len() as u64, self.mark_cum.len() as u64,
                s32.len() as u64, s64.len() as u64, self.pre.len() as u64, self.keys.len() as u64,
                self.table.len() as u64, self.table_shift as u64, self.piece_start.len() as u64,
            ];
            w.write_all(as_bytes(&head))?;
            w.write_all(as_bytes(&self.blocks))?;
            w.write_all(as_bytes(&self.sb))?;
            w.write_all(as_bytes(&self.marks))?;
            w.write_all(as_bytes(&self.mark_cum))?;
            w.write_all(as_bytes(s32))?;
            w.write_all(as_bytes(s64))?;
            w.write_all(as_bytes(&self.pre))?;
            w.write_all(as_bytes(&self.keys))?;
            w.write_all(as_bytes(&self.table))?;
            w.write_all(as_bytes(&self.piece_start))?;
            w.write_all(as_bytes(&self.piece_ref))?;
            w.write_all(as_bytes(&self.piece_pos))?;
            w.flush()?;
        }
        std::fs::rename(&tmp, path)
    }

    fn load(path: &str) -> std::io::Result<Self> {
        let mut r = BufReader::with_capacity(1 << 24, File::open(path)?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "not a repeat BWT"));
        }
        let mut head = [0u64; 23];
        r.read_exact(as_bytes_mut(&mut head))?;
        fn read<T: Copy + Default>(r: &mut impl Read, n: u64) -> std::io::Result<Vec<T>> {
            let mut v = vec![T::default(); n as usize];
            r.read_exact(as_bytes_mut(&mut v))?;
            Ok(v)
        }
        let blocks: Vec<Block> = read(&mut r, head[12])?;
        let sb: Vec<[u64; 4]> = read(&mut r, head[13])?;
        let marks: Vec<u64> = read(&mut r, head[14])?;
        let mark_cum: Vec<u64> = read(&mut r, head[15])?;
        let s32: Vec<u32> = read(&mut r, head[16])?;
        let s64: Vec<u64> = read(&mut r, head[17])?;
        let pre: Vec<[u64; 2]> = read(&mut r, head[18])?;
        let keys: Vec<u32> = read(&mut r, head[19])?;
        let table: Vec<u64> = read(&mut r, head[20])?;
        let piece_start: Vec<u64> = read(&mut r, head[22])?;
        let piece_ref: Vec<u32> = read(&mut r, head[22] - 1)?;
        let piece_pos: Vec<u64> = read(&mut r, head[22] - 1)?;
        Ok(Self {
            stamp: (head[0], head[1]), cap: head[2], read_len: head[3] as u32, min_dropped: head[4] as u32,
            n: head[5], primary: head[6], c: [head[7], head[8], head[9], head[10], head[11]],
            blocks, sb, marks, mark_cum,
            samples: if head[17] > 0 { Samples::U64(s64) } else { Samples::U32(s32) },
            pre, keys, table, table_shift: head[21] as u32, piece_start, piece_ref, piece_pos,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_of(s: &[u8]) -> Vec<u8> {
        s.iter().map(|&b| code(b)).collect()
    }

    /// Every occurrence of `pat` in `text`, the slow way.
    fn naive(text: &[u8], pat: &[u8]) -> Vec<u64> {
        (0..=text.len().saturating_sub(pat.len())).filter(|&i| &text[i..i + pat.len()] == pat).map(|i| i as u64).collect()
    }

    #[test]
    fn counts_and_locates_every_pattern() {
        // A repeat-rich text: copies of a unit with a few changes, past a block and a sample.
        let mut x = 7u64;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let unit: Vec<u8> = (0..300).map(|_| b"ACGT"[(rnd() % 4) as usize]).collect();
        let mut s = Vec::new();
        for _ in 0..40 {
            let mut u = unit.clone();
            let i = (rnd() % 300) as usize;
            u[i] = b"ACGT"[(rnd() % 4) as usize];
            s.extend(u);
            s.extend((0..(rnd() % 50)).map(|_| b"ACGT"[(rnd() % 4) as usize]));
        }
        let text = text_of(&s);
        let bwt = RepeatBwt::from_text(&text, 2);
        for start in [0usize, 5, 299, 1000, 4000, s.len() - 40] {
            for len in [1usize, 10, 11, 12, 31, 40] {
                if start + len > s.len() {
                    continue;
                }
                let pat = &text[start..start + len];
                let (mut lo, mut hi) = (0u64, bwt.n + 1);
                for &c in pat.iter().rev() {
                    (lo, hi) = bwt.step(c as usize, lo, hi);
                }
                let want = naive(&text, pat);
                assert_eq!(hi - lo, want.len() as u64, "count of {start}+{len}");
                let mut got: Vec<u64> = (lo..hi).map(|r| bwt.locate(r)).collect();
                got.sort_unstable();
                assert_eq!(got, want, "positions of {start}+{len}");
                if len == PRE {
                    let p = pat.iter().fold(0usize, |p, &c| p << 2 | c as usize);
                    assert_eq!(bwt.pre[p], [lo, hi]);
                }
            }
        }
    }

    #[test]
    fn lists_the_copies_of_a_read_across_an_error() {
        // Five copies of a 200-base unit, copy 2 differing in one base; a read from copy 2 with a
        // sequencing error: its stretches either side of the error find copy 2 (and the copies
        // sharing them), in both orientations, and the list is complete.
        let mut x = 5u64;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let unit: Vec<u8> = (0..200).map(|_| b"ACGT"[(rnd() % 4) as usize]).collect();
        let mut s = Vec::new();
        let mut starts = Vec::new();
        for copy in 0..5 {
            let mut u = unit.clone();
            if copy == 2 {
                u[100] = if u[100] == b'A' { b'C' } else { b'A' };
            }
            starts.push(s.len());
            s.extend(u);
            s.extend((0..41).map(|_| b"ACGT"[(rnd() % 4) as usize]));
        }
        let text = text_of(&s);
        let mut bwt = RepeatBwt::from_text(&text, 2);
        bwt.piece_ref = vec![7];
        bwt.piece_pos = vec![0];
        let mut read = s[starts[2] + 20..starts[2] + 170].to_vec();
        read[30] = if read[30] == b'G' { b'T' } else { b'G' };
        let mut scratch = BwtScratch::default();
        for rc in [false, true] {
            let r: Vec<u8> = if rc {
                read.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }).collect()
            } else {
                read.clone()
            };
            let mut out = Vec::new();
            assert!(bwt.read_copies(&r, 40, 16, &mut scratch, &mut out));
            assert!(out.iter().any(|c| c.rid == 7 && c.start == (starts[2] + 20) as i64 && c.rc == rc), "rc {rc}: {out:?}");
            // The stretch holding copy 2's own base (read position 80) finds copy 2 alone.
            let own: Vec<_> = out.iter().filter(|c| c.rc == rc && {
                let (a, b) = if rc { (150 - c.e as i64, 150 - c.j as i64) } else { (c.j as i64, c.e as i64) };
                a <= 80 && 80 < b
            }).collect();
            assert!(!own.is_empty() && own.iter().all(|c| c.start == (starts[2] + 20) as i64), "rc {rc}: {own:?}");
        }
    }

    #[test]
    fn resolves_a_read_from_one_copy() {
        // Ten copies of a 200-base repeat, each differing in one base; a read from copy 6, in both
        // orientations.
        let mut x = 11u64;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let unit: Vec<u8> = (0..200).map(|_| b"ACGT"[(rnd() % 4) as usize]).collect();
        let mut s = Vec::new();
        let mut starts = Vec::new();
        for copy in 0..10 {
            let mut u = unit.clone();
            u[20 * copy + 5] = if u[20 * copy + 5] == b'A' { b'C' } else { b'A' };
            starts.push(s.len());
            s.extend(u);
            s.extend((0..37).map(|_| b"ACGT"[(rnd() % 4) as usize]));
        }
        let text = text_of(&s);
        let mut bwt = RepeatBwt::from_text(&text, 2);
        bwt.piece_ref = vec![3];
        bwt.piece_pos = vec![1000];
        let read = &s[starts[6] + 30..starts[6] + 180];
        let qs = [0usize, 60, 119];
        let mut scratch = BwtScratch::default();
        for rc in [false, true] {
            let r: Vec<u8> = if rc {
                read.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }).collect()
            } else {
                read.to_vec()
            };
            let mut out = Vec::new();
            bwt.resolve(&r, 31, &qs, 4, &mut scratch, &mut out);
            out.sort_unstable();
            let want: Vec<(u32, u64)> = qs
                .iter()
                .enumerate()
                .map(|(i, &q)| {
                    let q = if rc { 150 - q - 31 } else { q };
                    (i as u32, VD::set(3, 1000 + (starts[6] + 30 + q) as u64))
                })
                .collect();
            // A K-mer that gets positions has the true copy among them. The K-mer nearest the read's
            // end reaches copy 6's own base, so it is resolved; those left of it extend only left,
            // where 7 copies agree (> max_locate), and stay unresolved.
            for (i, w) in want.iter().enumerate() {
                if out.iter().any(|h| h.0 == i as u32) {
                    assert!(out.contains(w), "rc {rc}: {w:?} missing from {out:?}");
                }
            }
            let last = if rc { 0 } else { 2 };
            assert!(out.contains(&want[last]), "rc {rc}: the K-mer at the read's end is not resolved: {out:?}");
            for i in 0..qs.len() as u32 {
                assert!(out.iter().filter(|h| h.0 == i).count() <= 4, "rc {rc}: core {i} over max_locate");
            }
        }
    }
}
