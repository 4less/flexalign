//! A whole-genome FM-index for reads whose cores the main index cannot resolve (`--fm`).
//!
//! The main index answers "where does this 31-mer occur" for cores of at most --max-range-size
//! copies, and seeding drops cores whose copies tie on the flank. A read from a young repeat then
//! has no seed that tells its copy from the others, though the whole read does. This index answers
//! "where does this stretch of the read occur" for a match of any length, so the copies sharing a
//! long stretch of the read (few) can be listed without listing those sharing 31 bases (many).
//!
//! **Pieces.** The references are cut, in id order and at reference boundaries, into pieces of at
//! most `max_piece` bases (default just under 2^32), each its own FM-index. A read pair lies on one
//! reference, so every candidate a pair needs is found within one piece, and a search over all
//! pieces finds every candidate a single index would (plus, where a piece's longest match is
//! shorter than the global one, a few more). Pieces keep positions in 32 bits, keep the build at
//! one piece's suffix array (wheat and the metagenomics DB do not fit one), and are the unit a
//! RAM-bounded (sharded) run can load one at a time. Human fits one piece.
//!
//! **Layout** (1 byte a base). Per 128 BWT rows one 64-byte cache line: the counts of A, C, G
//! before the line (T follows from the row), the rank of the sampled rows before it, the 128 bases
//! at 2 bits, and a bit per row: is its text position a multiple of SAMPLE. One backward step and
//! one locate step read one line. The text positions of the sampled rows, in row order, as u32.
//! Locating a row walks the BWT back to a sampled text position: at most SAMPLE - 1 steps.
//!
//! Text: a piece's references concatenated at 2 bits a base (a base outside ACGT is A: a match
//! running through it is checked against the reference afterwards anyway), plus the sentinel. A
//! match may run across two references; the caller checks. Kept next to the main index as
//! `<reference>.fmd`.

use std::io::{Read, Write};

const BLK: u64 = 128;
const SAMPLE: u64 = 8;
const MAGIC: [u8; 8] = *b"FLXFMD02";
/// The largest piece: rows (bases + sentinel) must fit a u32.
pub const MAX_PIECE: u64 = u32::MAX as u64 - 1;

#[repr(C, align(64))]
#[derive(Clone, Copy, Default)]
struct Line {
    /// Occurrences of A, C, G in the rows before this line (the sentinel not counted).
    cnt: [u32; 3],
    /// Sampled rows before this line.
    smp_rank: u32,
    bwt: [u64; 4],
    smp: [u64; 2],
}
const _: () = assert!(std::mem::size_of::<Line>() == 64);

#[inline]
fn match_mask(w: u64, a: u8) -> u64 {
    const PAT: [u64; 4] = [0, 0x5555_5555_5555_5555, 0xAAAA_AAAA_AAAA_AAAA, 0xFFFF_FFFF_FFFF_FFFF];
    let x = w ^ PAT[a as usize];
    !(x | (x >> 1)) & 0x5555_5555_5555_5555
}

/// One piece.
pub struct Fmd {
    /// Rows: the piece's text length plus the sentinel.
    pub n: u64,
    lines: Vec<Line>,
    c: [u64; 4],
    /// The row whose BWT entry is the sentinel (stored as A, corrected in `rank`).
    primary: u64,
    sa: Vec<u32>,
    /// Text offset of each of the piece's references; one past the end last.
    starts: Vec<u64>,
    /// The global id of the piece's first reference.
    first_id: usize,
}

impl Fmd {
    /// Occurrences of base `a` in BWT rows [0, i) of the line holding row `i` (`i` < n), or at
    /// the end (`i` == n).
    #[inline]
    pub fn rank(&self, a: u8, i: u64) -> u64 {
        let li = i / BLK;
        let l = &self.lines[li as usize];
        let start = li * BLK;
        let mut r = if a < 3 {
            l.cnt[a as usize] as u64
        } else {
            start - (l.cnt[0] as u64 + l.cnt[1] as u64 + l.cnt[2] as u64) - (self.primary < start) as u64
        };
        let off = i - start;
        let full = (off / 32) as usize;
        for w in 0..full {
            r += match_mask(l.bwt[w], a).count_ones() as u64;
        }
        let rem = off % 32;
        if rem > 0 {
            r += (match_mask(l.bwt[full], a) & ((1u64 << (2 * rem)) - 1)).count_ones() as u64;
        }
        if a == 0 && self.primary < i && self.primary >= start {
            r -= 1;
        }
        r
    }

    /// One backward step: the rows of `a` followed by the rows [lo, hi).
    #[inline]
    pub fn step(&self, a: u8, lo: u64, hi: u64) -> (u64, u64) {
        (self.c[a as usize] + self.rank(a, lo), self.c[a as usize] + self.rank(a, hi))
    }

    /// Prefetches the lines `step` will read for rows lo and hi.
    #[inline]
    pub fn prefetch(&self, lo: u64, hi: u64) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            for i in [lo, hi] {
                _mm_prefetch(self.lines.as_ptr().add((i / BLK) as usize) as *const i8, _MM_HINT_T0);
            }
        }
    }

    /// The text position of the suffix at `row`, and the BWT steps it took.
    #[inline]
    pub fn locate(&self, mut row: u64) -> (u64, u32) {
        let mut steps = 0u32;
        loop {
            let li = row / BLK;
            let l = &self.lines[li as usize];
            let off = row - li * BLK;
            let (w, b) = ((off / 64) as usize, off % 64);
            if (l.smp[w] >> b) & 1 == 1 {
                let below = if w == 1 { l.smp[0].count_ones() as u64 } else { 0 } + (l.smp[w] & ((1u64 << b) - 1)).count_ones() as u64;
                return (self.sa[(l.smp_rank as u64 + below) as usize] as u64 + steps as u64, steps);
            }
            if row == self.primary {
                return (steps as u64, steps);
            }
            let a = ((l.bwt[(off / 32) as usize] >> (2 * (off % 32))) & 3) as u8;
            row = self.c[a as usize] + self.rank(a, row);
            steps += 1;
        }
    }

    /// (global reference id, offset in it) of a text position of this piece.
    #[inline]
    pub fn to_ref(&self, pos: u64) -> (usize, u64) {
        let id = self.starts.partition_point(|&s| s <= pos) - 1;
        (self.first_id + id, pos - self.starts[id])
    }

    /// The length of a reference of this piece, by global id.
    #[inline]
    pub fn ref_len(&self, id: usize) -> u64 {
        let i = id - self.first_id;
        self.starts[i + 1] - self.starts[i]
    }

    fn build_piece(text: &[u8], starts: Vec<u64>, first_id: usize) -> Fmd {
        let len = text.len() as u64;
        assert!(len + 1 <= MAX_PIECE + 1, "an FM piece is limited to {MAX_PIECE} bases");
        let sa: Vec<i64> = libsais::SuffixArrayConstruction::for_text(text)
            .in_owned_buffer()
            .multi_threaded(libsais::ThreadCount::openmp_default())
            .run()
            .expect("suffix array")
            .into_vec();
        let n = len + 1;
        // Row 0 is the sentinel's suffix; row r + 1 is sa[r].
        let suffix = |row: u64| -> u64 { if row == 0 { len } else { sa[(row - 1) as usize] as u64 } };
        let mut lines = vec![Line::default(); (n / BLK + 1) as usize];
        let mut samples: Vec<u32> = Vec::with_capacity((n / SAMPLE + 2) as usize);
        let mut cnt = [0u64; 4];
        let mut primary = 0u64;
        for row in 0..n {
            let off = row % BLK;
            let l = &mut lines[(row / BLK) as usize];
            if off == 0 {
                l.cnt = [cnt[0] as u32, cnt[1] as u32, cnt[2] as u32];
                l.smp_rank = samples.len() as u32;
            }
            let p = suffix(row);
            if p % SAMPLE == 0 {
                l.smp[(off / 64) as usize] |= 1 << (off % 64);
                samples.push(p as u32);
            }
            if p == 0 {
                primary = row;
            } else {
                let a = text[(p - 1) as usize];
                l.bwt[(off / 32) as usize] |= (a as u64) << (2 * (off % 32));
                cnt[a as usize] += 1;
            }
        }
        if n % BLK == 0 {
            // The line past the end, read by rank(a, n).
            let l = &mut lines[(n / BLK) as usize];
            l.cnt = [cnt[0] as u32, cnt[1] as u32, cnt[2] as u32];
            l.smp_rank = samples.len() as u32;
        }
        let mut c = [0u64; 4];
        c[0] = 1;
        for a in 1..4 {
            c[a] = c[a - 1] + cnt[a - 1];
        }
        Fmd { n, lines, c, primary, sa: samples, starts, first_id }
    }

    fn write(&self, w: &mut dyn Write) -> std::io::Result<()> {
        for v in [self.n, self.primary, self.c[0], self.c[1], self.c[2], self.c[3], self.first_id as u64, self.lines.len() as u64, self.sa.len() as u64, self.starts.len() as u64] {
            w.write_all(&v.to_le_bytes())?;
        }
        // Lines and samples as their in-memory bytes (little-endian hosts only).
        let lb = unsafe { std::slice::from_raw_parts(self.lines.as_ptr() as *const u8, self.lines.len() * 64) };
        for ch in lb.chunks(1 << 26) {
            w.write_all(ch)?;
        }
        let sb = unsafe { std::slice::from_raw_parts(self.sa.as_ptr() as *const u8, self.sa.len() * 4) };
        for ch in sb.chunks(1 << 26) {
            w.write_all(ch)?;
        }
        for s in &self.starts {
            w.write_all(&s.to_le_bytes())?;
        }
        Ok(())
    }

    fn read(f: &mut dyn Read) -> std::io::Result<Fmd> {
        let mut h = [0u64; 10];
        for x in h.iter_mut() {
            let mut b = [0u8; 8];
            f.read_exact(&mut b)?;
            *x = u64::from_le_bytes(b);
        }
        let mut lines = vec![Line::default(); h[7] as usize];
        f.read_exact(unsafe { std::slice::from_raw_parts_mut(lines.as_mut_ptr() as *mut u8, lines.len() * 64) })?;
        let mut sa = vec![0u32; h[8] as usize];
        f.read_exact(unsafe { std::slice::from_raw_parts_mut(sa.as_mut_ptr() as *mut u8, sa.len() * 4) })?;
        let mut starts = vec![0u64; h[9] as usize];
        for s in starts.iter_mut() {
            let mut b = [0u8; 8];
            f.read_exact(&mut b)?;
            *s = u64::from_le_bytes(b);
        }
        Ok(Fmd { n: h[0], primary: h[1], c: [h[2], h[3], h[4], h[5]], first_id: h[6] as usize, lines, sa, starts })
    }
}

/// The index: its pieces in reference-id order.
pub struct FmdSet {
    pub pieces: Vec<Fmd>,
}

impl FmdSet {
    /// Builds the index of `refs` (each reference's bases, in id order), cut into pieces of at
    /// most `max_piece` bases at reference boundaries. A reference longer than a piece is an error.
    pub fn build(refs: &mut dyn Iterator<Item = Vec<u8>>, max_piece: u64) -> FmdSet {
        let max_piece = max_piece.min(MAX_PIECE);
        let mut pieces = Vec::new();
        let mut text: Vec<u8> = Vec::new();
        let mut starts = vec![0u64];
        let mut first_id = 0usize;
        let mut id = 0usize;
        let flush = |text: &mut Vec<u8>, starts: &mut Vec<u64>, first_id: usize, pieces: &mut Vec<Fmd>| {
            eprintln!("[fm] piece {}: {} references, {:.2} Gbp, building", pieces.len(), starts.len() - 1, text.len() as f64 / 1e9);
            pieces.push(Fmd::build_piece(text, std::mem::replace(starts, vec![0]), first_id));
            text.clear();
        };
        for r in refs {
            assert!(r.len() as u64 <= max_piece, "reference {id} ({} bp) is longer than an FM piece ({max_piece})", r.len());
            if !text.is_empty() && (text.len() + r.len()) as u64 > max_piece {
                flush(&mut text, &mut starts, first_id, &mut pieces);
                first_id = id;
            }
            text.extend(r.iter().map(|&b| match b { b'C' | b'c' => 1u8, b'G' | b'g' => 2, b'T' | b't' => 3, _ => 0 }));
            starts.push(text.len() as u64);
            id += 1;
        }
        if starts.len() > 1 {
            flush(&mut text, &mut starts, first_id, &mut pieces);
        }
        FmdSet { pieces }
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let tmp = format!("{path}.tmp");
        let mut w = std::io::BufWriter::with_capacity(1 << 24, std::fs::File::create(&tmp)?);
        w.write_all(&MAGIC)?;
        w.write_all(&(self.pieces.len() as u64).to_le_bytes())?;
        for p in &self.pieces {
            p.write(&mut w)?;
        }
        w.flush()?;
        drop(w);
        std::fs::rename(tmp, path)
    }

    pub fn load(path: &str) -> std::io::Result<FmdSet> {
        let mut f = std::io::BufReader::with_capacity(1 << 24, std::fs::File::open(path)?);
        let mut magic = [0u8; 8];
        f.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{path} is not an FM-index of this version (rebuild it with --fm-build)")));
        }
        let mut b = [0u8; 8];
        f.read_exact(&mut b)?;
        let np = u64::from_le_bytes(b);
        let pieces = (0..np).map(|_| Fmd::read(&mut f)).collect::<std::io::Result<Vec<_>>>()?;
        Ok(FmdSet { pieces })
    }
}

/// The loaded index of the run (`--fm`).
pub static FMD: std::sync::OnceLock<FmdSet> = std::sync::OnceLock::new();

#[cfg(test)]
mod tests {
    use super::*;

    fn code(b: u8) -> u8 {
        match b { b'C' => 1, b'G' => 2, b'T' => 3, _ => 0 }
    }

    /// Every (reference, offset) of `q` lying within one reference, by the index.
    fn find(set: &FmdSet, q: &[u8]) -> Vec<(usize, u64)> {
        let mut got = Vec::new();
        for fm in &set.pieces {
            let (mut lo, mut hi) = (0, fm.n);
            for &b in q.iter().rev() {
                (lo, hi) = fm.step(code(b), lo, hi);
            }
            for r in lo..hi {
                let (id, off) = fm.to_ref(fm.locate(r).0);
                if off + q.len() as u64 <= fm.ref_len(id) {
                    got.push((id, off));
                }
            }
        }
        got.sort_unstable();
        got
    }

    fn naive(refs: &[Vec<u8>], q: &[u8]) -> Vec<(usize, u64)> {
        let mut v = Vec::new();
        for (id, r) in refs.iter().enumerate() {
            if r.len() >= q.len() {
                v.extend((0..=r.len() - q.len()).filter(|&i| &r[i..i + q.len()] == q).map(|i| (id, i as u64)));
            }
        }
        v
    }

    #[test]
    fn finds_every_occurrence() {
        let mut x: u64 = 7;
        let mut refs = Vec::new();
        for len in [5000usize, 1, 3333, 777, 128, 2000] {
            let mut r = Vec::with_capacity(len);
            for _ in 0..len {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                r.push(b"ACGT"[(x >> 62) as usize]);
            }
            refs.push(r);
        }
        // Repeats: a stretch of the first reference in the third and the last.
        let rep = refs[0][100..400].to_vec();
        refs[2][1000..1300].copy_from_slice(&rep);
        refs[5][10..310].copy_from_slice(&rep);
        for max_piece in [MAX_PIECE, 5001, 6000] {
            let set = FmdSet::build(&mut refs.clone().into_iter(), max_piece);
            for (id, a, l) in [(0usize, 0usize, 12usize), (0, 150, 31), (2, 1200, 9), (0, 4000, 20), (2, 100, 40), (0, 100, 300), (5, 1990, 10), (3, 0, 5)] {
                let q = &refs[id][a..a + l];
                assert_eq!(find(&set, q), naive(&refs, q), "max_piece {max_piece}, query {id}:{a} len {l}");
            }
            let path = std::env::temp_dir().join(format!("fmd-test-{}-{max_piece}", std::process::id()));
            set.save(path.to_str().unwrap()).unwrap();
            let set2 = FmdSet::load(path.to_str().unwrap()).unwrap();
            std::fs::remove_file(path).unwrap();
            assert_eq!(set2.pieces.len(), set.pieces.len());
            assert_eq!(find(&set2, &rep), find(&set, &rep));
        }
        // Pieces of 5001 bases: the first reference alone, then 1, 3333, 777 (4111), then 128, 2000.
        assert_eq!(FmdSet::build(&mut refs.clone().into_iter(), 5001).pieces.len(), 3);
    }
}
