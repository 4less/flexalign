//! The flank side store (`--flankext-build`, `--flankext`): for the cores of the extended
//! multiplicity brackets (`--flank-brackets` / `--flank-lengths`), every copy's reference window
//! around its k-mer, packed at 2 bits a base, so `--flank-ext` compares a cut core's tied copies
//! with the read without touching the reference.
//!
//! Kept apart from the main index: a core outside the extended brackets never looks here, and the
//! main value blocks are unchanged. Cells follow the main range's copy order, so the index of a
//! copy in `flank_distances` addresses its cell directly.
//!
//! Per copy: the window [p - e, p + K + e) of the forward reference strand (p the copy's k-mer
//! start, K = 31, e = (flank length - 16) / 2 of the core's bracket), at 2 bits a base, base t at
//! bits 2 (t % 32) of word t / 32. A window that runs past its reference or holds a base outside
//! ACGT is an exception (a byte flag): the caller compares it against the reference as before.
//!
//! Key table: open addressing on the index key (already a mixed value), slots of (key, first
//! word, first flag); empty = u64::MAX.

use std::io::{Read, Write};

const MAGIC: [u8; 8] = *b"FLXEXT01";
pub const KMER: usize = 31;

pub struct FlankExt {
    pub brackets: Vec<usize>,
    pub lengths: Vec<usize>,
    /// Size of the main index file the store was built for (0: unknown).
    pub index_bytes: u64,
    table: Vec<[u64; 3]>,
    words: Vec<u64>,
    exc: Vec<u8>,
}

/// Extra bases a side for a core of `n` positions, by bracket.
pub fn ext_for(brackets: &[usize], lengths: &[usize], n: usize) -> usize {
    if lengths.is_empty() {
        return 0;
    }
    let bi = brackets.iter().position(|&b| n <= b).unwrap_or(brackets.len());
    let l = *lengths.get(bi).or(lengths.last()).unwrap_or(&16);
    l.saturating_sub(16) / 2
}

/// Words a packed window of `e` extra bases a side takes.
#[inline]
pub fn words_per_copy(e: usize) -> usize {
    (2 * (KMER + 2 * e)).div_ceil(64)
}

#[inline]
fn code(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Packs `w` (ACGT only) into `out`; false when a base is outside ACGT.
pub fn pack_window(w: &[u8], out: &mut [u64]) -> bool {
    out.iter_mut().for_each(|x| *x = 0);
    for (t, &b) in w.iter().enumerate() {
        let Some(c) = code(b) else { return false };
        out[t / 32] |= c << (2 * (t % 32));
    }
    true
}

/// The read side of a comparison: the read's bases at window offsets t (read position
/// kpos - e + t), packed, with a mask of the offsets inside the read and of those that are N.
pub struct ReadWindow {
    pub bases: Vec<u64>,
    pub valid: Vec<u64>,
    pub n: Vec<u64>,
}

impl ReadWindow {
    pub fn new(read: &[u8], kpos: i64, e: usize) -> ReadWindow {
        let w = KMER + 2 * e;
        let nw = words_per_copy(e);
        let (mut bases, mut valid, mut n) = (vec![0u64; nw], vec![0u64; nw], vec![0u64; nw]);
        for t in 0..w {
            let rp = kpos - e as i64 + t as i64;
            if rp < 0 || rp >= read.len() as i64 {
                continue;
            }
            let bit = 1u64 << (2 * (t % 32));
            valid[t / 32] |= bit;
            match code(read[rp as usize]) {
                Some(c) => bases[t / 32] |= c << (2 * (t % 32)),
                None => n[t / 32] |= bit,
            }
        }
        ReadWindow { bases, valid, n }
    }

    /// Mismatches against a packed reference window (ACGT only): a read N always mismatches.
    #[inline]
    pub fn mismatches(&self, cell: &[u64]) -> usize {
        let mut m = 0u32;
        for i in 0..cell.len() {
            let x = self.bases[i] ^ cell[i];
            let d = (x | (x >> 1)) & 0x5555_5555_5555_5555;
            m += (d & self.valid[i] & !self.n[i]).count_ones() + (self.n[i] & self.valid[i]).count_ones();
        }
        m as usize
    }
}

impl FlankExt {
    /// (first word, first flag) of a core's cells.
    #[inline]
    pub fn lookup(&self, key: u64) -> Option<(usize, usize)> {
        if self.table.is_empty() {
            return None;
        }
        let mask = self.table.len() - 1;
        let mut s = (key as usize) & mask;
        loop {
            let e = &self.table[s];
            if e[0] == key {
                return Some((e[1] as usize, e[2] as usize));
            }
            if e[0] == u64::MAX {
                return None;
            }
            s = (s + 1) & mask;
        }
    }

    #[inline]
    pub fn cell(&self, word: usize, nw: usize) -> &[u64] {
        &self.words[word..word + nw]
    }

    #[inline]
    pub fn is_exception(&self, flag: usize) -> bool {
        self.exc[flag] != 0
    }

    /// Builds the store from per-core cells: `cores` lists (key, copies, e) in fill order and
    /// `fill(core index, words, flags)` writes a core's cells (n * words_per_copy(e) words, n
    /// flags), called from `threads` threads on disjoint cores.
    pub fn assemble(
        brackets: Vec<usize>, lengths: Vec<usize>, index_bytes: u64, cores: &[(u64, usize, usize)], threads: usize,
        fill: &(dyn Fn(usize, &mut [u64], &mut [u8]) + Sync),
    ) -> FlankExt {
        let mut word_off = Vec::with_capacity(cores.len() + 1);
        let mut flag_off = Vec::with_capacity(cores.len() + 1);
        let (mut w, mut f) = (0usize, 0usize);
        for &(_, n, e) in cores {
            word_off.push(w);
            flag_off.push(f);
            w += n * words_per_copy(e);
            f += n;
        }
        word_off.push(w);
        flag_off.push(f);
        let mut words = vec![0u64; w];
        let mut exc = vec![0u8; f];
        // Disjoint slices per core, handed out in chunks.
        {
            let mut wslices: Vec<&mut [u64]> = Vec::with_capacity(cores.len());
            let mut rest: &mut [u64] = &mut words;
            for i in 0..cores.len() {
                let (a, b) = rest.split_at_mut(word_off[i + 1] - word_off[i]);
                wslices.push(a);
                rest = b;
            }
            let mut fslices: Vec<&mut [u8]> = Vec::with_capacity(cores.len());
            let mut rest: &mut [u8] = &mut exc;
            for i in 0..cores.len() {
                let (a, b) = rest.split_at_mut(flag_off[i + 1] - flag_off[i]);
                fslices.push(a);
                rest = b;
            }
            let jobs: Vec<(usize, &mut [u64], &mut [u8])> = wslices.into_iter().zip(fslices).enumerate().map(|(i, (a, b))| (i, a, b)).collect();
            let chunk = jobs.len().div_ceil(threads.max(1) * 16).max(1);
            let queue = std::sync::Mutex::new(jobs.into_iter());
            std::thread::scope(|sc| {
                for _ in 0..threads.max(1) {
                    sc.spawn(|| loop {
                        let batch: Vec<_> = queue.lock().unwrap().by_ref().take(chunk).collect();
                        if batch.is_empty() {
                            break;
                        }
                        for (i, ws, fs) in batch {
                            fill(i, ws, fs);
                        }
                    });
                }
            });
        }
        let cap = (cores.len() * 2).next_power_of_two().max(2);
        let mut table = vec![[u64::MAX, 0, 0]; cap];
        for (i, &(key, _, _)) in cores.iter().enumerate() {
            let mut s = (key as usize) & (cap - 1);
            while table[s][0] != u64::MAX {
                s = (s + 1) & (cap - 1);
            }
            table[s] = [key, word_off[i] as u64, flag_off[i] as u64];
        }
        FlankExt { brackets, lengths, index_bytes, table, words, exc }
    }

    pub fn bytes(&self) -> usize {
        self.table.len() * 24 + self.words.len() * 8 + self.exc.len()
    }

    pub fn exceptions(&self) -> usize {
        self.exc.iter().filter(|&&x| x != 0).count()
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let tmp = format!("{path}.tmp");
        let mut w = std::io::BufWriter::with_capacity(1 << 24, std::fs::File::create(&tmp)?);
        w.write_all(&MAGIC)?;
        for v in [self.index_bytes, self.brackets.len() as u64, self.lengths.len() as u64, self.table.len() as u64, self.words.len() as u64, self.exc.len() as u64] {
            w.write_all(&v.to_le_bytes())?;
        }
        for &v in self.brackets.iter().chain(self.lengths.iter()) {
            w.write_all(&(v as u64).to_le_bytes())?;
        }
        let tb = unsafe { std::slice::from_raw_parts(self.table.as_ptr() as *const u8, self.table.len() * 24) };
        let wb = unsafe { std::slice::from_raw_parts(self.words.as_ptr() as *const u8, self.words.len() * 8) };
        for part in [tb, wb, &self.exc[..]] {
            for ch in part.chunks(1 << 26) {
                w.write_all(ch)?;
            }
        }
        w.flush()?;
        drop(w);
        std::fs::rename(tmp, path)
    }

    pub fn load(path: &str) -> std::io::Result<FlankExt> {
        let mut f = std::io::BufReader::with_capacity(1 << 24, std::fs::File::open(path)?);
        let mut magic = [0u8; 8];
        f.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{path} is not a flank side store of this version")));
        }
        let mut u = || -> std::io::Result<u64> {
            let mut b = [0u8; 8];
            f.read_exact(&mut b)?;
            Ok(u64::from_le_bytes(b))
        };
        let (index_bytes, nb, nl, nt, nw, ne) = (u()?, u()?, u()?, u()?, u()?, u()?);
        let brackets = (0..nb).map(|_| u().map(|v| v as usize)).collect::<std::io::Result<Vec<_>>>()?;
        let lengths = (0..nl).map(|_| u().map(|v| v as usize)).collect::<std::io::Result<Vec<_>>>()?;
        let mut table = vec![[0u64; 3]; nt as usize];
        f.read_exact(unsafe { std::slice::from_raw_parts_mut(table.as_mut_ptr() as *mut u8, table.len() * 24) })?;
        let mut words = vec![0u64; nw as usize];
        f.read_exact(unsafe { std::slice::from_raw_parts_mut(words.as_mut_ptr() as *mut u8, words.len() * 8) })?;
        let mut exc = vec![0u8; ne as usize];
        f.read_exact(&mut exc)?;
        Ok(FlankExt { brackets, lengths, index_bytes, table, words, exc })
    }
}

/// The loaded store of the run (`--flankext`).
pub static FLANKEXT: std::sync::OnceLock<FlankExt> = std::sync::OnceLock::new();

#[cfg(test)]
mod tests {
    use super::*;

    fn naive(read: &[u8], kpos: i64, e: usize, win: &[u8]) -> usize {
        (0..win.len()).filter(|&t| {
            let rp = kpos - e as i64 + t as i64;
            rp >= 0 && rp < read.len() as i64 && read[rp as usize] != win[t]
        }).count()
    }

    #[test]
    fn packed_window_matches_naive_and_round_trips() {
        let mut x: u64 = 3;
        let mut rnd = |n: usize| -> Vec<u8> { (0..n).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[(x >> 62) as usize] }).collect() };
        // A toy reference: two contigs; cores at a few positions.
        let refs = vec![rnd(500), rnd(300)];
        let e = 8;
        let wl = KMER + 2 * e;
        let copies: Vec<(usize, usize)> = vec![(0, 20), (0, 100), (1, 3), (1, 290), (0, 460)];
        let cores = vec![(0xABCDu64, copies.len(), e), (0x1234u64, 1, e)];
        let fill = |i: usize, ws: &mut [u64], fs: &mut [u8]| {
            let list: Vec<(usize, usize)> = if i == 0 { copies.clone() } else { vec![(1, 50)] };
            let nw = words_per_copy(e);
            for (j, &(c, p)) in list.iter().enumerate() {
                let r = &refs[c];
                if p < e || p + KMER + e > r.len() {
                    fs[j] = 1;
                    continue;
                }
                if !pack_window(&r[p - e..p + KMER + e], &mut ws[j * nw..(j + 1) * nw]) {
                    fs[j] = 1;
                }
            }
        };
        let st = FlankExt::assemble(vec![256], vec![16, 32], 12345, &cores, 3, &fill);
        let path = std::env::temp_dir().join(format!("flankext-test-{}", std::process::id()));
        st.save(path.to_str().unwrap()).unwrap();
        let st = FlankExt::load(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!((st.index_bytes, st.brackets.clone(), st.lengths.clone()), (12345, vec![256], vec![16, 32]));
        let (w0, f0) = st.lookup(0xABCD).unwrap();
        assert!(st.lookup(0x9999).is_none());
        let nw = words_per_copy(e);
        // A read holding copy 1's window with a mismatch and an N, k-mer near the read start.
        for (kpos, read_from) in [(4i64, 96usize), (30, 70)] {
            let mut read = refs[0][read_from..read_from + 120].to_vec();
            read[10] = if read[10] == b'A' { b'C' } else { b'A' };
            read[40] = b'N';
            for (j, &(c, p)) in copies.iter().enumerate() {
                let is_exc = st.is_exception(f0 + j);
                assert_eq!(is_exc, p < e || p + KMER + e > refs[c].len());
                if is_exc {
                    continue;
                }
                let win = &refs[c][p - e..p + KMER + e];
                assert_eq!(win.len(), wl);
                let rw = ReadWindow::new(&read, kpos, e);
                assert_eq!(rw.mismatches(st.cell(w0 + j * nw, nw)), naive(&read, kpos, e, win), "copy {j}, kpos {kpos}");
            }
        }
        assert_eq!(ext_for(&[256], &[16, 32], 256), 0);
        assert_eq!(ext_for(&[256], &[16, 32], 257), 8);
        assert_eq!(ext_for(&[64, 512], &[16, 32, 64], 600), 24);
    }
}
