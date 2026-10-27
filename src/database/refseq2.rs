//! The reference sequences at 2 bits a base: the bundle's `REFSEQ2` section.
//!
//! A marker database is millions of short sequences (protal r226: 14.5 M, 1.1 kb on average, 15.4
//! GB as bytes); at 2 bits a base they take a quarter. The aligner still sees bytes: `get` decodes
//! a whole sequence into a per-thread scratch buffer, about 0.1 us for a marker, and hands out that
//! slice. The scratch is reset by [`reset_scratch`] at the start of every read (pair), so a slice
//! must not be kept past the read it was asked for -- which no caller does: reference slices are
//! locals of a read's processing.
//!
//! Only references whose sequences are all short are packed ([`MAX_PACKED_LEN`]); decoding a whole
//! multi-megabase genome per read would cost more than it saves, so a whole-genome reference keeps
//! its bytes (`REFSEQ`).
//!
//! Layout, little-endian, from the section start (a multiple of 8):
//!
//! ```text
//! MAGIC(8) | n u64 | runs u64 | packed_bytes u64
//! starts   (n + 1) u64   first base of sequence i in the packed bases; a multiple of 4, so every
//!                         sequence starts on a byte; starts[n] is the end
//! lens     n u32, zero-padded to 8 bytes
//! packed   packed_bytes  4 bases a byte, the first in the low bits: A 0, C 1, G 2, T 3
//! runs     runs x (pos u64 | len u32 | byte u8 | 0 0 0)
//!          every run of a byte that is not A, C, G or T (N, IUPAC, lowercase), by position in the
//!          packed bases; it is packed as A and put back when decoded
//! ```

use std::cell::RefCell;
use std::io::Write;

use flexmap::memory::SharedBytes;

pub const MAGIC: [u8; 8] = *b"FREFSQ2\0";
/// Sequences up to this long are packed; a reference with a longer one keeps its bytes.
pub const MAX_PACKED_LEN: usize = 1 << 16;
const RUN_BYTES: usize = 16;

fn code(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Whether the sequences `seq(0..n)` are all short enough to be packed.
pub fn packable<'s>(n: usize, seq: &dyn Fn(usize) -> &'s [u8]) -> bool {
    n > 0 && (0..n).all(|i| seq(i).len() <= MAX_PACKED_LEN)
}

/// Writes the section for sequences `seq(0..n)`; returns its length in bytes.
pub fn write_to<'s, W: Write>(f: &mut W, n: usize, seq: &dyn Fn(usize) -> &'s [u8]) -> std::io::Result<u64> {
    let mut starts = Vec::with_capacity(n + 1);
    let mut at = 0u64;
    for i in 0..n {
        starts.push(at);
        at += (seq(i).len() as u64 + 3) / 4 * 4;
    }
    starts.push(at);
    let packed_bytes = at / 4;

    // The runs of bytes that are not ACGT, by packed position.
    let mut runs: Vec<(u64, u32, u8)> = Vec::new();
    for i in 0..n {
        let s = seq(i);
        let mut j = 0;
        while j < s.len() {
            if code(s[j]).is_some() {
                j += 1;
                continue;
            }
            let b = s[j];
            let from = j;
            while j < s.len() && s[j] == b && (j - from) < u32::MAX as usize {
                j += 1;
            }
            runs.push((starts[i] + from as u64, (j - from) as u32, b));
        }
    }

    f.write_all(&MAGIC)?;
    for x in [n as u64, runs.len() as u64, packed_bytes] {
        f.write_all(&x.to_le_bytes())?;
    }
    for s in &starts {
        f.write_all(&s.to_le_bytes())?;
    }
    for i in 0..n {
        f.write_all(&(seq(i).len() as u32).to_le_bytes())?;
    }
    let lens_pad = (8 - (4 * n) % 8) % 8;
    f.write_all(&vec![0u8; lens_pad])?;
    let mut buf = Vec::with_capacity(1 << 20);
    for i in 0..n {
        for chunk in seq(i).chunks(4) {
            let mut byte = 0u8;
            for (k, &b) in chunk.iter().enumerate() {
                byte |= code(b).unwrap_or(0) << (2 * k);
            }
            buf.push(byte);
        }
        if buf.len() >= 1 << 20 {
            f.write_all(&buf)?;
            buf.clear();
        }
    }
    f.write_all(&buf)?;
    let packed_pad = (8 - packed_bytes % 8) % 8;
    f.write_all(&vec![0u8; packed_pad as usize])?;
    for &(pos, len, byte) in &runs {
        f.write_all(&pos.to_le_bytes())?;
        f.write_all(&len.to_le_bytes())?;
        f.write_all(&[byte, 0, 0, 0])?;
    }
    Ok(32 + 8 * (n as u64 + 1) + 4 * n as u64 + lens_pad as u64 + packed_bytes + packed_pad + (RUN_BYTES * runs.len()) as u64)
}

/// ASCII of the 4 bases of every packed byte, as one little-endian u32.
const LUT: [u32; 256] = {
    let acgt = [b'A', b'C', b'G', b'T'];
    let mut t = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = (acgt[i & 3] as u32)
            | (acgt[(i >> 2) & 3] as u32) << 8
            | (acgt[(i >> 4) & 3] as u32) << 16
            | (acgt[(i >> 6) & 3] as u32) << 24;
        i += 1;
    }
    t
};

/// A REFSEQ2 section, read in: views into bytes it keeps alive.
#[derive(Clone)]
pub struct RefSeq2 {
    _bytes: SharedBytes,
    n: usize,
    starts: *const u64,
    lens: *const u32,
    packed: *const u8,
    runs: *const u8,
    n_runs: usize,
}

// Read-only views into bytes the struct keeps alive.
unsafe impl Send for RefSeq2 {}
unsafe impl Sync for RefSeq2 {}

impl RefSeq2 {
    /// The section at `offset` (a multiple of 8) of `bytes`, `len` bytes long; None if it is not one.
    pub fn from_shared(bytes: &SharedBytes, offset: usize, len: usize) -> Option<Self> {
        let b = bytes.bytes().get(offset..offset + len)?;
        if b.len() < 32 || b[..8] != MAGIC {
            return None;
        }
        let u64at = |at: usize| u64::from_le_bytes(b[at..at + 8].try_into().unwrap());
        let (n, n_runs, packed_bytes) = (u64at(8) as usize, u64at(16) as usize, u64at(24) as usize);
        let starts = 32;
        let lens = starts + 8 * (n + 1);
        let packed = lens + 4 * n + (8 - (4 * n) % 8) % 8;
        let runs = packed + packed_bytes + (8 - packed_bytes % 8) % 8;
        if runs + RUN_BYTES * n_runs > b.len() {
            return None;
        }
        let base = b.as_ptr();
        Some(Self {
            _bytes: bytes.clone(),
            n,
            starts: unsafe { base.add(starts) } as *const u64,
            lens: unsafe { base.add(lens) } as *const u32,
            packed: unsafe { base.add(packed) },
            runs: unsafe { base.add(runs) },
            n_runs,
        })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    /// The length of sequence `id`.
    #[inline]
    pub fn seq_len(&self, id: usize) -> usize {
        unsafe { *self.lens.add(id) as usize }
    }

    #[inline]
    fn run(&self, i: usize) -> (u64, usize, u8) {
        let p = unsafe { self.runs.add(RUN_BYTES * i) };
        let pos = unsafe { (p as *const u64).read_unaligned() };
        let len = unsafe { (p.add(8) as *const u32).read_unaligned() } as usize;
        (pos, len, unsafe { *p.add(12) })
    }

    /// Decodes sequence `id` into `out` (exactly its length).
    pub fn decode(&self, id: usize, out: &mut [u8]) {
        let len = self.seq_len(id);
        debug_assert_eq!(out.len(), len);
        let start = unsafe { *self.starts.add(id) };
        let packed = unsafe { std::slice::from_raw_parts(self.packed.add((start / 4) as usize), (len + 3) / 4) };
        let full = len / 4;
        for (k, &byte) in packed[..full].iter().enumerate() {
            out[4 * k..4 * k + 4].copy_from_slice(&LUT[byte as usize].to_le_bytes());
        }
        if full < packed.len() {
            let tail = LUT[packed[full] as usize].to_le_bytes();
            out[4 * full..].copy_from_slice(&tail[..len - 4 * full]);
        }
        if self.n_runs == 0 {
            return;
        }
        // The runs inside [start, start + len): the first whose end is past start, on.
        let end = start + len as u64;
        let (mut lo, mut hi) = (0usize, self.n_runs);
        while lo < hi {
            let m = lo + (hi - lo) / 2;
            let (pos, l, _) = self.run(m);
            if pos + l as u64 <= start { lo = m + 1 } else { hi = m }
        }
        for i in lo..self.n_runs {
            let (pos, l, byte) = self.run(i);
            if pos >= end {
                break;
            }
            let from = pos.max(start) - start;
            let to = (pos + l as u64).min(end) - start;
            out[from as usize..to as usize].fill(byte);
        }
    }

    /// Sequence `id`, decoded into this thread's scratch (valid until the next [`reset_scratch`]
    /// on this thread). A sequence already decoded since then is handed out again, not decoded.
    #[inline]
    pub fn get(&self, id: usize) -> &[u8] {
        SCRATCH.with(|s| {
            let mut s = s.borrow_mut();
            if let Some(&(_, ptr, len)) = s.seen.iter().rev().find(|e| e.0 == id) {
                return unsafe { std::slice::from_raw_parts(ptr, len) };
            }
            let len = self.seq_len(id);
            let out = s.alloc(len);
            self.decode(id, unsafe { std::slice::from_raw_parts_mut(out, len) });
            s.seen.push((id, out as *const u8, len));
            unsafe { std::slice::from_raw_parts(out as *const u8, len) }
        })
    }
}

/// Decoded sequences of the current read: chunks that are never moved or shrunk while a read is
/// processed, so a slice handed out stays valid until [`reset_scratch`].
struct Scratch {
    chunks: Vec<Box<[u8]>>,
    chunk: usize,
    used: usize,
    seen: Vec<(usize, *const u8, usize)>,
}

const CHUNK: usize = 1 << 20;

impl Scratch {
    fn alloc(&mut self, len: usize) -> *mut u8 {
        loop {
            if self.chunk == self.chunks.len() {
                self.chunks.push(vec![0u8; len.max(CHUNK)].into_boxed_slice());
                self.used = 0;
            }
            let c = &mut self.chunks[self.chunk];
            if self.used + len <= c.len() {
                let p = unsafe { c.as_mut_ptr().add(self.used) };
                self.used += len;
                return p;
            }
            self.chunk += 1;
            self.used = 0;
        }
    }
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = const { RefCell::new(Scratch { chunks: Vec::new(), chunk: 0, used: 0, seen: Vec::new() }) };
}

/// Forgets this thread's decoded sequences: called when a read (pair) starts, after which no slice
/// of the previous read's is used. Keeps the first chunk for the next read; frees the rest.
#[inline]
pub fn reset_scratch() {
    SCRATCH.with(|s| {
        let mut s = s.borrow_mut();
        s.seen.clear();
        s.chunks.truncate(1);
        s.chunk = 0;
        s.used = 0;
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(seqs: &[Vec<u8>]) {
        let mut buf = Vec::new();
        let n = write_to(&mut buf, seqs.len(), &|i| &seqs[i]).unwrap();
        assert_eq!(n as usize, buf.len());
        // At offset 8 of a file, as at a section offset inside a bundle.
        let path = std::env::temp_dir().join(format!("refseq2-{}-{}", std::process::id(), seqs.len()));
        let mut file_bytes = vec![0u8; 8];
        file_bytes.extend_from_slice(&buf);
        std::fs::write(&path, &file_bytes).unwrap();
        let bytes = SharedBytes::read(&std::fs::File::open(&path).unwrap(), 0, file_bytes.len(), 2).unwrap();
        std::fs::remove_file(&path).unwrap();
        let r = RefSeq2::from_shared(&bytes, 8, buf.len()).expect("a REFSEQ2 section");
        assert_eq!(r.len(), seqs.len());
        reset_scratch();
        for (i, s) in seqs.iter().enumerate() {
            assert_eq!(r.get(i), &s[..], "sequence {i}");
        }
        // Asked again within the read: the same bytes, not a second copy.
        let a = r.get(0).as_ptr();
        assert_eq!(r.get(0).as_ptr(), a);
        reset_scratch();
    }

    #[test]
    fn decodes_every_sequence_with_its_runs() {
        let mut seqs = vec![
            b"ACGT".to_vec(),
            b"A".to_vec(),
            b"".to_vec(),
            b"NNNNACGTNNACGTTTGCAN".to_vec(),
            b"GATTACAGATTACAGAT".to_vec(),
            b"acgtNRYacgt".to_vec(),
            vec![b'N'; 9],
        ];
        // Long random ones crossing chunk sizes, with sparse Ns.
        let mut x = 12345u64;
        for len in [1000usize, 4095, 65536] {
            let s: Vec<u8> = (0..len)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    if x % 997 == 0 { b'N' } else { b"ACGT"[(x % 4) as usize] }
                })
                .collect();
            seqs.push(s);
        }
        roundtrip(&seqs);
        roundtrip(&[b"TTTT".to_vec()]);
    }
}
