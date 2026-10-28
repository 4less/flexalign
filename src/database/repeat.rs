//! The repeat index: every position of the cores the main index drops, sorted by flank.
//!
//! The index build sets a core with more than `--max-range-size` positions to zero (flexmap
//! `FMKeys::build`), so at query time it reads like a core absent from the genome. On a human
//! genome that is 82 k cores holding 108 M positions -- young LINE-1, alpha satellite -- and for a
//! read from such a repeat most of its cores return nothing; the rare ones left point to some
//! copies, often not the true one.
//!
//! Here those positions are kept, each with its 16-base flank, sorted by `(key, flank)`. A read's
//! core then looks up the run of positions whose flank equals the read's: the positions where its
//! whole 31-mer occurs. That is what the main index's flank resolution would keep at distance 0,
//! found by a binary search instead of a scan. A run of up to `--repeat-bucket-max` positions is
//! handed to the seed extractor as an ordinary headered range (`AnyRange::Plain`); a longer one, or
//! none, leaves the core unresolved.
//!
//! The set of keys is exactly the complement of the main index's: the same reference pass, the same
//! core selection (canonical core of each K-mer, closed syncmer, `mask_coremer` key) and the same
//! cap. It is kept next to the index as `<index>.repeat`, stamped with the index file's size and
//! modification time and the cap; a mismatch rebuilds it. Deleting it is always safe.

use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::UNIX_EPOCH;

use bioreader::fasta_byte_reader::FastaByteReader;
use bioreader::fasta_reader::FastaReader;
use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::memory::prefetch;
use flexmap::range::AnyRange;
use flexmap::values::{HeaderSeq, VCell, VRange};
use flexmap::VD;
use kmerrs::consecutive::kmer::KmerIter;
use kmerrs::kmer_utils::mask_coremer;
use kmerrs::minimizer::context_free::Minimizer;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

const MAGIC: [u8; 8] = *b"FAREPT02";
/// A sequence is cut into pieces of this many K-mer starts, so one chromosome spreads over threads.
const PIECE: usize = 8 << 20;

/// What the repeat index says about a core the main index does not have.
pub enum RepeatLookup<'a> {
    /// Not a dropped core: absent from the reference (or no repeat index).
    NotRepeat,
    /// The positions whose flank equals the read's, as a headered range.
    Range(AnyRange<'a>),
    /// A dropped core whose exact-flank run is empty or longer than the bucket limit.
    Unresolved,
}

impl RepeatLookup<'_> {
    /// A range of `n` positions, or (for n = 0) not a range.
    pub fn is_range_of(&self, n: usize) -> bool {
        match self {
            RepeatLookup::Range(r) => r.len() == n,
            _ => n == 0,
        }
    }
}

pub struct RepeatIndex {
    cap: u64,
    stamp: (u64, u64),
    /// Sorted dropped keys and, for key `i`, its entries `starts[i]..starts[i + 1]`.
    keys: Vec<u32>,
    starts: Vec<u64>,
    /// Per entry, sorted by flank within a key: the flank cell and the value (`VD::set(ref, pos)`).
    flanks: Vec<u32>,
    cells: Vec<VCell>,
    /// Open addressing over `keys`: `(key + 1) << 32 | index`, 0 empty.
    table: Vec<u64>,
    table_shift: u32,
    /// Per key, `dir_first << 8 | d`: its directory is `dir[dir_first ..= dir_first + 2^d]`, over
    /// its distinct flanks hashed into 2^d slots; slot `t`'s runs are `runs[dir[t] .. dir[t + 1]]`.
    /// The copies of a young repeat mostly share their flanks, so the flank values themselves
    /// cluster; their hash does not, and a slot holds a few distinct flanks however alike the copies
    /// are. A lookup is then the key, its slot and a scan of a few runs: a handful of cache misses,
    /// which `lookup_batch` overlaps across a read's cores.
    meta: Vec<u64>,
    dir: Vec<u32>,
    /// Per distinct (key, flank): the flank, the first of its entries and how many there are.
    runs: Vec<[u32; 3]>,
}

/// Distinct flanks per directory slot a key's directory aims at.
const DIR_SLOT: usize = 4;

fn stamp_of(index_file: &str) -> std::io::Result<(u64, u64)> {
    let meta = std::fs::metadata(index_file)?;
    let modified = meta.modified()?.duration_since(UNIX_EPOCH).map(|d| d.as_nanos() as u64).unwrap_or(0);
    Ok((meta.len(), modified))
}

pub fn path_for(index_file: &str) -> String {
    format!("{index_file}.repeat")
}

impl RepeatIndex {
    /// The dropped core `key`'s positions whose flank is `flank`, as a range of at most `max`.
    #[inline]
    pub fn lookup(&self, key: u64, flank: u32, max: usize) -> RepeatLookup<'_> {
        let Some(i) = self.find(key as u32) else { return RepeatLookup::NotRepeat };
        let (ua, ub) = self.slot_runs(i, flank);
        self.run(ua, ub, flank, max)
    }

    /// `f(j, lookup(queries[j]))` for every query that is a repeat key, in stages over the queries
    /// so that their cache misses (key table, key, directory slot, runs) overlap instead of
    /// following one another: a read's repeat cores are looked up together.
    pub fn lookup_batch<'s>(&'s self, queries: &[(u64, u32)], max: usize, f: &mut dyn FnMut(usize, RepeatLookup<'s>)) {
        const B: usize = 32;
        for (c, chunk) in queries.chunks(B).enumerate() {
            let mut idx = [usize::MAX; B];
            let mut win = [(0usize, 0usize); B];
            for &(k, _) in chunk {
                prefetch(self.table.as_ptr().wrapping_add(Self::home(k as u32, self.table_shift)) as *const u8);
            }
            for (j, &(k, _)) in chunk.iter().enumerate() {
                if let Some(i) = self.find(k as u32) {
                    idx[j] = i;
                    prefetch(self.meta.as_ptr().wrapping_add(i) as *const u8);
                }
            }
            for (j, &(_, flank)) in chunk.iter().enumerate() {
                if idx[j] != usize::MAX {
                    prefetch(self.dir.as_ptr().wrapping_add(self.slot_of(idx[j], flank)) as *const u8);
                }
            }
            for (j, &(_, flank)) in chunk.iter().enumerate() {
                if idx[j] != usize::MAX {
                    win[j] = self.slot_runs(idx[j], flank);
                    prefetch(self.runs.as_ptr().wrapping_add(win[j].0) as *const u8);
                }
            }
            for (j, &(_, flank)) in chunk.iter().enumerate() {
                if idx[j] != usize::MAX {
                    let r = self.run(win[j].0, win[j].1, flank, max);
                    if let RepeatLookup::Range(AnyRange::Plain(v)) = &r {
                        prefetch(v.positions.as_ptr() as *const u8);
                        prefetch(v.header.unwrap_or_default().as_ptr() as *const u8);
                    }
                    f(c * B + j, r);
                }
            }
        }
    }

    /// The index in `dir` of key `i`'s slot for `flank`.
    #[inline(always)]
    fn slot_of(&self, i: usize, flank: u32) -> usize {
        let m = self.meta[i];
        (m >> 8) as usize + Self::slot(flank, (m & 0xFF) as u32)
    }

    #[inline(always)]
    fn slot(flank: u32, d: u32) -> usize {
        if d == 0 { 0 } else { ((flank as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - d)) as usize }
    }

    /// The runs `[ua, ub)` of key `i`'s slot for `flank`.
    #[inline(always)]
    fn slot_runs(&self, i: usize, flank: u32) -> (usize, usize) {
        let t = self.slot_of(i, flank);
        (self.dir[t] as usize, self.dir[t + 1] as usize)
    }

    /// The entries of the run with `flank` among runs `[ua, ub)`, if there are 1 to `max`.
    #[inline(always)]
    fn run(&self, ua: usize, ub: usize, flank: u32, max: usize) -> RepeatLookup<'_> {
        let Some(&[_, start, len]) = self.runs[ua..ub].iter().find(|r| r[0] == flank) else {
            return RepeatLookup::Unresolved;
        };
        let (lo, n) = (start as usize, len as usize);
        if n > max {
            return RepeatLookup::Unresolved;
        }
        // HeaderSeq is a transparent u32.
        let header = unsafe { std::slice::from_raw_parts(self.flanks[lo..lo + n].as_ptr() as *const HeaderSeq, n) };
        RepeatLookup::Range(AnyRange::Plain(VRange::new(Some(header), &self.cells[lo..lo + n])))
    }

    #[inline]
    fn find(&self, key: u32) -> Option<usize> {
        let mask = self.table.len() - 1;
        let mut slot = Self::home(key, self.table_shift);
        loop {
            let e = self.table[slot];
            if e == 0 {
                return None;
            }
            if (e >> 32) as u32 == key.wrapping_add(1) {
                return Some(e as u32 as usize);
            }
            slot = (slot + 1) & mask;
        }
    }

    #[inline(always)]
    fn home(key: u32, shift: u32) -> usize {
        ((key as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> shift) as usize
    }

    /// The per-key directories over the runs of equal flanks of the sorted entries (see `meta`).
    fn make_dirs(starts: &[u64], flanks: &[u32]) -> (Vec<u64>, Vec<u32>, Vec<[u32; 3]>) {
        let mut meta = Vec::with_capacity(starts.len());
        let (mut dir, mut runs) = (Vec::new(), Vec::new());
        let mut key_runs: Vec<(usize, [u32; 3])> = Vec::new();
        for w in starts.windows(2) {
            let (s, e) = (w[0] as usize, w[1] as usize);
            key_runs.clear();
            let mut j = s;
            while j < e {
                let mut k = j + 1;
                while k < e && flanks[k] == flanks[j] {
                    k += 1;
                }
                key_runs.push((0, [flanks[j], j as u32, (k - j) as u32]));
                j = k;
            }
            let u = key_runs.len();
            let d = if u <= 2 * DIR_SLOT { 0 } else { (u / DIR_SLOT).ilog2().min(26) };
            for r in key_runs.iter_mut() {
                r.0 = Self::slot(r.1[0], d);
            }
            key_runs.sort_unstable_by_key(|r| (r.0, r.1[0]));
            meta.push((dir.len() as u64) << 8 | d as u64);
            let base = runs.len();
            let mut q = 0usize;
            for t in 0..=(1usize << d) {
                while q < u && key_runs[q].0 < t {
                    q += 1;
                }
                dir.push((base + q) as u32);
            }
            runs.extend(key_runs.iter().map(|r| r.1));
        }
        (meta, dir, runs)
    }

    fn make_table(keys: &[u32]) -> (Vec<u64>, u32) {
        let size = (keys.len() * 2).next_power_of_two().max(16);
        let shift = 64 - size.trailing_zeros();
        let mut table = vec![0u64; size];
        for (i, &k) in keys.iter().enumerate() {
            let mut slot = Self::home(k, shift);
            while table[slot] != 0 {
                slot = (slot + 1) & (size - 1);
            }
            table[slot] = ((k as u64 + 1) << 32) | i as u64;
        }
        (table, shift)
    }

    pub fn entries(&self) -> usize {
        self.cells.len()
    }

    pub fn key_count(&self) -> usize {
        self.keys.len()
    }

    pub fn memory_bytes(&self) -> usize {
        self.keys.len() * 4 + self.starts.len() * 8 + self.flanks.len() * 4 + self.cells.len() * 8 + self.table.len() * 8
            + self.meta.len() * 8 + self.dir.len() * 4 + self.runs.len() * 12
    }

    /// Counts every core of the reference, then keeps the positions of those over `cap`.
    /// `rid_of` maps a FASTA header (up to the first space) to the index's reference id.
    pub fn build<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
        reference: &Path,
        rid_of: &HashMap<String, usize>,
        cap: u64,
        threads: usize,
    ) -> Self {
        Self::build_in_pieces::<K, C, F, S, L>(reference, rid_of, cap, threads, PIECE)
    }

    fn build_in_pieces<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
        reference: &Path,
        rid_of: &HashMap<String, usize>,
        cap: u64,
        threads: usize,
        piece: usize,
    ) -> Self {
        let threads = threads.max(1);
        let seqs = read_reference(reference, rid_of);
        let pieces: Vec<(usize, usize, usize)> = seqs
            .iter()
            .enumerate()
            .flat_map(|(r, (_, s))| (0..s.len()).step_by(piece).map(move |a| (r, a, min(a + piece, s.len()))))
            .collect();

        // Pass 1: positions per key.
        let size = 1usize << (2 * C);
        let mut counts = vec![0u32; size];
        // AtomicU32 has the layout of u32.
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

        // Pass 2: the dropped keys' entries, partitioned by key range so each part sorts alone.
        let parts = threads * 4;
        let part_of = |key: u32| (key as usize * parts) >> (2 * C);
        let shared: Vec<Mutex<Vec<(u64, u64)>>> = (0..parts).map(|_| Mutex::new(Vec::new())).collect();
        let next = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    let mut cs = ClosedSyncmer::<C, S, L>::new();
                    let mut local: Vec<Vec<(u64, u64)>> = vec![Vec::new(); parts];
                    while let Some(&(r, a, b)) = pieces.get(next.fetch_add(1, Ordering::Relaxed)) {
                        let rid = seqs[r].0 as u64;
                        each_core::<K, C, F, S, L>(&seqs[r].1, a, b, &mut cs, |pos, key, flank| {
                            if counts[key as usize] as u64 > cap {
                                local[part_of(key)].push((((key as u64) << 32) | flank as u64, VD::set(rid, pos as u64)));
                            }
                        });
                    }
                    for (p, v) in local.into_iter().enumerate() {
                        shared[p].lock().unwrap().extend(v);
                    }
                });
            }
        });
        drop(counts);
        drop(seqs);

        let mut parts: Vec<Vec<(u64, u64)>> = shared.into_iter().map(|m| m.into_inner().unwrap()).collect();
        std::thread::scope(|scope| {
            let chunk = parts.len().div_ceil(threads);
            for group in parts.chunks_mut(chunk) {
                scope.spawn(move || group.iter_mut().for_each(|v| v.sort_unstable()));
            }
        });

        let n: usize = parts.iter().map(|v| v.len()).sum();
        let (mut keys, mut starts) = (Vec::new(), Vec::new());
        let (mut flanks, mut cells) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for v in parts {
            for (sort_key, cell) in v {
                let key = (sort_key >> 32) as u32;
                if keys.last() != Some(&key) {
                    keys.push(key);
                    starts.push(flanks.len() as u64);
                }
                flanks.push(sort_key as u32);
                cells.push(VCell(cell));
            }
        }
        starts.push(flanks.len() as u64);
        let (table, table_shift) = Self::make_table(&keys);
        let (meta, dir, runs) = Self::make_dirs(&starts, &flanks);
        Self { cap, stamp: (0, 0), keys, starts, flanks, cells, table, table_shift, meta, dir, runs }
    }

    /// Loads `<index_file>.repeat` if it was built from this index with this cap; otherwise builds it
    /// and writes it (when the directory takes it). Returns the index and what was done.
    pub fn load_or_build<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
        index_file: &str,
        reference: &Path,
        rid_of: impl FnOnce() -> HashMap<String, usize>,
        cap: u64,
        threads: usize,
    ) -> (Arc<Self>, &'static str) {
        let stamp = stamp_of(index_file).unwrap_or((0, 0));
        let path = path_for(index_file);
        if let Ok(r) = Self::load(&path) {
            if r.stamp == stamp && r.cap == cap {
                return (Arc::new(r), "read");
            }
        }
        let mut r = Self::build::<K, C, F, S, L>(reference, &rid_of(), cap, threads);
        r.stamp = stamp;
        if let Err(e) = r.save(&path) {
            log::warn!("could not write the repeat index {path}: {e}");
        }
        (Arc::new(r), "built")
    }

    fn save(&self, path: &str) -> std::io::Result<()> {
        let tmp = format!("{path}.tmp{}", std::process::id());
        {
            let mut w = BufWriter::with_capacity(1 << 24, File::create(&tmp)?);
            w.write_all(&MAGIC)?;
            for v in [self.cap, self.stamp.0, self.stamp.1, self.keys.len() as u64, self.cells.len() as u64, self.dir.len() as u64, self.runs.len() as u64] {
                w.write_all(&v.to_le_bytes())?;
            }
            w.write_all(as_bytes(&self.keys))?;
            w.write_all(as_bytes(&self.starts))?;
            w.write_all(as_bytes(&self.flanks))?;
            w.write_all(as_bytes(&self.cells))?;
            w.write_all(as_bytes(&self.meta))?;
            w.write_all(as_bytes(&self.dir))?;
            w.write_all(as_bytes(&self.runs))?;
            w.flush()?;
        }
        std::fs::rename(&tmp, path)
    }

    fn load(path: &str) -> std::io::Result<Self> {
        let mut r = BufReader::with_capacity(1 << 24, File::open(path)?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "not a repeat index"));
        }
        let mut head = [0u64; 7];
        for v in &mut head {
            let mut b = [0u8; 8];
            r.read_exact(&mut b)?;
            *v = u64::from_le_bytes(b);
        }
        let [cap, len, modified, n_keys, n, n_dir, n_runs] = head;
        let mut keys = vec![0u32; n_keys as usize];
        let mut starts = vec![0u64; n_keys as usize + 1];
        let mut flanks = vec![0u32; n as usize];
        let mut cells = vec![VCell(0); n as usize];
        r.read_exact(as_bytes_mut(&mut keys))?;
        r.read_exact(as_bytes_mut(&mut starts))?;
        r.read_exact(as_bytes_mut(&mut flanks))?;
        r.read_exact(as_bytes_mut(&mut cells))?;
        let mut meta = vec![0u64; n_keys as usize];
        let mut dir = vec![0u32; n_dir as usize];
        let mut runs = vec![[0u32; 3]; n_runs as usize];
        r.read_exact(as_bytes_mut(&mut meta))?;
        r.read_exact(as_bytes_mut(&mut dir))?;
        r.read_exact(as_bytes_mut(&mut runs))?;
        let (table, table_shift) = Self::make_table(&keys);
        Ok(Self { cap, stamp: (len, modified), keys, starts, flanks, cells, table, table_shift, meta, dir, runs })
    }
}

/// The plain-old-data arrays above as bytes (u32, u64 and VCell, which is a repr(C) u64).
fn as_bytes<T: Copy>(v: &[T]) -> &[u8] {
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

fn as_bytes_mut<T: Copy>(v: &mut [T]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, std::mem::size_of_val(v)) }
}

/// `f(position, key, flank)` for every selected core whose K-mer starts in `[a, b)` of `seq`, as
/// flexmap's build selects and keys them.
fn each_core<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
    seq: &[u8],
    a: usize,
    b: usize,
    cs: &mut ClosedSyncmer<C, S, L>,
    mut f: impl FnMut(usize, u32, u32),
) {
    let piece = &seq[a..min(seq.len(), b + K - 1)];
    for (pos, kmer_fwd, kmer_rev) in KmerIter::<K, true>::new(piece) {
        let cmer_fwd = kmer_fwd.middle::<C>();
        let cmer_rev = kmer_rev.middle::<C>();
        let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
        let cmer = min(cmer_fwd, cmer_rev);
        if !cs.is_minimizer(cmer.0) {
            continue;
        }
        let key = mask_coremer::<C>(cmer.0) as u32;
        f(a + pos, key, kmer.flanks::<F>().0 as u32);
    }
}

/// The reference's sequences with their index ids, read as the index build reads them.
fn read_reference(path: &Path, rid_of: &HashMap<String, usize>) -> Vec<(usize, Vec<u8>)> {
    let file = File::open(path).unwrap_or_else(|e| panic!("reference {}: {e}", path.display()));
    let buffer_size = 1 << 24;
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size).expect("reference reader")));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);
    let mut record = OwnedFastaRecord::new();
    let mut out = Vec::new();
    while let Some(()) = fasta_reader.load_batch_par(&mut byte_reader).expect("reference batch") {
        while fasta_reader.next(&mut record).is_some() {
            let header = String::from_utf8_lossy(&record.head()[1..]).split(' ').next().unwrap().to_string();
            let rid = *rid_of.get(&header).unwrap_or_else(|| panic!("reference {header} is not in the index"));
            out.push((rid, record.seq().to_vec()));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use flexmap::flexmap::VRangeGetter;

    /// Keys with many distinct flanks (multi-slot directories), some flanks repeated many times:
    /// every flank finds exactly its run, absent flanks find nothing, runs over `max` are unresolved.
    #[test]
    fn directories_find_every_run() {
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let mut entries: Vec<(u64, u64)> = Vec::new();
        for key in [3u32, 77, 1 << 20] {
            for e in 0..2000u64 {
                // Few top-bit patterns and a skewed spread of repeats, as in a young repeat.
                let flank = ((next() % 40) as u32) | (((e % 3) as u32) << 30);
                let copies = if e % 50 == 0 { 30 } else { 1 };
                for c in 0..copies {
                    entries.push((((key as u64) << 32) | flank as u64, e * 100 + c));
                }
            }
        }
        entries.sort_unstable();
        let (mut keys, mut starts, mut flanks, mut cells) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for &(sk, cell) in &entries {
            let key = (sk >> 32) as u32;
            if keys.last() != Some(&key) {
                keys.push(key);
                starts.push(flanks.len() as u64);
            }
            flanks.push(sk as u32);
            cells.push(VCell(cell));
        }
        starts.push(flanks.len() as u64);
        let (table, table_shift) = RepeatIndex::make_table(&keys);
        let (meta, dir, runs) = RepeatIndex::make_dirs(&starts, &flanks);
        assert!(meta.iter().any(|m| m & 0xFF > 0), "a directory with slots");
        let rep = RepeatIndex { cap: 0, stamp: (0, 0), keys, starts, flanks, cells, table, table_shift, meta, dir, runs };
        for (i, &key) in rep.keys.iter().enumerate() {
            let (s, e) = (rep.starts[i] as usize, rep.starts[i + 1] as usize);
            for q in (0..64u32).flat_map(|f| [f, f | 1 << 30, f | 2 << 30, f | 3 << 30]) {
                let n = (s..e).filter(|&j| rep.flanks[j] == q).count();
                for max in [1usize, 29, 30, 1000] {
                    let r = rep.lookup(key as u64, q, max);
                    assert!(r.is_range_of(if n <= max { n } else { 0 }), "key {key} flank {q:x} max {max}: {n} entries");
                }
            }
        }
        assert!(matches!(rep.lookup(5, 0, 10), RepeatLookup::NotRepeat));
    }

    /// Every key's entries in the repeat index (built with cap 0, so it keeps every key) equal the
    /// main index's values and flank cells for that key: same positions, ids and flanks.
    #[test]
    fn repeat_entries_match_the_main_index() {
        const K: usize = 31;
        const C: usize = 15;
        const F: usize = 16;
        const S: usize = 4;
        const L: usize = C - S + 1;
        let dir = std::env::temp_dir().join(format!("flexalign_repeat_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("ref.fna");
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut base = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x & 3) as usize]
        };
        let unit: Vec<u8> = (0..300).map(|_| base()).collect();
        let mut text = String::new();
        for r in 0..3 {
            // Random sequence with copies of one unit (some with a change) and a run of Ns.
            let mut s: Vec<u8> = (0..2000).map(|_| base()).collect();
            for c in 0..4 {
                let mut u = unit.clone();
                if c == 1 {
                    u[150] = if u[150] == b'A' { b'C' } else { b'A' };
                }
                s.extend_from_slice(&u);
            }
            s.extend_from_slice(b"NNNNN");
            s.extend((0..500).map(|_| base()));
            text.push_str(&format!(">seq{r} test\n{}\n", String::from_utf8(s).unwrap()));
        }
        std::fs::write(&fasta, text).unwrap();

        let (main, reference2id, _) =
            flexmap::build::default_build::<K, C, F, S, L, 16, 2>(&fasta, usize::MAX).expect("build");
        let rep = RepeatIndex::build::<K, C, F, S, L>(&fasta, &reference2id, 0, 3);
        assert!(rep.key_count() > 0);
        // Sequences cut into pieces of 97 K-mer starts give the same index.
        let cut = RepeatIndex::build_in_pieces::<K, C, F, S, L>(&fasta, &reference2id, 0, 3, 97);
        assert_eq!((&cut.keys, &cut.starts, &cut.flanks), (&rep.keys, &rep.starts, &rep.flanks));
        assert!(cut.cells.iter().zip(&rep.cells).all(|(a, b)| a.0 == b.0));
        let mut checked = 0;
        for (i, &key) in rep.keys.iter().enumerate() {
            let (s, e) = (rep.starts[i] as usize, rep.starts[i + 1] as usize);
            let mut want: Vec<(u32, u64)> = Vec::new();
            let range = main.get_vrange(key as u64).expect("a key of the reference is in the main index");
            for j in 0..range.positions.len() {
                let flank = range.header.map(|h| h[j].get()).unwrap_or(0);
                want.push((flank, range.positions[j].0));
            }
            let mut got: Vec<(u32, u64)> = (s..e).map(|j| (rep.flanks[j], rep.cells[j].0)).collect();
            if range.header.is_none() {
                // Headerless keys store no flank; compare the positions only.
                got.iter_mut().for_each(|g| g.0 = 0);
            }
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(got, want, "key {key}");
            checked += 1;
        }
        assert!(checked > 100, "{checked} keys checked");

        // A lookup by the flank of an entry returns exactly the entries with that flank.
        let i = rep.keys.len() / 2;
        let (s, e) = (rep.starts[i] as usize, rep.starts[i + 1] as usize);
        let flank = rep.flanks[s];
        match rep.lookup(rep.keys[i] as u64, flank, usize::MAX) {
            RepeatLookup::Range(r) => assert_eq!(r.len(), (s..e).filter(|&j| rep.flanks[j] == flank).count()),
            _ => panic!("expected a range"),
        }
        assert!(matches!(rep.lookup(rep.keys[i] as u64, flank, 0), RepeatLookup::Unresolved));
        // Every key, with each of its flanks and a few that are not there: the run the directory
        // finds is the run a scan of all its entries finds.
        for (i, &key) in rep.keys.iter().enumerate() {
            let (s, e) = (rep.starts[i] as usize, rep.starts[i + 1] as usize);
            let probes: Vec<u32> = rep.flanks[s..e].iter().flat_map(|&f| [f, f ^ 1, f.wrapping_add(1 << 30)]).collect();
            for q in probes {
                let n = (s..e).filter(|&j| rep.flanks[j] == q).count();
                let mut batch = Vec::new();
                rep.lookup_batch(&[(key as u64, q), (u64::MAX >> 34, q)], 1000, &mut |j, r| batch.push((j, r.is_range_of(n))));
                assert_eq!(batch, vec![(0, true)], "key {key} flank {q:x}: batch");
                for max in [1usize, 2, 1000] {
                    match rep.lookup(key as u64, q, max) {
                        RepeatLookup::Range(r) => assert!(n > 0 && n <= max && r.len() == n, "key {key} flank {q:x}"),
                        RepeatLookup::Unresolved => assert!(n == 0 || n > max, "key {key} flank {q:x}"),
                        RepeatLookup::NotRepeat => panic!("key {key} is a repeat key"),
                    }
                }
            }
        }

        let path = dir.join("idx.repeat");
        rep.save(path.to_str().unwrap()).unwrap();
        let back = RepeatIndex::load(path.to_str().unwrap()).unwrap();
        assert_eq!((back.keys.clone(), back.starts.clone(), back.flanks.clone()), (rep.keys.clone(), rep.starts.clone(), rep.flanks.clone()));
        assert!(back.cells.iter().zip(&rep.cells).all(|(a, b)| a.0 == b.0));
        std::fs::remove_dir_all(&dir).ok();
    }
}
