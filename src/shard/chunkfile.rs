//! Chunked, self-indexed blob files -- the container behind `reads.bin` and `evidence.s{N}.bin`.
//!
//! A chunk is one unit of work (one bioreader batch, ~48k pairs). Each file reserves an index at
//! its head, sized `n_chunks * 12 B`, and fills it in at pass end. Writers reserve space with an
//! atomic bump and `write_at`, so blobs land in **arbitrary order** and threads never hold a lock
//! during IO -- the index names each blob explicitly, so their file order is irrelevant.
//!
//! ```text
//! +----------------+-----------+-----------+-----+-------------+
//! | header + index | chunk k   | chunk 0   | ... | chunk 12    |   (any order)
//! | 32 + n*12 B    |           |           |     |             |
//! +----------------+-----------+-----------+-----+-------------+
//!      index[c] = (offset u64, len u32)
//! ```
//!
//! This is what makes the rejoin a set of independent positioned reads rather than a k-way merge
//! over synchronised streams (ProjectShard.md §9). Readers use `read_at` (pread): workers read
//! different offsets of the same fd concurrently, and `seek`+`read` would race on the shared file
//! cursor.

use std::fs::{File, OpenOptions};
use std::io;
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const MAGIC: u64 = 0x464C58_5348_4400; // "FLXSHD"
const VERSION: u32 = 1;
const HEADER_BYTES: u64 = 32;
const INDEX_ENTRY_BYTES: u64 = 12;

/// Where one chunk's blob lives. `len == 0` means the chunk produced no evidence (every read in
/// it missed this shard), which is normal, not an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ChunkExtent {
    pub offset: u64,
    pub len: u32,
}

impl ChunkExtent {
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

fn index_bytes(n_chunks: u64) -> u64 {
    n_chunks * INDEX_ENTRY_BYTES
}

fn data_start(n_chunks: u64) -> u64 {
    HEADER_BYTES + index_bytes(n_chunks)
}

/// Writes a chunk file. Cloneable handles share one atomic tail, so every thread in a pass can
/// write concurrently without coordinating.
pub struct ChunkFileWriter {
    file: File,
    n_chunks: u64,
    tail: AtomicU64,
    index: Vec<AtomicU64>,
    lens: Vec<AtomicU64>,
}

impl ChunkFileWriter {
    /// `n_chunks` must be known up front -- it sizes the reserved index. For evidence files it
    /// comes from `reads.bin`, which phase 0 wrote first.
    pub fn create(path: impl AsRef<Path>, n_chunks: u64) -> io::Result<Self> {
        let file = OpenOptions::new().create(true).write(true).truncate(true).open(path)?;
        let start = data_start(n_chunks);
        file.set_len(start)?;
        Ok(Self {
            file,
            n_chunks,
            tail: AtomicU64::new(start),
            index: (0..n_chunks).map(|_| AtomicU64::new(0)).collect(),
            lens: (0..n_chunks).map(|_| AtomicU64::new(0)).collect(),
        })
    }

    pub fn n_chunks(&self) -> u64 {
        self.n_chunks
    }

    /// Reserves space, writes the blob, and records the extent. Safe to call from any thread; no
    /// lock is held during the write.
    pub fn write_chunk(&self, chunk_id: u64, buf: &[u8]) -> io::Result<ChunkExtent> {
        if chunk_id >= self.n_chunks {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("chunk_id {} out of range (n_chunks {})", chunk_id, self.n_chunks),
            ));
        }
        if buf.is_empty() {
            self.lens[chunk_id as usize].store(1, Ordering::Release); // mark written, len 0
            return Ok(ChunkExtent::default());
        }

        let offset = self.tail.fetch_add(buf.len() as u64, Ordering::AcqRel);
        self.file.write_at(buf, offset)?;

        self.index[chunk_id as usize].store(offset, Ordering::Release);
        self.lens[chunk_id as usize].store(buf.len() as u64 + 1, Ordering::Release);
        Ok(ChunkExtent { offset, len: buf.len() as u32 })
    }

    /// Writes the header and index, then fsyncs. Fails if any chunk was never written -- a pass
    /// is complete only when its index covers every chunk (ProjectShard.md §8).
    pub fn finish(self) -> io::Result<()> {
        let mut missing = Vec::new();
        for (i, l) in self.lens.iter().enumerate() {
            if l.load(Ordering::Acquire) == 0 {
                missing.push(i);
            }
        }
        if !missing.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{} chunk(s) never written, first {:?} -- pass is incomplete",
                    missing.len(),
                    &missing[..missing.len().min(8)]
                ),
            ));
        }

        // Data before index, so a crash leaves a short-but-valid file rather than an index
        // pointing at bytes that never landed.
        self.file.sync_data()?;

        let mut head = Vec::with_capacity(data_start(self.n_chunks) as usize);
        head.extend_from_slice(&MAGIC.to_le_bytes());
        head.extend_from_slice(&VERSION.to_le_bytes());
        head.extend_from_slice(&0u32.to_le_bytes()); // reserved
        head.extend_from_slice(&self.n_chunks.to_le_bytes());
        head.extend_from_slice(&0u64.to_le_bytes()); // reserved
        debug_assert_eq!(head.len() as u64, HEADER_BYTES);

        for i in 0..self.n_chunks as usize {
            let len = self.lens[i].load(Ordering::Acquire) - 1;
            head.extend_from_slice(&self.index[i].load(Ordering::Acquire).to_le_bytes());
            head.extend_from_slice(&(len as u32).to_le_bytes());
        }

        self.file.write_at(&head, 0)?;
        self.file.sync_all()
    }
}

/// Reads a chunk file. Shared across workers; every read is positioned.
pub struct ChunkFileReader {
    file: File,
    index: Vec<ChunkExtent>,
}

impl ChunkFileReader {
    pub fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let file = File::open(path)?;

        let mut head = [0u8; HEADER_BYTES as usize];
        file.read_exact_at(&mut head, 0)?;

        let magic = u64::from_le_bytes(head[0..8].try_into().unwrap());
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a flexalign chunk file"));
        }
        let version = u32::from_le_bytes(head[8..12].try_into().unwrap());
        if version != VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("chunk file version {}, expected {}", version, VERSION),
            ));
        }
        let n_chunks = u64::from_le_bytes(head[16..24].try_into().unwrap());

        let mut raw = vec![0u8; index_bytes(n_chunks) as usize];
        file.read_exact_at(&mut raw, HEADER_BYTES)?;

        let index = raw
            .chunks_exact(INDEX_ENTRY_BYTES as usize)
            .map(|e| ChunkExtent {
                offset: u64::from_le_bytes(e[0..8].try_into().unwrap()),
                len: u32::from_le_bytes(e[8..12].try_into().unwrap()),
            })
            .collect();

        Ok(Self { file, index })
    }

    pub fn n_chunks(&self) -> usize {
        self.index.len()
    }

    pub fn extent(&self, chunk_id: usize) -> Option<ChunkExtent> {
        self.index.get(chunk_id).copied()
    }

    /// Reads one chunk into `out`. This is the single IO the rejoin does per (chunk, shard).
    pub fn read_chunk(&self, chunk_id: usize, out: &mut Vec<u8>) -> io::Result<()> {
        let e = self.extent(chunk_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("chunk {} out of range ({} chunks)", chunk_id, self.index.len()),
            )
        })?;
        out.clear();
        if e.is_empty() {
            return Ok(());
        }
        out.resize(e.len as usize, 0);
        self.file.read_exact_at(out, e.offset)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn tmp(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("flexalign-chunkfile-{}-{}", std::process::id(), name));
        p
    }

    #[test]
    fn roundtrips_chunks_written_out_of_order() {
        let path = tmp("out-of-order");
        let w = ChunkFileWriter::create(&path, 4).unwrap();
        // Deliberately not in chunk order -- this is what the atomic bump produces.
        w.write_chunk(2, b"chunk-two").unwrap();
        w.write_chunk(0, b"zero").unwrap();
        w.write_chunk(3, b"three!!").unwrap();
        w.write_chunk(1, b"one").unwrap();
        w.finish().unwrap();

        let r = ChunkFileReader::open(&path).unwrap();
        assert_eq!(r.n_chunks(), 4);
        let mut buf = Vec::new();
        for (id, want) in [(0, &b"zero"[..]), (1, b"one"), (2, b"chunk-two"), (3, b"three!!")] {
            r.read_chunk(id, &mut buf).unwrap();
            assert_eq!(buf, want, "chunk {}", id);
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn empty_chunks_are_valid() {
        let path = tmp("empty");
        let w = ChunkFileWriter::create(&path, 3).unwrap();
        w.write_chunk(0, b"data").unwrap();
        w.write_chunk(1, b"").unwrap(); // no read in this chunk hit this shard
        w.write_chunk(2, b"more").unwrap();
        w.finish().unwrap();

        let r = ChunkFileReader::open(&path).unwrap();
        let mut buf = Vec::new();
        r.read_chunk(1, &mut buf).unwrap();
        assert!(buf.is_empty());
        assert!(r.extent(1).unwrap().is_empty());
        r.read_chunk(2, &mut buf).unwrap();
        assert_eq!(buf, b"more");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn incomplete_pass_is_rejected() {
        let path = tmp("incomplete");
        let w = ChunkFileWriter::create(&path, 3).unwrap();
        w.write_chunk(0, b"a").unwrap();
        w.write_chunk(2, b"c").unwrap();
        let err = w.finish().unwrap_err();
        assert!(err.to_string().contains("never written"), "{}", err);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn concurrent_writers_do_not_overlap() {
        let path = tmp("concurrent");
        let n = 64u64;
        let w = Arc::new(ChunkFileWriter::create(&path, n).unwrap());

        std::thread::scope(|s| {
            for t in 0..8u64 {
                let w = Arc::clone(&w);
                s.spawn(move || {
                    for c in (t..n).step_by(8) {
                        // Varying sizes so a bad bump would misalign rather than coincide.
                        let blob = vec![c as u8; 1 + (c as usize % 97)];
                        w.write_chunk(c, &blob).unwrap();
                    }
                });
            }
        });

        Arc::try_unwrap(w).ok().unwrap().finish().unwrap();

        let r = ChunkFileReader::open(&path).unwrap();
        let mut buf = Vec::new();
        for c in 0..n {
            r.read_chunk(c as usize, &mut buf).unwrap();
            let want = vec![c as u8; 1 + (c as usize % 97)];
            assert_eq!(buf, want, "chunk {} corrupted by a concurrent writer", c);
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rejects_foreign_and_out_of_range() {
        let path = tmp("foreign");
        std::fs::write(&path, vec![0u8; 64]).unwrap();
        assert!(ChunkFileReader::open(&path).is_err());
        std::fs::remove_file(&path).ok();

        let path = tmp("range");
        let w = ChunkFileWriter::create(&path, 1).unwrap();
        assert!(w.write_chunk(5, b"x").is_err());
        w.write_chunk(0, b"x").unwrap();
        w.finish().unwrap();
        let r = ChunkFileReader::open(&path).unwrap();
        let mut buf = Vec::new();
        assert!(r.read_chunk(9, &mut buf).is_err());
        std::fs::remove_file(&path).ok();
    }
}
