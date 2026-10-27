//! One file for everything a packed-index run loads: the k-mer index, the reference sequences and
//! the reference names. It is read in one parallel pass into huge pages (or mapped once with
//! `--index-mmap`), and each part is a view into those bytes -- nothing is copied or deserialised.
//!
//! Layout, little-endian:
//!
//! ```text
//! header (4096 bytes)   MAGIC(8) | VERSION u32 | sections u32 | per section: kind u32, 0 u32, offset u64, len u64
//! sections              each at a multiple of ALIGN (2 MiB), zero-padded between
//!   INDEX               flexmap's packed index, as PackedFlexmap::write_to writes it
//!   REFSEQ              REFSEQ_MAGIC(8) | n u64 | offsets (n + 1) u64 | sequences
//!   NAMES               n u64 | offsets (n + 1) u64 | names (UTF-8, no separators)
//! ```
//!
//! ALIGN is the huge-page size, so in a read-in buffer (huge-page aligned) every section starts on
//! a huge page, and in a mapping on a page; the index needs 4096 and the tables 8.

use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};

use flexmap::memory::SharedBytes;

const MAGIC: [u8; 8] = *b"FLXBNDL1";
const VERSION: u32 = 1;
const HEADER: usize = 4096;
const ALIGN: u64 = 2 << 20;

pub const INDEX: u32 = 1;
pub const REFSEQ: u32 = 2;
pub const NAMES: u32 = 3;

/// Writes a bundle to `path` (through `path.tmp`, renamed when complete, so a reader never sees a
/// half-written bundle). `index` writes the packed index and returns its length; `refseq` writes
/// the reference section.
pub fn write<'s>(
    path: &str,
    index: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<u64>,
    refseq: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<u64>,
    names: &dyn Fn(usize) -> &'s str,
    n_names: usize,
) -> std::io::Result<()> {
    let tmp = format!("{path}.tmp");
    let mut out = BufWriter::with_capacity(8 << 20, File::create(&tmp)?);
    out.write_all(&vec![0u8; HEADER])?;
    let mut table = Vec::new();

    pad_to(&mut out, ALIGN)?;
    let at = out.stream_position()?;
    table.push((INDEX, at, index(&mut out)?));

    pad_to(&mut out, ALIGN)?;
    let at = out.stream_position()?;
    table.push((REFSEQ, at, refseq(&mut out)?));

    pad_to(&mut out, ALIGN)?;
    let at = out.stream_position()?;
    table.push((NAMES, at, write_names(&mut out, names, n_names)?));

    let mut header = Vec::with_capacity(HEADER);
    header.extend_from_slice(&MAGIC);
    header.extend_from_slice(&VERSION.to_le_bytes());
    header.extend_from_slice(&(table.len() as u32).to_le_bytes());
    for (kind, offset, len) in &table {
        header.extend_from_slice(&kind.to_le_bytes());
        header.extend_from_slice(&0u32.to_le_bytes());
        header.extend_from_slice(&offset.to_le_bytes());
        header.extend_from_slice(&len.to_le_bytes());
    }
    out.seek(std::io::SeekFrom::Start(0))?;
    out.write_all(&header)?;
    out.flush()?;
    drop(out);
    std::fs::rename(&tmp, path)
}

fn pad_to(out: &mut BufWriter<File>, align: u64) -> std::io::Result<()> {
    let at = out.stream_position()?;
    let end = (at + align - 1) / align * align;
    std::io::copy(&mut std::io::repeat(0).take(end - at), out)?;
    Ok(())
}


fn write_names<'s>(out: &mut BufWriter<File>, names: &dyn Fn(usize) -> &'s str, n: usize) -> std::io::Result<u64> {
    out.write_all(&(n as u64).to_le_bytes())?;
    let mut offset = 0u64;
    out.write_all(&offset.to_le_bytes())?;
    for i in 0..n {
        offset += names(i).len() as u64;
        out.write_all(&offset.to_le_bytes())?;
    }
    for i in 0..n {
        out.write_all(names(i).as_bytes())?;
    }
    Ok(8 * (n as u64 + 2) + offset)
}

/// A bundle file, read in or mapped.
pub struct Bundle {
    pub bytes: SharedBytes,
    /// Whether the bytes are a mapping (`--index-mmap` on a local filesystem) or were read in.
    pub mapped: bool,
    sections: Vec<(u32, usize, usize)>,
}

impl Bundle {
    /// `mmap`: map the file and populate it on `threads` threads (the pages are the page cache's,
    /// shared by every run on the machine). Otherwise -- and always on a network filesystem, see
    /// [`on_network_fs`] -- read it all on `threads` threads into huge pages of this process.
    pub fn open(path: &str, mmap: bool, threads: usize) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len() as usize;
        // Mapped, the pages are populated at once on all threads: cached, that only fills page
        // tables (~25 ms for 1.4 GB); uncached it reads the file as fast as the parallel read does.
        // Left to fault in lazily, k-mer lookups would read it page by page in random order (3-4x
        // slower runs from a cold cache).
        let mmap = mmap && !on_network_fs(&file);
        let bytes = if mmap { SharedBytes::mmap(&file)? } else { SharedBytes::read(&file, 0, len, threads)? };
        if mmap {
            populate(&bytes, threads);
        }
        let sections = Self::parse(&bytes.bytes()[..HEADER.min(bytes.len())])
            .ok_or_else(|| std::io::Error::other(format!("{path}: not a flexalign bundle")))?;
        if sections.iter().any(|&(_, offset, len)| offset + len > bytes.len()) {
            return Err(std::io::Error::other(format!("{path}: truncated bundle")));
        }
        Ok(Self { bytes, mapped: mmap, sections })
    }

    fn parse(header: &[u8]) -> Option<Vec<(u32, usize, usize)>> {
        if header.len() < 16 || header[..8] != MAGIC || u32le(header, 8) != VERSION {
            return None;
        }
        let count = u32le(header, 12) as usize;
        if 16 + count * 24 > header.len() {
            return None;
        }
        Some(
            (0..count)
                .map(|i| {
                    let at = 16 + i * 24;
                    (u32le(header, at), u64le(header, at + 8) as usize, u64le(header, at + 16) as usize)
                })
                .collect(),
        )
    }

    /// The header and the first bytes of the index of the bundle at `path`, read without loading
    /// it: for a compatibility check.
    pub fn peek_index(path: &str) -> Option<Vec<u8>> {
        let mut file = File::open(path).ok()?;
        let mut header = vec![0u8; HEADER];
        file.read_exact(&mut header).ok()?;
        let (_, offset, _) = *Self::parse(&header)?.iter().find(|s| s.0 == INDEX)?;
        let mut head = vec![0u8; 64];
        file.seek(std::io::SeekFrom::Start(offset as u64)).ok()?;
        file.read_exact(&mut head).ok()?;
        Some(head)
    }

    /// The sections from the one of `from` to the end of the file, read alone on `threads` threads
    /// (a sharded run's join needs the reference and the names, not the index before them).
    /// Section offsets are then relative to what was read.
    pub fn open_from(path: &str, from: u32, threads: usize) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len() as usize;
        let mut header = vec![0u8; HEADER];
        std::os::unix::fs::FileExt::read_exact_at(&file, &mut header, 0)?;
        let sections = Self::parse(&header).ok_or_else(|| std::io::Error::other(format!("{path}: not a flexalign bundle")))?;
        let (_, base, _) = *sections
            .iter()
            .find(|s| s.0 == from)
            .ok_or_else(|| std::io::Error::other(format!("{path}: no section {from}")))?;
        if sections.iter().any(|&(_, offset, n)| offset + n > len) {
            return Err(std::io::Error::other(format!("{path}: truncated bundle")));
        }
        let bytes = SharedBytes::read(&file, base as u64, len - base, threads)?;
        let sections = sections.into_iter().filter(|s| s.1 >= base).map(|(k, offset, n)| (k, offset - base, n)).collect();
        Ok(Self { bytes, mapped: false, sections })
    }

    /// (offset, length) in the file of the section of `kind`, from the header alone.
    pub fn section_in_file(path: &str, kind: u32) -> std::io::Result<(u64, u64)> {
        let file = File::open(path)?;
        let mut header = vec![0u8; HEADER];
        std::os::unix::fs::FileExt::read_exact_at(&file, &mut header, 0)?;
        let sections = Self::parse(&header).ok_or_else(|| std::io::Error::other(format!("{path}: not a flexalign bundle")))?;
        sections
            .iter()
            .find(|s| s.0 == kind)
            .map(|&(_, offset, n)| (offset as u64, n as u64))
            .ok_or_else(|| std::io::Error::other(format!("{path}: no section {kind}")))
    }

    /// (offset, length) of the section of `kind`.
    pub fn section(&self, kind: u32) -> Option<(usize, usize)> {
        self.sections.iter().find(|s| s.0 == kind).map(|&(_, offset, len)| (offset, len))
    }
}

/// Whether `file` lives on a network filesystem, where page faults of a mapping go through the
/// filesystem client one page at a time: on Lustre a cached 6.56 GB bundle took 14.97 s mapped and
/// populated on 32 threads, 0.24 s read in parallel. There the bundle is always read.
pub fn on_network_fs(file: &File) -> bool {
    use std::os::fd::AsRawFd;
    const LUSTRE: i64 = 0x0BD0_0BD0;
    const NFS: i64 = 0x6969;
    const GPFS: i64 = 0x4750_4653;
    const CIFS: i64 = 0xFF53_4D42;
    const SMB2: i64 = 0xFE53_4D42;
    const BEEGFS: i64 = 0x1983_0326;
    let mut st: libc::statfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::fstatfs(file.as_raw_fd(), &mut st) } != 0 {
        return false;
    }
    matches!(st.f_type as i64, LUSTRE | NFS | GPFS | CIFS | SMB2 | BEEGFS)
}

/// Faults every page of a mapped bundle in now, on `threads` threads (`MADV_POPULATE_READ`, one
/// range each): from the page cache that only fills page tables; from disk it is one parallel read
/// instead of page-by-page faults in lookup order later.
pub fn populate(bytes: &SharedBytes, threads: usize) {
    let Some(map) = bytes.mapping() else { return };
    const PIECE: usize = 64 << 20;
    let pieces = (map.len() + PIECE - 1) / PIECE;
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads.min(pieces).max(1) {
            scope.spawn(|| loop {
                let piece = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if piece >= pieces {
                    return;
                }
                let from = piece * PIECE;
                let _ = map.advise_range(memmap2::Advice::PopulateRead, from, PIECE.min(map.len() - from));
            });
        }
    });
}

fn u32le(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

fn u64le(b: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(b[at..at + 8].try_into().unwrap())
}

/// The reference names in a bundle: an offset table and the names back to back, used where they
/// lie. Nothing is deserialised, so a load costs nothing however many names there are.
#[derive(Clone)]
pub struct NameTable {
    _bytes: SharedBytes,
    offsets: *const u64,
    data: *const u8,
    n: usize,
}

// Read-only views into bytes the table keeps alive.
unsafe impl Send for NameTable {}
unsafe impl Sync for NameTable {}

impl NameTable {
    pub fn from_bundle(bundle: &Bundle) -> Option<Self> {
        let (offset, len) = bundle.section(NAMES)?;
        let b = &bundle.bytes.bytes()[offset..offset + len];
        let n = u64le(b, 0) as usize;
        if 8 * (n + 2) > len || u64le(b, 8 * (n + 1)) as usize + 8 * (n + 2) > len {
            return None;
        }
        let base = b.as_ptr();
        Some(Self { _bytes: bundle.bytes.clone(), offsets: unsafe { base.add(8) } as *const u64, data: unsafe { base.add(8 * (n + 2)) }, n })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn get(&self, id: usize) -> Option<&str> {
        if id >= self.n {
            return None;
        }
        // The names were written from Strings, so they are UTF-8; the offsets are 8-aligned (the
        // section starts on 2 MiB and the count is 8 bytes).
        unsafe {
            let start = *self.offsets.add(id) as usize;
            let end = *self.offsets.add(id + 1) as usize;
            Some(std::str::from_utf8_unchecked(std::slice::from_raw_parts(self.data.add(start), end - start)))
        }
    }
}

/// Reference names: deserialised from a name file, or a [`NameTable`] in a bundle.
#[derive(Clone)]
pub enum Names {
    Owned(Vec<String>),
    Table(NameTable),
}

impl Names {
    #[inline]
    pub fn get(&self, id: usize) -> Option<&str> {
        match self {
            Names::Owned(v) => v.get(id).map(|s| s.as_str()),
            Names::Table(t) => t.get(id),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Names::Owned(v) => v.len(),
            Names::Table(t) => t.len(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_temp_dir_is_not_a_network_filesystem() {
        let path = std::env::temp_dir().join(format!("flexalign_fs_{}.bin", std::process::id()));
        std::fs::write(&path, b"x").unwrap();
        assert!(!on_network_fs(&File::open(&path).unwrap()));
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn names_and_sections_round_trip() {
        let path = std::env::temp_dir().join(format!("flexalign_bundle_{}.bin", std::process::id()));
        let path = path.to_string_lossy().to_string();
        let names: Vec<String> = vec!["dummy".into(), "chr1A".into(), "".into(), "NZ_CP0000.1".into()];
        write(
            &path,
            |w| { w.write_all(b"INDEXBYTES")?; Ok(10) },
            |w| { w.write_all(b"REFSEQ")?; Ok(6) },
            &|i| names[i].as_str(),
            names.len(),
        )
        .unwrap();
        for mmap in [false, true] {
            let bundle = Bundle::open(&path, mmap, 2).unwrap();
            let (offset, len) = bundle.section(INDEX).unwrap();
            assert_eq!(offset as u64 % ALIGN, 0);
            assert_eq!(&bundle.bytes.bytes()[offset..offset + len], b"INDEXBYTES");
            let (offset, len) = bundle.section(REFSEQ).unwrap();
            assert_eq!(&bundle.bytes.bytes()[offset..offset + len], b"REFSEQ");
            let table = NameTable::from_bundle(&bundle).unwrap();
            drop(bundle);
            assert_eq!(table.len(), names.len());
            for (i, name) in names.iter().enumerate() {
                assert_eq!(table.get(i), Some(name.as_str()));
            }
            assert_eq!(table.get(names.len()), None);
        }
        assert_eq!(Bundle::peek_index(&path).unwrap()[..10], *b"INDEXBYTES");
        std::fs::remove_file(&path).unwrap();
    }
}
