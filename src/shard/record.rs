//! Evidence records written by a shard pass and replayed by the rejoin.
//!
//! A shard pass stops after the range lookup (see ProjectShard.md §1) and writes what it found
//! for each (read, mate). It never builds a `Seed`: the payload is the raw `VCell` u64 verbatim,
//! because `Seed::from_flexmer(qpos, rpos, rval, dist)` is a pure function of the cell plus the
//! range's `qpos` and `min_dist`, both of which are per-range rather than per-seed.
//!
//! Layout of one group (one read+mate that had at least one range in this shard). All positions
//! are chunk-relative, so `chunk_id` is never stored -- it is the index slot the group lives in.
//!
//! ```text
//! idx_delta    varint   read index within the chunk, delta from the previous group
//! mate         u8       0 | 1
//! n_ranges     varint
//! range[n_ranges]:
//!     qpos_delta   varint   delta from the previous range in this group (qpos ascends, unique)
//!     n_positions  varint   == vrange.positions.len(), the global sort key for the rejoin
//!     flags_dist   u8       bit0: headered; bits1..2: reserved; bits3..7: min_dist (0..=16)
//!     n_cells      varint
//!     cells        n_cells x u64 little-endian raw VCell
//! ```
//!
//! `min_dist` needs 5 bits, not 4: `HeaderSeq::dist` is a hamming distance over F=16 flanking
//! bases, so it spans 0..=16 -- 17 values. It is not just a flag for `dist == 0`; it lands in
//! `Seed::mismatch` and feeds anchor scoring, so an all-mismatch flank (16) must not alias to 0.
//!
//! Typical size is ~20 B per range and ~39 B per group (1.8 ranges per mate per shard at N=10).

use std::io;

pub const FLAG_HEADERED: u8 = 0b0000_0001;
const DIST_SHIFT: u8 = 3;
const DIST_MAX: u8 = 0x1F;

/// One range's evidence: the sort key, the seed-construction metadata, and the raw cells.
///
/// `cells` are already flex-filtered by the shard (the filter is key-local once the threshold is
/// fixed -- ProjectShard.md §5). Ranges failing the filter are not emitted at all.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RangeRecord {
    pub qpos: u32,
    /// `vrange.positions.len()`. The rejoin orders ranges by this, ascending.
    pub n_positions: u32,
    /// True when the value block carried flank headers (`size > HEADER_THRESHOLD`). Drives
    /// `Seed::from_flexmer` vs `from_coremer`, and whether the range consumes the `--ranges`
    /// budget.
    pub headered: bool,
    /// Hamming distance of the read's flanks to the block's closest header. Meaningless when
    /// `!headered`.
    pub min_dist: u8,
    pub cells: Vec<u64>,
}

/// One (read, mate) that had at least one surviving range in this shard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupRecord {
    /// Read index within the chunk.
    pub idx: u32,
    pub mate: u8,
    pub ranges: Vec<RangeRecord>,
}

#[derive(Debug, thiserror::Error)]
pub enum RecordError {
    #[error("truncated record at byte {0}")]
    Truncated(usize),
    #[error("varint at byte {0} does not terminate within 5 bytes")]
    VarintOverflow(usize),
    #[error("group indices must ascend within a chunk, got {got} after {prev}")]
    NonAscending { prev: u32, got: u32 },
    #[error("min_dist {0} does not fit in 5 bits")]
    DistOverflow(u8),
}

impl From<RecordError> for io::Error {
    fn from(e: RecordError) -> Self {
        io::Error::new(io::ErrorKind::InvalidData, e)
    }
}

fn put_varint(out: &mut Vec<u8>, mut v: u32) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn get_varint(buf: &[u8], pos: &mut usize) -> Result<u32, RecordError> {
    let start = *pos;
    let mut v: u32 = 0;
    let mut shift = 0;
    loop {
        let b = *buf.get(*pos).ok_or(RecordError::Truncated(*pos))?;
        *pos += 1;
        v |= ((b & 0x7F) as u32) << shift;
        if b & 0x80 == 0 {
            return Ok(v);
        }
        shift += 7;
        if shift >= 35 {
            return Err(RecordError::VarintOverflow(start));
        }
    }
}

/// Appends groups to a chunk buffer, tracking the deltas.
///
/// Groups must be pushed in ascending `(idx, mate)` -- the rejoin relies on it to advance its
/// cursors monotonically, and nothing downstream re-sorts.
pub struct GroupWriter {
    buf: Vec<u8>,
    prev_idx: u32,
    started: bool,
    groups: u32,
}

impl GroupWriter {
    pub fn new() -> Self {
        Self { buf: Vec::new(), prev_idx: 0, started: false, groups: 0 }
    }

    pub fn with_capacity(cap: usize) -> Self {
        Self { buf: Vec::with_capacity(cap), prev_idx: 0, started: false, groups: 0 }
    }

    pub fn groups(&self) -> u32 {
        self.groups
    }

    pub fn is_empty(&self) -> bool {
        self.groups == 0
    }

    pub fn bytes(&self) -> &[u8] {
        &self.buf
    }

    /// Hands out the encoded chunk and resets for the next one.
    pub fn take(&mut self) -> Vec<u8> {
        self.prev_idx = 0;
        self.started = false;
        self.groups = 0;
        std::mem::take(&mut self.buf)
    }

    pub fn push(&mut self, group: &GroupRecord) -> Result<(), RecordError> {
        if self.started && group.idx < self.prev_idx {
            return Err(RecordError::NonAscending { prev: self.prev_idx, got: group.idx });
        }
        let idx_delta = if self.started { group.idx - self.prev_idx } else { group.idx };
        self.prev_idx = group.idx;
        self.started = true;
        self.groups += 1;

        put_varint(&mut self.buf, idx_delta);
        self.buf.push(group.mate);
        put_varint(&mut self.buf, group.ranges.len() as u32);

        let mut prev_qpos = 0u32;
        for r in &group.ranges {
            if r.min_dist > DIST_MAX {
                return Err(RecordError::DistOverflow(r.min_dist));
            }
            put_varint(&mut self.buf, r.qpos - prev_qpos);
            prev_qpos = r.qpos;

            put_varint(&mut self.buf, r.n_positions);

            let mut flags = if r.headered { FLAG_HEADERED } else { 0 };
            flags |= r.min_dist << DIST_SHIFT;
            self.buf.push(flags);

            put_varint(&mut self.buf, r.cells.len() as u32);
            for cell in &r.cells {
                self.buf.extend_from_slice(&cell.to_le_bytes());
            }
        }
        Ok(())
    }
}

impl Default for GroupWriter {
    fn default() -> Self {
        Self::new()
    }
}

/// A monotone cursor over one shard's chunk blob.
///
/// The rejoin holds N of these -- one per shard -- and walks the chunk's reads in order, taking
/// from each cursor while its head matches the current read. A read with no hit in this shard
/// simply leaves the cursor's head above it (ProjectShard.md §9). Nothing synchronises the
/// cursors: they are independent readers over byte ranges whose extents came from the index.
pub struct GroupReader<'a> {
    buf: &'a [u8],
    pos: usize,
    prev_idx: u32,
    started: bool,
    head: Option<(u32, u8)>,
}

impl<'a> GroupReader<'a> {
    pub fn new(buf: &'a [u8]) -> Result<Self, RecordError> {
        let mut r = Self { buf, pos: 0, prev_idx: 0, started: false, head: None };
        r.head = r.peek_key()?;
        Ok(r)
    }

    /// `(idx, mate)` of the next group, without consuming it.
    pub fn head(&self) -> Option<(u32, u8)> {
        self.head
    }

    pub fn exhausted(&self) -> bool {
        self.head.is_none()
    }

    fn peek_key(&mut self) -> Result<Option<(u32, u8)>, RecordError> {
        if self.pos >= self.buf.len() {
            return Ok(None);
        }
        let mut probe = self.pos;
        let delta = get_varint(self.buf, &mut probe)?;
        let idx = if self.started { self.prev_idx + delta } else { delta };
        let mate = *self.buf.get(probe).ok_or(RecordError::Truncated(probe))?;
        Ok(Some((idx, mate)))
    }

    /// Consumes and returns the next group.
    pub fn next_group(&mut self) -> Result<Option<GroupRecord>, RecordError> {
        if self.pos >= self.buf.len() {
            self.head = None;
            return Ok(None);
        }

        let delta = get_varint(self.buf, &mut self.pos)?;
        let idx = if self.started { self.prev_idx + delta } else { delta };
        self.prev_idx = idx;
        self.started = true;

        let mate = *self.buf.get(self.pos).ok_or(RecordError::Truncated(self.pos))?;
        self.pos += 1;

        let n_ranges = get_varint(self.buf, &mut self.pos)?;
        let mut ranges = Vec::with_capacity(n_ranges as usize);

        let mut qpos = 0u32;
        for _ in 0..n_ranges {
            qpos += get_varint(self.buf, &mut self.pos)?;
            let n_positions = get_varint(self.buf, &mut self.pos)?;

            let flags = *self.buf.get(self.pos).ok_or(RecordError::Truncated(self.pos))?;
            self.pos += 1;

            let n_cells = get_varint(self.buf, &mut self.pos)? as usize;
            let need = n_cells * 8;
            if self.pos + need > self.buf.len() {
                return Err(RecordError::Truncated(self.pos));
            }
            let mut cells = Vec::with_capacity(n_cells);
            for _ in 0..n_cells {
                let mut b = [0u8; 8];
                b.copy_from_slice(&self.buf[self.pos..self.pos + 8]);
                cells.push(u64::from_le_bytes(b));
                self.pos += 8;
            }

            ranges.push(RangeRecord {
                qpos,
                n_positions,
                headered: flags & FLAG_HEADERED != 0,
                min_dist: (flags >> DIST_SHIFT) & DIST_MAX,
                cells,
            });
        }

        self.head = self.peek_key()?;
        Ok(Some(GroupRecord { idx, mate, ranges }))
    }

    /// Takes every group belonging to `idx`, in mate order. Returns nothing when this shard had
    /// no hit for that read -- the common case for most (read, shard) pairs.
    pub fn take_read(&mut self, idx: u32, out: &mut Vec<GroupRecord>) -> Result<(), RecordError> {
        while let Some((head_idx, _)) = self.head {
            if head_idx != idx {
                break;
            }
            match self.next_group()? {
                Some(g) => out.push(g),
                None => break,
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn range(qpos: u32, n_positions: u32, headered: bool, min_dist: u8, cells: &[u64]) -> RangeRecord {
        RangeRecord { qpos, n_positions, headered, min_dist, cells: cells.to_vec() }
    }

    fn roundtrip(groups: &[GroupRecord]) -> Vec<GroupRecord> {
        let mut w = GroupWriter::new();
        for g in groups {
            w.push(g).expect("push");
        }
        let buf = w.take();
        let mut r = GroupReader::new(&buf).expect("reader");
        let mut out = Vec::new();
        while let Some(g) = r.next_group().expect("decode") {
            out.push(g);
        }
        out
    }

    #[test]
    fn roundtrips_a_headered_range() {
        let g = GroupRecord {
            idx: 7,
            mate: 0,
            ranges: vec![range(22, 5, true, 0, &[0xDEAD_BEEF_0000_0001, 2])],
        };
        assert_eq!(roundtrip(&[g.clone()]), vec![g]);
    }

    #[test]
    fn roundtrips_the_documented_example() {
        // ProjectShard.md: idx 7, two ranges, 43 bytes.
        let g = GroupRecord {
            idx: 7,
            mate: 0,
            ranges: vec![
                range(22, 5, true, 0, &[1, 2]),
                range(63, 2, false, 0, &[3, 4]),
            ],
        };
        let mut w = GroupWriter::new();
        w.push(&g).unwrap();
        assert_eq!(w.bytes().len(), 43, "documented size drifted");
        assert_eq!(roundtrip(&[g.clone()]), vec![g]);
    }

    #[test]
    fn roundtrips_max_min_dist() {
        // F=16 flanking bases, so min_dist spans 0..=16 -- 17 values, needing 5 bits. dist==16
        // (an all-mismatch flank) is a real value that reaches Seed::mismatch; it must not be
        // rejected and must not alias to 0.
        for d in 0..=16u8 {
            let g = GroupRecord { idx: 0, mate: 1, ranges: vec![range(0, 3, true, d, &[9])] };
            assert_eq!(roundtrip(&[g.clone()]), vec![g], "min_dist {} did not survive", d);
        }
    }

    #[test]
    fn min_dist_does_not_collide_with_flags() {
        // dist occupies bits3..7, so a dist that sets its low bit must not leak into
        // FLAG_HEADERED, and a headered range must not inflate dist.
        for d in 0..=16u8 {
            for headered in [false, true] {
                let g = GroupRecord {
                    idx: 1,
                    mate: 0,
                    ranges: vec![range(0, 3, headered, d, &[9])],
                };
                let out = roundtrip(&[g.clone()]);
                assert_eq!(out[0].ranges[0].min_dist, d, "dist {} headered {}", d, headered);
                assert_eq!(out[0].ranges[0].headered, headered, "dist {} headered {}", d, headered);
            }
        }
    }

    #[test]
    fn rejects_dist_beyond_five_bits() {
        let mut w = GroupWriter::new();
        let g = GroupRecord { idx: 0, mate: 0, ranges: vec![range(0, 3, true, 32, &[9])] };
        assert!(matches!(w.push(&g), Err(RecordError::DistOverflow(32))));
    }

    #[test]
    fn roundtrips_empty_ranges_and_empty_cells() {
        let g = GroupRecord { idx: 3, mate: 1, ranges: vec![] };
        assert_eq!(roundtrip(&[g.clone()]), vec![g]);

        // n_cells == 0 is how a discarded range would ride along if the recovery path (§4a)
        // is ever reinstated.
        let g = GroupRecord { idx: 4, mate: 0, ranges: vec![range(10, 200, true, 3, &[])] };
        assert_eq!(roundtrip(&[g.clone()]), vec![g]);
    }

    #[test]
    fn roundtrips_many_groups_with_large_deltas() {
        let groups: Vec<GroupRecord> = (0..64)
            .map(|i| GroupRecord {
                idx: i * 1000,
                mate: (i % 2) as u8,
                ranges: vec![range(i, 300 + i, i % 2 == 0, (i % 16) as u8, &[i as u64, 1 << 40])],
            })
            .collect();
        assert_eq!(roundtrip(&groups), groups);
    }

    #[test]
    fn rejects_descending_idx() {
        let mut w = GroupWriter::new();
        w.push(&GroupRecord { idx: 5, mate: 0, ranges: vec![] }).unwrap();
        let err = w.push(&GroupRecord { idx: 4, mate: 0, ranges: vec![] });
        assert!(matches!(err, Err(RecordError::NonAscending { prev: 5, got: 4 })));
    }

    #[test]
    fn head_exposes_next_key_without_consuming() {
        let groups = vec![
            GroupRecord { idx: 2, mate: 0, ranges: vec![range(1, 4, true, 1, &[7])] },
            GroupRecord { idx: 2, mate: 1, ranges: vec![range(1, 4, true, 1, &[8])] },
            GroupRecord { idx: 9, mate: 0, ranges: vec![] },
        ];
        let mut w = GroupWriter::new();
        for g in &groups {
            w.push(g).unwrap();
        }
        let buf = w.take();
        let mut r = GroupReader::new(&buf).unwrap();

        assert_eq!(r.head(), Some((2, 0)));
        assert_eq!(r.head(), Some((2, 0)), "head must not consume");
        r.next_group().unwrap();
        assert_eq!(r.head(), Some((2, 1)));
        r.next_group().unwrap();
        assert_eq!(r.head(), Some((9, 0)));
        r.next_group().unwrap();
        assert_eq!(r.head(), None);
        assert!(r.exhausted());
    }

    #[test]
    fn take_read_gathers_both_mates_and_skips_absent_reads() {
        let groups = vec![
            GroupRecord { idx: 2, mate: 0, ranges: vec![range(1, 4, true, 1, &[7])] },
            GroupRecord { idx: 2, mate: 1, ranges: vec![range(5, 4, true, 1, &[8])] },
            GroupRecord { idx: 9, mate: 0, ranges: vec![range(3, 2, false, 0, &[9])] },
        ];
        let mut w = GroupWriter::new();
        for g in &groups {
            w.push(g).unwrap();
        }
        let buf = w.take();
        let mut r = GroupReader::new(&buf).unwrap();

        let mut out = Vec::new();
        r.take_read(0, &mut out).unwrap();
        assert!(out.is_empty(), "read 0 has no evidence in this shard");

        r.take_read(2, &mut out).unwrap();
        assert_eq!(out.len(), 2, "both mates of read 2");
        assert_eq!(out[0].mate, 0);
        assert_eq!(out[1].mate, 1);

        out.clear();
        for idx in 3..9 {
            r.take_read(idx, &mut out).unwrap();
        }
        assert!(out.is_empty(), "reads 3..9 absent");

        r.take_read(9, &mut out).unwrap();
        assert_eq!(out.len(), 1);
        assert!(r.exhausted());
    }

    #[test]
    fn empty_buffer_reads_as_no_groups() {
        let mut r = GroupReader::new(&[]).unwrap();
        assert!(r.exhausted());
        assert_eq!(r.next_group().unwrap(), None);
    }

    #[test]
    fn truncation_is_an_error_not_a_panic() {
        let mut w = GroupWriter::new();
        w.push(&GroupRecord { idx: 1, mate: 0, ranges: vec![range(4, 9, true, 2, &[1, 2, 3])] })
            .unwrap();
        let full = w.take();
        for cut in 1..full.len() {
            let mut r = match GroupReader::new(&full[..cut]) {
                Ok(r) => r,
                Err(_) => continue,
            };
            // Must not panic; either decodes or errors.
            let _ = r.next_group();
        }
    }
}
