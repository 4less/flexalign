//! A reference sequence as alignment sees it: its bytes, or a sequence of a 2-bit store decoded a
//! window at a time.
//!
//! A whole genome at a byte a base is the largest thing the sharded join holds (3.1 GB of human);
//! at 2 bits a base it is a quarter. But decoding a whole chromosome per read would cost far more
//! than the read, so a long packed sequence hands out windows instead: `&view[a..b]` decodes those
//! bases alone into this thread's scratch (valid until the next `refseq2::reset_scratch`, which
//! runs when a read starts -- the same contract as a short sequence's `RefSeq2::get`). Callers
//! index it as they indexed the byte slice: `len()`, `[a..b]`, `[i]`. Not `[a..]` or `[..b]`: on a
//! chromosome those decode megabases, so a caller says where its window ends.

use std::ops::{Index, Range};

use super::refseq2::RefSeq2;

#[derive(Clone, Copy)]
pub enum RefView<'a> {
    Bytes(&'a [u8]),
    /// Sequence `id` of a 2-bit store, `len` bases.
    Packed { store: &'a RefSeq2, id: usize, len: usize },
}

/// Every byte value, so `view[i]` can hand out a reference to a decoded base.
static BYTES: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    t
};

impl<'a> RefView<'a> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            RefView::Bytes(b) => b.len(),
            RefView::Packed { len, .. } => *len,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bases `[r.start, r.end)`; panics past the end, as a slice does.
    #[inline]
    pub fn window(&self, r: Range<usize>) -> &'a [u8] {
        match *self {
            RefView::Bytes(b) => &b[r],
            RefView::Packed { store, id, len } => {
                assert!(r.start <= r.end && r.end <= len, "range {r:?} out of a reference of {len}");
                store.get_range(id, r.start, r.end)
            }
        }
    }

    /// Asks for the cache line holding base `i` (of the bytes, or of the packed bases).
    #[inline]
    pub fn prefetch(&self, i: usize) -> *const u8 {
        match *self {
            RefView::Bytes(b) => b[i..].as_ptr(),
            RefView::Packed { store, id, .. } => store.packed_ptr(id, i),
        }
    }

    /// Base `i`.
    #[inline]
    pub fn at(&self, i: usize) -> u8 {
        match *self {
            RefView::Bytes(b) => b[i],
            RefView::Packed { .. } => self.window(i..i + 1)[0],
        }
    }
}

impl Index<Range<usize>> for RefView<'_> {
    type Output = [u8];
    #[inline]
    fn index(&self, r: Range<usize>) -> &[u8] {
        self.window(r)
    }
}

impl Index<usize> for RefView<'_> {
    type Output = u8;
    #[inline]
    fn index(&self, i: usize) -> &u8 {
        &BYTES[self.at(i) as usize]
    }
}
