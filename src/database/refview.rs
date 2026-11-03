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

    /// Mismatches of `ascii` (the read) against bases `[from, from + ascii.len())`, stopping past
    /// `bound`; `packed` is the same read at 2 bits a base (see `RefSeq2::mismatches_packed`), used
    /// on a packed sequence so nothing is decoded.
    #[inline]
    pub fn mismatches(&self, from: usize, ascii: &[u8], packed: &[u64], bound: usize) -> usize {
        if let RefView::Packed { store, id, .. } = *self {
            if let Some(m) = store.mismatches_packed(id, from, packed, ascii.len(), bound) {
                // FA_TIE_CHECK: the decoded comparison must agree (both at or under the bound and
                // equal, or both over it), ignoring windows with bases outside ACGT.
                static CHECK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
                if *CHECK.get_or_init(|| std::env::var_os("FA_TIE_CHECK").is_some()) {
                    let w = self.window(from..from + ascii.len());
                    if w.iter().all(|b| b"ACGT".contains(b)) && ascii.iter().all(|b| b"ACGT".contains(b)) {
                        let full = ascii.iter().zip(w).filter(|(a, b)| a != b).count();
                        let ok = if m <= bound { full == m } else { full > bound };
                        assert!(ok, "packed mismatches {m} vs decoded {full} (bound {bound}) at {id}:{from}");
                    }
                }
                return m;
            }
        }
        let w = self.window(from..from + ascii.len());
        let mut c = 0;
        for (a, b) in ascii.iter().zip(w) {
            if a != b {
                c += 1;
                if c > bound {
                    break;
                }
            }
        }
        c
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
