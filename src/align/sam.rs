use std::fmt::Display;


pub struct Flag(u16);

impl Flag {
    pub fn new() -> Self {
        Self(0)
    }

    /// The raw bitfield, for column 2 of a SAM record.
    pub fn value(&self) -> u16 {
        self.0
    }

    pub fn paired_end(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x1u16,   // Use |= to set the bit
            false => self.0 &= !0x1u16, // Use &= with a negated mask to clear the bit
        };
        self
    }

    pub fn proper_pair(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x2u16,
            false => self.0 &= !0x2u16,
        };
        self
    }

    pub fn unmapped(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x4u16,
            false => self.0 &= !0x4u16,
        };
        self
    }

    pub fn mate_unmapped(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x8u16,
            false => self.0 &= !0x8u16,
        };
        self
    }

    pub fn reverse(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x10u16,
            false => self.0 &= !0x10u16,
        };
        self
    }

    pub fn mate_reverse(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x20u16,
            false => self.0 &= !0x20u16,
        };
        self
    }

    pub fn first_in_pair(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x40u16,
            false => self.0 &= !0x40u16,
        };
        self
    }

    pub fn last_in_pair(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x80u16,
            false => self.0 &= !0x80u16,
        };
        self
    }

    pub fn not_primary(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x100u16,
            false => self.0 &= !0x100u16,
        };
        self
    }

    pub fn failed_qc(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x200u16,
            false => self.0 &= !0x200u16,
        };
        self
    }

    pub fn duplicate(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x400u16,
            false => self.0 &= !0x400u16,
        };
        self
    }

    pub fn supplementary(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x800u16,
            false => self.0 &= !0x800u16,
        };
        self
    }

    pub fn is_paired_end(&self) -> bool {
        (self.0 & 0x1u16) != 0
    }

    pub fn is_proper_pair(&self) -> bool {
        (self.0 & 0x2u16) != 0
    }

    pub fn is_unmapped(&self) -> bool {
        (self.0 & 0x4u16) != 0
    }

    pub fn is_mate_unmapped(&self) -> bool {
        (self.0 & 0x8u16) != 0
    }

    pub fn is_reverse(&self) -> bool {
        (self.0 & 0x10u16) != 0
    }

    pub fn is_mate_reverse(&self) -> bool {
        (self.0 & 0x20u16) != 0
    }

    pub fn is_first_in_pair(&self) -> bool {
        (self.0 & 0x40u16) != 0
    }

    pub fn is_last_in_pair(&self) -> bool {
        (self.0 & 0x80u16) != 0
    }

    pub fn is_not_primary(&self) -> bool {
        (self.0 & 0x100u16) != 0
    }

    pub fn is_failed_qc(&self) -> bool {
        (self.0 & 0x200u16) != 0
    }

    pub fn is_duplicate(&self) -> bool {
        (self.0 & 0x400u16) != 0
    }

    pub fn is_supplementary(&self) -> bool {
        (self.0 & 0x800u16) != 0
    }
}

type CigarOp = u8;


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cigar(pub Vec<CigarOp>);
pub struct CigarRef<'a>(pub &'a [CigarOp]);

impl AsRef::<Cigar> for Cigar {
    fn as_ref(&self) -> &Cigar {
        &self
    }
}

// \*|([0-9]+[MIDNSHP=X])+
impl Cigar {
    pub fn add_softclip(&mut self, n: usize) {
        self.0.extend(std::iter::repeat(b'S').take(n));
    }

    pub fn add_matches(&mut self, n: usize) {
        self.0.extend(std::iter::repeat(b'M').take(n));
    }

    pub fn count_leading_chars(&self, c: u8) -> usize {
        self.0.iter()
            .take_while(|&ch| *ch == c)
            .count()
    }

    pub fn count_trailing_chars(&self, c: u8) -> usize {
        self.0.iter()
            .rev()
            .take_while(|&ch| *ch == c)
            .count()
    }

    pub fn valid(&self) -> bool {
        true
    }

    pub fn new() -> Self {
        Self { 0: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self { 0: Vec::with_capacity(capacity) }
    }

    /// How many query bases this CIGAR consumes.
    ///
    /// Note the operation convention in use here is WFA2's, not SAM's: `D` is a gap in the
    /// *reference* and therefore consumes query, while `I` is a gap in the *query* and consumes
    /// reference. That is the exact opposite of SAM's meaning for the same two letters, which is why
    /// [`to_sam`](Self::to_sam) swaps them. `print_alignment` and `is_alignment_valid` in
    /// `align::common` walk the CIGAR under the same convention.
    pub fn query_consumed(&self) -> usize {
        self.0.iter().filter(|c| matches!(**c, b'M' | b'X' | b'D' | b'S')).count()
    }

    /// How many reference bases this CIGAR consumes.
    pub fn reference_consumed(&self) -> usize {
        self.0.iter().filter(|c| matches!(**c, b'M' | b'X' | b'I')).count()
    }

    /// Renders the expanded, WFA-convention CIGAR as a run-length-encoded SAM CIGAR string.
    ///
    /// Two translations happen here:
    /// - `I` and `D` are swapped, per the convention note on [`query_consumed`](Self::query_consumed).
    /// - `M` and `X` are both emitted as SAM `M` ("alignment match", which the spec defines as
    ///   either a sequence match or mismatch). Keeping `=`/`X` would be more precise but is far less
    ///   widely consumed; the mismatch count is reported in the `NM` tag instead.
    ///
    /// Returns `None` if the CIGAR consumes more query bases than the read has, which would make an
    /// invalid record -- callers emit nothing rather than corrupt output. If it consumes fewer, the
    /// remainder is soft-clipped so that query length and CIGAR agree, as SAM requires.
    pub fn to_sam(&self, query_length: usize) -> Option<String> {
        let consumed = self.query_consumed();
        if consumed > query_length {
            return None;
        }

        let mut out = String::with_capacity(16);
        let mut run_op = 0u8;
        let mut run_len = 0usize;

        let mut flush = |op: u8, len: usize, out: &mut String| {
            if len == 0 {
                return;
            }
            let sam_op = match op {
                b'M' | b'X' => 'M',
                b'D' => 'I', // gap in reference == insertion, in SAM's vocabulary
                b'I' => 'D', // gap in query == deletion
                b'S' => 'S',
                other => panic!("Unknown cigar op {}", other as char),
            };
            out.push_str(&len.to_string());
            out.push(sam_op);
        };

        for op in self.0.iter() {
            // M and X collapse into one SAM run, so treat X as M for run-breaking purposes.
            let op = if *op == b'X' { b'M' } else { *op };
            if op == run_op {
                run_len += 1;
            } else {
                flush(run_op, run_len, &mut out);
                run_op = op;
                run_len = 1;
            }
        }
        flush(run_op, run_len, &mut out);

        if consumed < query_length {
            out.push_str(&(query_length - consumed).to_string());
            out.push('S');
        }

        if out.is_empty() { None } else { Some(out) }
    }

    /// Edit distance (SAM's `NM` tag): mismatches plus gap bases.
    pub fn edit_distance(&self) -> u32 {
        self.0.iter().filter(|c| matches!(**c, b'X' | b'I' | b'D')).count() as u32
    }
}

impl Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(s: &str) -> Cigar {
        Cigar(s.as_bytes().to_vec())
    }

    #[test]
    fn to_sam_run_length_encodes_and_merges_matches_with_mismatches() {
        // M and X are both SAM 'M', so they merge into a single run.
        assert_eq!(cigar("MMMXMMM").to_sam(7).unwrap(), "7M");
        assert_eq!(cigar("MMMMM").to_sam(5).unwrap(), "5M");
    }

    #[test]
    fn to_sam_swaps_the_indel_convention() {
        // Internally D consumes query (a gap in the reference) == SAM insertion, and I consumes
        // reference == SAM deletion. Getting this backwards silently shifts every downstream
        // coordinate, so pin it.
        assert_eq!(cigar("MMDDMM").to_sam(6).unwrap(), "2M2I2M");
        assert_eq!(cigar("MMIIMM").to_sam(4).unwrap(), "2M2D2M");
    }

    #[test]
    fn to_sam_pads_short_cigars_with_soft_clips() {
        // SAM requires the CIGAR to account for every base of SEQ.
        assert_eq!(cigar("MMMM").to_sam(10).unwrap(), "4M6S");
        assert_eq!(cigar("SSMMMM").to_sam(6).unwrap(), "2S4M");
    }

    #[test]
    fn to_sam_refuses_a_cigar_longer_than_the_read() {
        assert!(cigar("MMMMMMMM").to_sam(4).is_none());
    }

    #[test]
    fn consumed_counts_follow_the_internal_convention() {
        let c = cigar("MMXDDIIS");
        assert_eq!(c.query_consumed(), 6);     // M M X D D S
        assert_eq!(c.reference_consumed(), 5); // M M X I I
        assert_eq!(c.edit_distance(), 5);      // X D D I I
    }

    #[test]
    fn flag_bits_are_the_sam_spec_values() {
        let mut f = Flag::new();
        f.paired_end(true).proper_pair(true).reverse(true).first_in_pair(true);
        // 1 + 2 + 16 + 64: the flag samtools shows for "read1, reverse strand, properly paired".
        assert_eq!(f.value(), 0x1 | 0x2 | 0x10 | 0x40);
        assert_eq!(f.value(), 83);
        assert!(f.is_reverse() && f.is_first_in_pair() && !f.is_mate_reverse());

        // Bits must not collide: setting then clearing one leaves the others untouched.
        f.reverse(false);
        assert_eq!(f.value(), 0x1 | 0x2 | 0x40);
    }
}
