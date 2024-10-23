use std::fmt::Display;


pub struct Flag(u16);

impl Flag {
    pub fn new() -> Self {
        Self(0)
    }

    pub fn paired_end(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x1u16,   // Use |= to set the bit
            false => self.0 &= !0x1u16, // Use &= with a negated mask to clear the bit
        };
        self
    }

    pub fn both_aligned(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x2u16,
            false => self.0 &= !0x2u16,
        };
        self
    }

    pub fn read1_mapped(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x4u16,
            false => self.0 &= !0x4u16,
        };
        self
    }

    pub fn read2_mapped(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x8u16,
            false => self.0 &= !0x8u16,
        };
        self
    }

    pub fn read1_rc(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x16u16,
            false => self.0 &= !0x16u16,
        };
        self
    }

    pub fn read2_rc(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x32u16,
            false => self.0 &= !0x32u16,
        };
        self
    }

    pub fn read1(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x64u16,
            false => self.0 &= !0x64u16,
        };
        self
    }

    pub fn read2(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x128u16,
            false => self.0 &= !0x128u16,
        };
        self
    }

    pub fn not_primary(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x256u16,
            false => self.0 &= !0x256u16,
        };
        self
    }

    pub fn alignment_failed_qc(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x512u16,
            false => self.0 &= !0x512u16,
        };
        self
    }

    pub fn duplicate(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x1024u16,
            false => self.0 &= !0x1024u16,
        };
        self
    }

    pub fn supplementary(&mut self, on: bool) -> &mut Self {
        match on {
            true => self.0 |= 0x2048u16,
            false => self.0 &= !0x2048u16,
        };
        self
    }

    pub fn is_paired_end(&self) -> bool {
        (self.0 & 0x1u16) != 0
    }

    pub fn is_both_aligned(&self) -> bool {
        (self.0 & 0x2u16) != 0
    }

    pub fn is_read1_mapped(&self) -> bool {
        (self.0 & 0x4u16) != 0
    }

    pub fn is_read2_mapped(&self) -> bool {
        (self.0 & 0x8u16) != 0
    }

    pub fn is_read1_rc(&self) -> bool {
        (self.0 & 0x16u16) != 0
    }

    pub fn is_read2_rc(&self) -> bool {
        (self.0 & 0x32u16) != 0
    }

    pub fn is_read1(&self) -> bool {
        (self.0 & 0x64u16) != 0
    }

    pub fn is_read2(&self) -> bool {
        (self.0 & 0x128u16) != 0
    }

    pub fn is_not_primary(&self) -> bool {
        (self.0 & 0x256u16) != 0
    }

    pub fn is_alignmend_failed_qc(&self) -> bool {
        (self.0 & 0x512u16) != 0
    }

    pub fn is_duplicate(&self) -> bool {
        (self.0 & 0x1024u16) != 0
    }

    pub fn is_supplementary(&self) -> bool {
        (self.0 & 0x2048u16) != 0
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
}

impl Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}
