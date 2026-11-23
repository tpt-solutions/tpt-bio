//! Nucleotide types: the four unambiguous bases and the full IUPAC alphabet.
//!
//! [`Iupac`] uses the 4-bit code of the BAM specification (`=ACMGRSVTWYHKDBN`),
//! so [`Iupac::to_nibble`] is the value stored in BAM `SEQ` fields.

use alloc::vec::Vec;

use crate::error::InvalidSymbol;

/// One of the four unambiguous DNA bases, with its 2-bit code (A=0, C=1, G=2, T=3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
///
/// # Examples
///
/// ```
/// use tpt_bio_core::nucleotide::Base;
///
/// assert_eq!(Base::from_ascii(b'g'), Some(Base::G));
/// assert_eq!(Base::G.complement(), Base::C);
/// assert_eq!(Base::T.to_2bit(), 3);
/// ```
pub enum Base {
    /// Adenine.
    A = 0,
    /// Cytosine.
    C = 1,
    /// Guanine.
    G = 2,
    /// Thymine (uracil in RNA input is read as thymine).
    T = 3,
}

impl Base {
    /// Parses an ASCII base, case-insensitively. `U` is read as `T`, for RNA input.
    #[must_use]
    pub const fn from_ascii(b: u8) -> Option<Self> {
        match b {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' | b'U' | b'u' => Some(Self::T),
            _ => None,
        }
    }

    /// Uppercase ASCII representation.
    #[must_use]
    pub const fn to_ascii(self) -> u8 {
        b"ACGT"[self as usize]
    }

    /// Watson-Crick complement.
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    /// The 2-bit code (0..=3).
    #[must_use]
    pub const fn to_2bit(self) -> u8 {
        self as u8
    }

    /// Decodes a 2-bit code. Returns `None` for values above 3.
    #[must_use]
    pub const fn from_2bit(v: u8) -> Option<Self> {
        match v {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            _ => None,
        }
    }
}

/// Full IUPAC nucleotide code, with the BAM 4-bit encoding as the discriminant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
///
/// # Examples
///
/// ```
/// use tpt_bio_core::nucleotide::Iupac;
///
/// let n = Iupac::from_ascii(b'N').unwrap();
/// assert_eq!(n.to_ascii(), b'N');
/// assert_eq!(n.to_base(), None);
/// assert_eq!(Iupac::from_ascii(b'R').unwrap().complement().to_ascii(), b'Y');
/// ```
pub enum Iupac {
    /// `=`: identical to the reference (BAM only).
    Eq = 0,
    /// A.
    A = 1,
    /// C.
    C = 2,
    /// M: A or C.
    M = 3,
    /// G.
    G = 4,
    /// R: A or G.
    R = 5,
    /// S: C or G.
    S = 6,
    /// V: A, C or G.
    V = 7,
    /// T (or U).
    T = 8,
    /// W: A or T.
    W = 9,
    /// Y: C or T.
    Y = 10,
    /// H: A, C or T.
    H = 11,
    /// K: G or T.
    K = 12,
    /// D: A, G or T.
    D = 13,
    /// B: C, G or T.
    B = 14,
    /// N: any base.
    N = 15,
}

const NIBBLE_TO_ASCII: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

impl Iupac {
    /// Parses an ASCII IUPAC code, case-insensitively. `U` is read as `T`.
    #[must_use]
    pub const fn from_ascii(b: u8) -> Option<Self> {
        let upper = b.to_ascii_uppercase();
        let upper = if upper == b'U' { b'T' } else { upper };
        let mut i = 1;
        while i < 16 {
            if NIBBLE_TO_ASCII[i] == upper {
                return Self::from_nibble(i as u8);
            }
            i += 1;
        }
        None
    }

    /// Uppercase ASCII representation.
    #[must_use]
    pub const fn to_ascii(self) -> u8 {
        NIBBLE_TO_ASCII[self as usize]
    }

    /// The 4-bit BAM code (0..=15).
    #[must_use]
    pub const fn to_nibble(self) -> u8 {
        self as u8
    }

    /// Decodes a 4-bit BAM code. Returns `None` for values above 15.
    #[must_use]
    pub const fn from_nibble(v: u8) -> Option<Self> {
        Some(match v {
            0 => Self::Eq,
            1 => Self::A,
            2 => Self::C,
            3 => Self::M,
            4 => Self::G,
            5 => Self::R,
            6 => Self::S,
            7 => Self::V,
            8 => Self::T,
            9 => Self::W,
            10 => Self::Y,
            11 => Self::H,
            12 => Self::K,
            13 => Self::D,
            14 => Self::B,
            15 => Self::N,
            _ => return None,
        })
    }

    /// Complement of the code; ambiguity codes complement to their mirror set.
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::Eq => Self::Eq,
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::M => Self::K,
            Self::K => Self::M,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::S => Self::S,
            Self::W => Self::W,
            Self::V => Self::B,
            Self::B => Self::V,
            Self::H => Self::D,
            Self::D => Self::H,
            Self::N => Self::N,
        }
    }

    /// Returns the base if this code is one of A, C, G or T.
    #[must_use]
    pub const fn to_base(self) -> Option<Base> {
        match self {
            Self::A => Some(Base::A),
            Self::C => Some(Base::C),
            Self::G => Some(Base::G),
            Self::T => Some(Base::T),
            _ => None,
        }
    }
}

impl From<Base> for Iupac {
    fn from(b: Base) -> Self {
        // A=1, C=2, G=4, T=8 in the BAM nibble encoding.
        match b {
            Base::A => Self::A,
            Base::C => Self::C,
            Base::G => Self::G,
            Base::T => Self::T,
        }
    }
}

/// Complement of one ASCII IUPAC byte, preserving case. `None` if the byte is not IUPAC.
#[must_use]
pub const fn complement_ascii(b: u8) -> Option<u8> {
    match Iupac::from_ascii(b) {
        Some(code) => {
            let c = code.complement().to_ascii();
            Some(if b.is_ascii_lowercase() { c.to_ascii_lowercase() } else { c })
        }
        None => None,
    }
}

/// Reverse complement of an ASCII IUPAC sequence. Case is preserved.
///
/// # Errors
/// Returns [`InvalidSymbol`] for the first byte that is not an IUPAC code.
pub fn reverse_complement(seq: &[u8]) -> Result<Vec<u8>, InvalidSymbol> {
    let mut out = Vec::with_capacity(seq.len());
    for (i, &b) in seq.iter().enumerate().rev() {
        let c = complement_ascii(b).ok_or(InvalidSymbol { position: i, byte: b })?;
        out.push(c);
    }
    Ok(out)
}

/// In-place reverse complement of an ASCII IUPAC sequence. Case is preserved.
///
/// The slice is left unchanged if an invalid byte is found.
///
/// # Errors
/// Returns [`InvalidSymbol`] for the first byte that is not an IUPAC code.
pub fn reverse_complement_in_place(seq: &mut [u8]) -> Result<(), InvalidSymbol> {
    if let Some((position, &byte)) =
        seq.iter().enumerate().find(|(_, b)| complement_ascii(**b).is_none())
    {
        return Err(InvalidSymbol { position, byte });
    }
    let n = seq.len();
    for i in 0..n / 2 {
        seq.swap(i, n - 1 - i);
    }
    for b in seq.iter_mut() {
        // Every byte was validated above.
        *b = complement_ascii(*b).unwrap_or(*b);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_roundtrip() {
        for b in [b'A', b'C', b'G', b'T'] {
            let base = Base::from_ascii(b).unwrap();
            assert_eq!(base.to_ascii(), b);
            assert_eq!(Base::from_2bit(base.to_2bit()), Some(base));
        }
        assert_eq!(Base::from_ascii(b'u'), Some(Base::T));
        assert_eq!(Base::from_ascii(b'N'), None);
        assert_eq!(Base::from_2bit(4), None);
    }

    #[test]
    fn iupac_nibble_matches_bam_spec() {
        let encoded: Vec<u8> = (0..16).map(|v| Iupac::from_nibble(v).unwrap().to_ascii()).collect();
        assert_eq!(&encoded, b"=ACMGRSVTWYHKDBN");
        for v in 0..16 {
            assert_eq!(Iupac::from_nibble(v).unwrap().to_nibble(), v);
        }
        assert_eq!(Iupac::from_nibble(16), None);
    }

    #[test]
    fn base_to_iupac_nibble_agrees() {
        for b in [Base::A, Base::C, Base::G, Base::T] {
            assert_eq!(Iupac::from(b).to_ascii(), b.to_ascii());
            assert_eq!(Iupac::from(b).to_base(), Some(b));
        }
    }

    #[test]
    fn iupac_complement_is_involution() {
        for v in 0..16 {
            let c = Iupac::from_nibble(v).unwrap();
            assert_eq!(c.complement().complement(), c);
        }
    }

    #[test]
    fn reverse_complement_preserves_case_and_ambiguity() {
        assert_eq!(reverse_complement(b"ACGTN").unwrap(), b"NACGT");
        assert_eq!(reverse_complement(b"acgRYm").unwrap(), b"kRYcgt");
        assert_eq!(reverse_complement(b"").unwrap(), b"");
    }

    #[test]
    fn reverse_complement_rejects_invalid() {
        let err = reverse_complement(b"AC?T").unwrap_err();
        assert_eq!(err, InvalidSymbol { position: 2, byte: b'?' });
    }

    #[test]
    fn in_place_matches_allocating_version() {
        let mut buf = *b"GATTACAnn";
        reverse_complement_in_place(&mut buf).unwrap();
        assert_eq!(&buf[..], &reverse_complement(b"GATTACAnn").unwrap()[..]);

        let mut bad = *b"AC!";
        assert!(reverse_complement_in_place(&mut bad).is_err());
        assert_eq!(&bad, b"AC!");
    }
}
