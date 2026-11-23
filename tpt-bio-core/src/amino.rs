//! Amino acids and codon translation.

use crate::nucleotide::Base;

/// A standard amino acid, or a stop codon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AminoAcid {
    /// Alanine (A).
    Ala,
    /// Arginine (R).
    Arg,
    /// Asparagine (N).
    Asn,
    /// Aspartic acid (D).
    Asp,
    /// Cysteine (C).
    Cys,
    /// Glutamine (Q).
    Gln,
    /// Glutamic acid (E).
    Glu,
    /// Glycine (G).
    Gly,
    /// Histidine (H).
    His,
    /// Isoleucine (I).
    Ile,
    /// Leucine (L).
    Leu,
    /// Lysine (K).
    Lys,
    /// Methionine (M).
    Met,
    /// Phenylalanine (F).
    Phe,
    /// Proline (P).
    Pro,
    /// Serine (S).
    Ser,
    /// Threonine (T).
    Thr,
    /// Tryptophan (W).
    Trp,
    /// Tyrosine (Y).
    Tyr,
    /// Valine (V).
    Val,
    /// Stop codon (*).
    Stop,
}

impl AminoAcid {
    /// Parses a one-letter code. `*` is [`AminoAcid::Stop`]. Case-insensitive.
    #[must_use]
    pub const fn from_one_letter(b: u8) -> Option<Self> {
        Some(match b.to_ascii_uppercase() {
            b'A' => Self::Ala,
            b'R' => Self::Arg,
            b'N' => Self::Asn,
            b'D' => Self::Asp,
            b'C' => Self::Cys,
            b'Q' => Self::Gln,
            b'E' => Self::Glu,
            b'G' => Self::Gly,
            b'H' => Self::His,
            b'I' => Self::Ile,
            b'L' => Self::Leu,
            b'K' => Self::Lys,
            b'M' => Self::Met,
            b'F' => Self::Phe,
            b'P' => Self::Pro,
            b'S' => Self::Ser,
            b'T' => Self::Thr,
            b'W' => Self::Trp,
            b'Y' => Self::Tyr,
            b'V' => Self::Val,
            b'*' => Self::Stop,
            _ => return None,
        })
    }

    /// The one-letter code (`*` for stop).
    #[must_use]
    pub const fn to_one_letter(self) -> u8 {
        b"ARNDCQEGHILKMFPSTWYV*"[self as usize]
    }
}

/// Maps a codon to an amino acid. Built from a 64-entry table in TCAG order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
///
/// # Examples
///
/// ```
/// use tpt_bio_core::{AminoAcid, CodonTable};
/// use tpt_bio_core::nucleotide::Base;
///
/// let table = CodonTable::standard();
/// assert_eq!(table.translate([Base::A, Base::T, Base::G]), AminoAcid::Met);
/// assert_eq!(table.translate_ascii(b"TAA"), Some(AminoAcid::Stop));
/// ```
pub struct CodonTable {
    table: [u8; 64],
}

/// NCBI translation table 1 (the standard genetic code), indexed by codons in
/// TCAG order: index = 16·i1 + 4·i2 + i3, with T=0, C=1, A=2, G=3.
const STANDARD_TCAG: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

impl CodonTable {
    /// The standard genetic code (NCBI table 1).
    #[must_use]
    pub const fn standard() -> Self {
        let mut table = [0u8; 64];
        let mut i = 0;
        while i < 64 {
            table[i] = STANDARD_TCAG[i];
            i += 1;
        }
        Self { table }
    }

    /// Translates a codon of unambiguous bases.
    #[must_use]
    pub fn translate(&self, codon: [Base; 3]) -> AminoAcid {
        let idx = tcag_index(codon[0]) * 16 + tcag_index(codon[1]) * 4 + tcag_index(codon[2]);
        // The table only contains valid one-letter codes.
        AminoAcid::from_one_letter(self.table[idx]).unwrap_or(AminoAcid::Stop)
    }

    /// Translates an ASCII codon. Returns `None` if any byte is not an unambiguous base.
    #[must_use]
    pub fn translate_ascii(&self, codon: &[u8; 3]) -> Option<AminoAcid> {
        let b =
            [Base::from_ascii(codon[0])?, Base::from_ascii(codon[1])?, Base::from_ascii(codon[2])?];
        Some(self.translate(b))
    }
}

impl Default for CodonTable {
    fn default() -> Self {
        Self::standard()
    }
}

const fn tcag_index(b: Base) -> usize {
    match b {
        Base::T => 0,
        Base::C => 1,
        Base::A => 2,
        Base::G => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_table_known_codons() {
        let t = CodonTable::standard();
        assert_eq!(t.translate_ascii(b"ATG"), Some(AminoAcid::Met));
        assert_eq!(t.translate_ascii(b"TGG"), Some(AminoAcid::Trp));
        assert_eq!(t.translate_ascii(b"TAA"), Some(AminoAcid::Stop));
        assert_eq!(t.translate_ascii(b"TAG"), Some(AminoAcid::Stop));
        assert_eq!(t.translate_ascii(b"TGA"), Some(AminoAcid::Stop));
        assert_eq!(t.translate_ascii(b"GCT"), Some(AminoAcid::Ala));
        assert_eq!(t.translate_ascii(b"aaa"), Some(AminoAcid::Lys));
        assert_eq!(t.translate_ascii(b"ANA"), None);
    }

    #[test]
    fn every_codon_maps_to_a_known_amino_acid() {
        let t = CodonTable::standard();
        let bases = [Base::T, Base::C, Base::A, Base::G];
        let mut stops = 0;
        for &a in &bases {
            for &b in &bases {
                for &c in &bases {
                    if t.translate([a, b, c]) == AminoAcid::Stop {
                        stops += 1;
                    }
                }
            }
        }
        assert_eq!(stops, 3);
    }

    #[test]
    fn one_letter_roundtrip() {
        for b in b"ARNDCQEGHILKMFPSTWYV*" {
            let aa = AminoAcid::from_one_letter(*b).unwrap();
            assert_eq!(aa.to_one_letter(), *b);
        }
        assert_eq!(AminoAcid::from_one_letter(b'X'), None);
    }
}
