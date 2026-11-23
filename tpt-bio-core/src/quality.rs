//! Phred quality scores.
//!
//! Qualities are stored as the integer Q, where the error probability is
//! `10^(-Q/10)`. The FASTQ/SAM text encoding is Sanger Phred+33.

use crate::error::InvalidSymbol;

/// Highest Q representable in Phred+33 text (`~`, ASCII 126).
pub const MAX_Q: u8 = 93;

/// A Phred-scaled quality score in `0..=93`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
///
/// # Examples
///
/// ```
/// use tpt_bio_core::quality::Phred;
///
/// // 'I' is ASCII 73, so Phred+33 gives Q40.
/// let q = Phred::from_phred33(b'I').unwrap();
/// assert_eq!(q.value(), 40);
/// assert!((q.error_probability() - 1e-4).abs() < 1e-12);
/// ```
pub struct Phred(u8);

impl Phred {
    /// Creates a score. Returns `None` above [`MAX_Q`].
    #[must_use]
    pub const fn new(q: u8) -> Option<Self> {
        if q <= MAX_Q { Some(Self(q)) } else { None }
    }

    /// Decodes a Phred+33 (Sanger) byte. Returns `None` outside `!..=~`.
    #[must_use]
    pub const fn from_phred33(b: u8) -> Option<Self> {
        if b < 33 { None } else { Self::new(b - 33) }
    }

    /// Encodes as a Phred+33 byte.
    #[must_use]
    pub const fn to_phred33(self) -> u8 {
        self.0 + 33
    }

    /// The integer score Q.
    #[must_use]
    pub const fn value(self) -> u8 {
        self.0
    }

    /// The error probability `10^(-Q/10)`.
    #[must_use]
    pub fn error_probability(self) -> f64 {
        libm::pow(10.0, -f64::from(self.0) / 10.0)
    }

    /// The score for an error probability, rounded to the nearest integer and
    /// clamped to [`MAX_Q`]. Returns `None` if `p` is not in `(0, 1]`.
    #[must_use]
    pub fn from_error_probability(p: f64) -> Option<Self> {
        if !(p > 0.0 && p <= 1.0) {
            return None;
        }
        let q = libm::round(-10.0 * libm::log10(p));
        Some(Self(q.clamp(0.0, f64::from(MAX_Q)) as u8))
    }
}

/// Decodes a Phred+33 quality string into scores.
///
/// # Errors
/// Returns [`InvalidSymbol`] for the first byte outside `!..=~`.
pub fn decode_phred33(text: &[u8]) -> Result<alloc::vec::Vec<Phred>, InvalidSymbol> {
    text.iter()
        .enumerate()
        .map(|(position, &byte)| Phred::from_phred33(byte).ok_or(InvalidSymbol { position, byte }))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred33_roundtrip() {
        let q = Phred::from_phred33(b'I').unwrap();
        assert_eq!(q.value(), 40);
        assert_eq!(q.to_phred33(), b'I');
        assert_eq!(Phred::from_phred33(b' '), None);
        assert_eq!(Phred::new(94), None);
    }

    #[test]
    fn error_probability_known_values() {
        assert!((Phred::new(10).unwrap().error_probability() - 0.1).abs() < 1e-12);
        assert!((Phred::new(20).unwrap().error_probability() - 0.01).abs() < 1e-12);
        assert!((Phred::new(30).unwrap().error_probability() - 0.001).abs() < 1e-12);
    }

    #[test]
    fn from_error_probability_inverts() {
        assert_eq!(Phred::from_error_probability(0.001).unwrap().value(), 30);
        assert_eq!(Phred::from_error_probability(0.0), None);
        assert_eq!(Phred::from_error_probability(1.0).unwrap().value(), 0);
        assert_eq!(Phred::from_error_probability(1e-30).unwrap().value(), MAX_Q);
    }

    #[test]
    fn decode_string_reports_position() {
        let err = decode_phred33(b"II I").unwrap_err();
        assert_eq!(err, InvalidSymbol { position: 2, byte: b' ' });
        assert_eq!(decode_phred33(b"!I").unwrap()[1].value(), 40);
    }
}
