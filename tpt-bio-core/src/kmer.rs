//! K-mer encoding and iteration.
//!
//! A k-mer of up to 32 bases packs into a `u64` at 2 bits per base, with the
//! first base in the most significant position. Windows that contain a
//! non-ACGT byte are skipped.

use crate::error::InvalidSymbol;
use crate::nucleotide::Base;

/// Largest k that fits in a `u64` at 2 bits per base.
pub const MAX_K: usize = 32;

/// Errors from k-mer APIs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerError {
    /// k is zero or greater than [`MAX_K`].
    InvalidK(usize),
    /// The sequence contains a byte that is not an unambiguous base.
    Invalid(InvalidSymbol),
}

impl core::fmt::Display for KmerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidK(k) => write!(f, "k must be between 1 and {MAX_K}, got {k}"),
            Self::Invalid(e) => e.fmt(f),
        }
    }
}

impl core::error::Error for KmerError {}

/// Encodes a short sequence of unambiguous bases as a packed `u64`.
///
/// # Errors
/// Returns [`KmerError::InvalidK`] when the length is 0 or above 32, and
/// [`KmerError::Invalid`] for the first byte that is not ACGT.
pub fn encode(seq: &[u8]) -> Result<u64, KmerError> {
    if seq.is_empty() || seq.len() > MAX_K {
        return Err(KmerError::InvalidK(seq.len()));
    }
    let mut v = 0u64;
    for (i, &b) in seq.iter().enumerate() {
        let base = Base::from_ascii(b)
            .ok_or(KmerError::Invalid(InvalidSymbol { position: i, byte: b }))?;
        v = (v << 2) | u64::from(base.to_2bit());
    }
    Ok(v)
}

/// Reverse complement of a packed k-mer.
#[must_use]
pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    debug_assert!((1..=MAX_K).contains(&k));
    let mut x = kmer;
    let mut out = 0u64;
    for _ in 0..k {
        out = (out << 2) | (3 - (x & 3));
        x >>= 2;
    }
    out
}

/// The smaller of a k-mer and its reverse complement.
#[must_use]
pub fn canonical(kmer: u64, k: usize) -> u64 {
    kmer.min(reverse_complement(kmer, k))
}

/// Iterator over the k-mers of an ASCII sequence.
///
/// Yields `(offset, kmer)` where `offset` is the 0-based start of the window.
/// Windows containing a non-ACGT byte are skipped. Call [`Kmers::canonical`] to
/// yield the strand-independent form instead.
#[derive(Debug, Clone)]
pub struct Kmers<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
    valid: usize,
    fwd: u64,
    rev: u64,
    mask: u64,
    canonical: bool,
}

impl<'a> Kmers<'a> {
    /// Creates an iterator over the k-mers of `seq`.
    ///
    /// # Errors
    /// Returns [`KmerError::InvalidK`] when `k` is not in `1..=32`.
    pub fn new(seq: &'a [u8], k: usize) -> Result<Self, KmerError> {
        if !(1..=MAX_K).contains(&k) {
            return Err(KmerError::InvalidK(k));
        }
        let mask = if k == MAX_K { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        Ok(Self { seq, k, pos: 0, valid: 0, fwd: 0, rev: 0, mask, canonical: false })
    }

    /// Yields the canonical form of each k-mer.
    #[must_use]
    pub fn canonical(mut self) -> Self {
        self.canonical = true;
        self
    }
}

impl Iterator for Kmers<'_> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let i = self.pos;
            self.pos += 1;
            let Some(base) = Base::from_ascii(self.seq[i]) else {
                self.valid = 0;
                continue;
            };
            let code = u64::from(base.to_2bit());
            self.fwd = ((self.fwd << 2) | code) & self.mask;
            self.rev = (self.rev >> 2) | ((3 - code) << (2 * (self.k - 1)));
            self.valid += 1;
            if self.valid >= self.k {
                let kmer = if self.canonical { self.fwd.min(self.rev) } else { self.fwd };
                return Some((i + 1 - self.k, kmer));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;

    #[test]
    fn encode_and_revcomp() {
        let v = encode(b"ACGT").unwrap();
        assert_eq!(v, 0b00_01_10_11);
        assert_eq!(reverse_complement(v, 4), encode(b"ACGT").unwrap());
        assert_eq!(reverse_complement(encode(b"AAAC").unwrap(), 4), encode(b"GTTT").unwrap());
        assert_eq!(
            encode(b"ACNT").unwrap_err(),
            KmerError::Invalid(InvalidSymbol { position: 2, byte: b'N' })
        );
        assert_eq!(encode(b""), Err(KmerError::InvalidK(0)));
    }

    #[test]
    fn iterates_windows_and_skips_ambiguous() {
        let got: Vec<(usize, u64)> = Kmers::new(b"ACGNACG", 3).unwrap().collect();
        let expect = vec![(0, encode(b"ACG").unwrap()), (4, encode(b"ACG").unwrap())];
        assert_eq!(got, expect);
    }

    #[test]
    fn canonical_matches_min_of_strands() {
        let seq = b"TTTAAAGGG";
        for (off, kmer) in Kmers::new(seq, 5).unwrap().canonical() {
            let window = &seq[off..off + 5];
            let rc = crate::nucleotide::reverse_complement(window).unwrap();
            let rc_code = encode(&rc).unwrap();
            assert_eq!(kmer, kmer.min(rc_code));
        }
    }

    #[test]
    fn k_bounds() {
        assert!(Kmers::new(b"ACGT", 0).is_err());
        assert!(Kmers::new(b"ACGT", 33).is_err());
        assert!(Kmers::new(&[b'A'; 40], 32).unwrap().count() == 9);
    }
}
