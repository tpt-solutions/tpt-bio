//! Error types shared across the tpt-bio crates.

use core::fmt;

/// A byte that is not a valid symbol for the alphabet being processed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidSymbol {
    /// Zero-based offset of the offending byte in the input.
    pub position: usize,
    /// The offending byte.
    pub byte: u8,
}

impl fmt::Display for InvalidSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid symbol 0x{:02X} ({:?}) at offset {}",
            self.byte, self.byte as char, self.position
        )
    }
}

impl core::error::Error for InvalidSymbol {}

/// Failure to parse a region string such as `chr1:1000-2000`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegionParseError {
    /// The contig name is empty.
    EmptyContig,
    /// A coordinate is not a non-negative integer.
    BadNumber,
    /// The range is reversed or starts before 1 (coordinates are 1-based).
    BadRange,
}

impl fmt::Display for RegionParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyContig => f.write_str("region has an empty contig name"),
            Self::BadNumber => f.write_str("region coordinate is not a valid integer"),
            Self::BadRange => f.write_str("region range is reversed or starts before 1"),
        }
    }
}

impl core::error::Error for RegionParseError {}
