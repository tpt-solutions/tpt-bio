//! Genomic coordinates and regions.
//!
//! Internally every [`Interval`] is **0-based and half-open**, `[start, end)`,
//! which is the convention of BAM, BED and the in-memory APIs. Text regions such
//! as `chr1:1000-2000` are **1-based and inclusive**, following samtools. Use
//! [`Interval::from_one_based_inclusive`] when converting between the two.

use alloc::string::{String, ToString};

use crate::error::RegionParseError;

/// A 0-based, half-open interval `[start, end)` on one contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval {
    start: u64,
    end: u64,
}

impl Interval {
    /// Creates `[start, end)`. Returns `None` if `end < start`.
    #[must_use]
    pub const fn new(start: u64, end: u64) -> Option<Self> {
        if end < start { None } else { Some(Self { start, end }) }
    }

    /// Creates the interval for 1-based inclusive coordinates `start..=end`.
    /// Returns `None` if `start` is 0 or `end < start`.
    #[must_use]
    pub const fn from_one_based_inclusive(start: u64, end: u64) -> Option<Self> {
        if start == 0 || end < start { None } else { Some(Self { start: start - 1, end }) }
    }

    /// 0-based inclusive start.
    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    /// 0-based exclusive end.
    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    /// Number of positions covered.
    #[must_use]
    pub const fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the interval covers no positions.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// Whether the 0-based position is inside the interval.
    #[must_use]
    pub const fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos < self.end
    }

    /// Whether two intervals share at least one position.
    #[must_use]
    pub const fn overlaps(&self, other: &Self) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// The shared part of two intervals, if any.
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Option<Self> {
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        Self::new(start, end).filter(|i| !i.is_empty())
    }
}

/// A contig name plus an optional interval. `None` means the whole contig.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    contig: String,
    interval: Option<Interval>,
}

impl Region {
    /// Creates a region covering a whole contig.
    #[must_use]
    pub fn whole_contig(contig: &str) -> Self {
        Self { contig: contig.to_string(), interval: None }
    }

    /// Creates a region from a contig and a 0-based half-open interval.
    #[must_use]
    pub fn new(contig: &str, interval: Interval) -> Self {
        Self { contig: contig.to_string(), interval: Some(interval) }
    }

    /// Parses `contig`, `contig:start`, or `contig:start-end`.
    ///
    /// Coordinates are 1-based and inclusive. Commas are allowed as thousands
    /// separators. If the text after the last `:` is not a coordinate range, the
    /// whole string is treated as a contig name. As in samtools, a numeric suffix
    /// is always read as a coordinate, so a name such as `HLA-A*01:01` must be
    /// requested with [`Region::whole_contig`].
    ///
    /// # Errors
    /// Returns [`RegionParseError`] when the contig is empty, a coordinate is not
    /// an integer, or the range is reversed or starts before 1.
    pub fn parse(text: &str) -> Result<Self, RegionParseError> {
        let Some(colon) = text.rfind(':') else {
            return Self::checked_contig(text).map(Self::whole_contig);
        };
        let (contig, range) = (&text[..colon], &text[colon + 1..]);
        if !looks_like_range(range) {
            return Self::checked_contig(text).map(Self::whole_contig);
        }
        Self::checked_contig(contig)?;
        let (start_text, end_text) = match range.split_once('-') {
            Some((s, e)) => (s, Some(e)),
            None => (range, None),
        };
        let start = parse_coord(start_text)?;
        let end = match end_text {
            Some(e) => parse_coord(e)?,
            None => u64::MAX,
        };
        let interval =
            Interval::from_one_based_inclusive(start, end).ok_or(RegionParseError::BadRange)?;
        Ok(Self::new(contig, interval))
    }

    fn checked_contig(name: &str) -> Result<&str, RegionParseError> {
        if name.is_empty() { Err(RegionParseError::EmptyContig) } else { Ok(name) }
    }

    /// The contig name.
    #[must_use]
    pub fn contig(&self) -> &str {
        &self.contig
    }

    /// The interval, or `None` for a whole contig.
    #[must_use]
    pub const fn interval(&self) -> Option<Interval> {
        self.interval
    }
}

fn looks_like_range(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit() || b == b'-' || b == b',')
}

fn parse_coord(s: &str) -> Result<u64, RegionParseError> {
    let digits: String = s.chars().filter(|&c| c != ',').collect();
    if digits.is_empty() {
        return Err(RegionParseError::BadNumber);
    }
    digits.parse::<u64>().map_err(|_| RegionParseError::BadNumber)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_conversions() {
        let i = Interval::from_one_based_inclusive(1, 10).unwrap();
        assert_eq!((i.start(), i.end(), i.len()), (0, 10, 10));
        assert!(i.contains(0) && i.contains(9) && !i.contains(10));
        assert_eq!(Interval::from_one_based_inclusive(0, 5), None);
        assert_eq!(Interval::new(5, 4), None);
    }

    #[test]
    fn interval_overlap_and_intersection() {
        let a = Interval::new(0, 10).unwrap();
        let b = Interval::new(5, 20).unwrap();
        let c = Interval::new(10, 20).unwrap();
        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c), "half-open intervals touching at 10 do not overlap");
        assert_eq!(a.intersect(&b), Interval::new(5, 10));
        assert_eq!(a.intersect(&c), None);
    }

    #[test]
    fn parse_regions() {
        let r = Region::parse("chr1:1,000-2000").unwrap();
        assert_eq!(r.contig(), "chr1");
        assert_eq!(r.interval(), Interval::new(999, 2000));

        let single = Region::parse("chrX:5").unwrap();
        assert_eq!(single.interval(), Interval::new(4, u64::MAX));

        let whole = Region::parse("chrM").unwrap();
        assert_eq!(whole.interval(), None);

        let named = Region::parse("chrUn:scaffold_1").unwrap();
        assert_eq!(named.contig(), "chrUn:scaffold_1");
        assert_eq!(named.interval(), None);
    }

    #[test]
    fn parse_region_errors() {
        assert_eq!(Region::parse(":100").unwrap_err(), RegionParseError::EmptyContig);
        assert_eq!(Region::parse("chr1:200-100").unwrap_err(), RegionParseError::BadRange);
        assert_eq!(Region::parse("chr1:0-100").unwrap_err(), RegionParseError::BadRange);
    }
}
