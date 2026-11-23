//! CRAI index: one line per slice, gzip-compressed text.
//!
//! Each line has six tab-separated fields: reference id, 1-based alignment start,
//! alignment span, container byte offset, slice byte offset (relative to the
//! container's first block) and slice length in bytes. Unmapped slices use
//! reference id -1 and multi-reference slices use -2.
//!
//! A region query `[beg, end)` (0-based, half-open) on reference `r` returns the
//! slices on `r` whose 1-based span overlaps it, in file order.

use std::fmt::Write as _;
use std::io::{Read, Write};

use flate2::{Compression, read::GzDecoder, write::GzEncoder};

use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAI";

/// One index line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CraiEntry {
    /// Reference id, or -1 (unmapped) or -2 (multi-reference).
    pub ref_id: i32,
    /// 1-based alignment start, or 0 when none.
    pub start: i32,
    /// Alignment span on the reference.
    pub span: i32,
    /// Byte offset of the container header in the CRAM file.
    pub container_offset: u64,
    /// Byte offset of the slice, relative to the container's first block.
    pub slice_offset: u64,
    /// Slice length in bytes.
    pub slice_len: u64,
}

/// A CRAI index.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Crai {
    entries: Vec<CraiEntry>,
}

impl Crai {
    /// Builds an index from entries, in the order they appear in the CRAM file.
    #[must_use]
    pub fn from_entries(entries: Vec<CraiEntry>) -> Self {
        Self { entries }
    }

    /// All index entries, in file order.
    #[must_use]
    pub fn entries(&self) -> &[CraiEntry] {
        &self.entries
    }

    /// Slices on `ref_id` that overlap the 0-based half-open region `[beg, end)`.
    #[must_use]
    pub fn query(&self, ref_id: i32, beg: i64, end: i64) -> Vec<CraiEntry> {
        self.entries
            .iter()
            .filter(|e| e.ref_id == ref_id && overlaps(e, beg, end))
            .copied()
            .collect()
    }

    /// Serializes to gzip-compressed text.
    ///
    /// # Errors
    /// Returns an error only if compression fails.
    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut text = String::new();
        for e in &self.entries {
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}",
                e.ref_id, e.start, e.span, e.container_offset, e.slice_offset, e.slice_len
            );
        }
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(text.as_bytes())?;
        Ok(enc.finish()?)
    }

    /// Parses gzip-compressed index bytes.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for bad gzip data or a line that does not
    /// have six integer fields.
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut text = String::new();
        GzDecoder::new(data)
            .read_to_string(&mut text)
            .map_err(|e| SeqError::malformed(FORMAT, format!("gzip: {e}")))?;
        let mut entries = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 6 {
                return Err(SeqError::malformed(FORMAT, format!("expected 6 fields: {line}")));
            }
            let bad = |_| SeqError::malformed(FORMAT, format!("bad integer in line: {line}"));
            entries.push(CraiEntry {
                ref_id: f[0].parse().map_err(bad)?,
                start: f[1].parse().map_err(bad)?,
                span: f[2].parse().map_err(bad)?,
                container_offset: f[3].parse().map_err(bad)?,
                slice_offset: f[4].parse().map_err(bad)?,
                slice_len: f[5].parse().map_err(bad)?,
            });
        }
        Ok(Self { entries })
    }
}

/// Whether a slice's 1-based span `[start, start+span-1]` overlaps `[beg, end)` (0-based).
fn overlaps(e: &CraiEntry, beg: i64, end: i64) -> bool {
    if e.start <= 0 {
        return false;
    }
    let s = i64::from(e.start) - 1;
    let t = s + i64::from(e.span.max(1));
    s < end && beg < t
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(ref_id: i32, start: i32, span: i32, n: u64) -> CraiEntry {
        CraiEntry {
            ref_id,
            start,
            span,
            container_offset: n * 100,
            slice_offset: 5,
            slice_len: 90 + n,
        }
    }

    #[test]
    fn roundtrip_text_format() {
        let idx = Crai::from_entries(vec![entry(0, 1, 100, 1), entry(-1, 0, 0, 2)]);
        let bytes = idx.to_bytes().unwrap();
        assert_eq!(Crai::from_bytes(&bytes).unwrap(), idx);
    }

    #[test]
    fn query_matches_brute_force_scan() {
        // Deterministic pseudo-random slices on two references.
        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = move |m: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % m
        };
        let mut entries = Vec::new();
        for n in 0..200u64 {
            let ref_id = (next(2)) as i32;
            let start = 1 + next(10_000) as i32;
            let span = 1 + next(500) as i32;
            entries.push(entry(ref_id, start, span, n));
        }
        let idx = Crai::from_entries(entries.clone());
        for _ in 0..200 {
            let ref_id = next(2) as i32;
            let beg = next(10_000) as i64;
            let end = beg + 1 + next(800) as i64;
            let expected: Vec<CraiEntry> = entries
                .iter()
                .filter(|e| e.ref_id == ref_id)
                .filter(|e| {
                    // Brute force on 1-based inclusive coordinates.
                    let first = i64::from(e.start);
                    let last = first + i64::from(e.span) - 1;
                    first <= end && last > beg
                })
                .copied()
                .collect();
            assert_eq!(idx.query(ref_id, beg, end), expected, "ref {ref_id} [{beg},{end})");
        }
    }

    #[test]
    fn malformed_line_is_rejected() {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(b"0\t1\t2\n").unwrap();
        let bytes = enc.finish().unwrap();
        assert!(matches!(Crai::from_bytes(&bytes), Err(SeqError::Malformed { .. })));
    }
}
