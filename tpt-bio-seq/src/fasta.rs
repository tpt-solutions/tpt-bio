//! FASTA reading and writing, and the samtools-compatible `.fai` index.
//!
//! Sequence lines may be wrapped at any width. The index records, per sequence,
//! the offset of its first base and the line geometry, which is enough to fetch
//! any subsequence without scanning the file.

use std::io::{BufRead, Write};

use crate::error::{Result, SeqError};

const FORMAT: &str = "FASTA";

/// One FASTA entry with its sequence joined into a single buffer.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FastaRecord {
    /// Name: the header up to the first whitespace.
    pub name: String,
    /// Free-text description after the name, if any.
    pub description: Option<String>,
    /// Sequence bytes with line breaks removed.
    pub seq: Vec<u8>,
}

/// Streaming FASTA reader.
#[derive(Debug)]
pub struct FastaReader<R> {
    inner: R,
    pending_header: Option<String>,
    line: Vec<u8>,
    started: bool,
}

impl<R: BufRead> FastaReader<R> {
    /// Wraps a buffered reader.
    pub fn new(inner: R) -> Self {
        Self { inner, pending_header: None, line: Vec::new(), started: false }
    }

    fn read_line(&mut self) -> Result<bool> {
        self.line.clear();
        if self.inner.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(true)
    }

    /// The next record, or `None` at end of input.
    ///
    /// # Errors
    /// Returns an error if the input does not start with `>` or on I/O failure.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>> {
        let header = match self.pending_header.take() {
            Some(h) => h,
            None => loop {
                if !self.read_line()? {
                    return Ok(None);
                }
                if self.line.is_empty() {
                    continue;
                }
                match self.line.strip_prefix(b">") {
                    Some(rest) => break String::from_utf8_lossy(rest).into_owned(),
                    None if !self.started => {
                        return Err(SeqError::malformed(FORMAT, "input does not start with '>'"));
                    }
                    None => {
                        return Err(SeqError::malformed(FORMAT, "sequence line outside a record"));
                    }
                }
            },
        };
        self.started = true;
        let (name, description) = split_header(&header);
        let mut seq = Vec::new();
        loop {
            if !self.read_line()? {
                break;
            }
            if self.line.first() == Some(&b'>') {
                self.pending_header = Some(String::from_utf8_lossy(&self.line[1..]).into_owned());
                break;
            }
            seq.extend_from_slice(&self.line);
        }
        Ok(Some(FastaRecord { name, description, seq }))
    }
}

fn split_header(header: &str) -> (String, Option<String>) {
    let header = header.trim_end();
    match header.split_once(|c: char| c.is_ascii_whitespace()) {
        Some((name, desc)) => {
            let desc = desc.trim_start();
            (name.to_string(), (!desc.is_empty()).then(|| desc.to_string()))
        }
        None => (header.to_string(), None),
    }
}

/// FASTA writer with a fixed line width.
#[derive(Debug)]
pub struct FastaWriter<W: Write> {
    inner: W,
    width: usize,
}

impl<W: Write> FastaWriter<W> {
    /// Creates a writer that wraps sequence lines at `width` bases.
    ///
    /// # Panics
    /// Panics if `width` is zero.
    pub fn new(inner: W, width: usize) -> Self {
        assert!(width > 0, "line width must be non-zero");
        Self { inner, width }
    }

    /// Writes one record.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn write_record(&mut self, rec: &FastaRecord) -> Result<()> {
        self.inner.write_all(b">")?;
        self.inner.write_all(rec.name.as_bytes())?;
        if let Some(desc) = &rec.description {
            self.inner.write_all(b" ")?;
            self.inner.write_all(desc.as_bytes())?;
        }
        self.inner.write_all(b"\n")?;
        for line in rec.seq.chunks(self.width) {
            self.inner.write_all(line)?;
            self.inner.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Flushes and returns the inner writer.
    ///
    /// # Errors
    /// Returns an error if the flush fails.
    pub fn finish(mut self) -> Result<W> {
        self.inner.flush()?;
        Ok(self.inner)
    }
}

/// One line of a `.fai` index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    /// Sequence name.
    pub name: String,
    /// Sequence length in bases.
    pub length: u64,
    /// Byte offset of the first base.
    pub offset: u64,
    /// Bases per line.
    pub line_bases: u64,
    /// Bytes per line, including the line terminator.
    pub line_width: u64,
}

/// A `.fai` index, in samtools format.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Fai {
    entries: Vec<FaiEntry>,
}

impl Fai {
    /// Builds an index over an in-memory FASTA file.
    ///
    /// # Errors
    /// Returns an error if a sequence has irregular line lengths, which the
    /// `.fai` format cannot express, or if the input is not FASTA.
    pub fn build(data: &[u8]) -> Result<Self> {
        let mut entries: Vec<FaiEntry> = Vec::new();
        let mut pos = 0usize;
        // Current sequence state.
        let mut current: Option<FaiEntry> = None;
        let mut line_lens: Vec<usize> = Vec::new();
        let mut terminator_len = 1u64;

        let mut finish = |cur: &mut Option<FaiEntry>, lens: &mut Vec<usize>| -> Result<()> {
            if let Some(mut e) = cur.take() {
                if let Some((last, body)) = lens.split_last()
                    && (body.iter().any(|&l| l as u64 != e.line_bases)
                        || *last as u64 > e.line_bases)
                {
                    return Err(SeqError::malformed(
                        FORMAT,
                        format!("sequence {} has irregular line lengths", e.name),
                    ));
                }
                e.length = lens.iter().map(|&l| l as u64).sum();
                entries.push(e);
            }
            lens.clear();
            Ok(())
        };

        while pos < data.len() {
            let end = data[pos..].iter().position(|&b| b == b'\n').map_or(data.len(), |i| pos + i);
            let mut line = &data[pos..end];
            let had_newline = end < data.len();
            let line_start = pos;
            pos = if had_newline { end + 1 } else { end };
            if line.last() == Some(&b'\r') {
                line = &line[..line.len() - 1];
                terminator_len = 2;
            }
            if let Some(header) = line.strip_prefix(b">") {
                finish(&mut current, &mut line_lens)?;
                let name = split_header(&String::from_utf8_lossy(header)).0;
                current = Some(FaiEntry {
                    name,
                    length: 0,
                    offset: pos as u64,
                    line_bases: 0,
                    line_width: 0,
                });
                continue;
            }
            if line.is_empty() {
                continue;
            }
            let Some(cur) = current.as_mut() else {
                return Err(SeqError::malformed_at(
                    FORMAT,
                    line_start as u64,
                    "sequence before header",
                ));
            };
            if cur.line_bases == 0 {
                cur.line_bases = line.len() as u64;
                cur.line_width = (line.len() as u64) + terminator_len;
            }
            line_lens.push(line.len());
        }
        finish(&mut current, &mut line_lens)?;
        Ok(Self { entries })
    }

    /// The entries in file order.
    #[must_use]
    pub fn entries(&self) -> &[FaiEntry] {
        &self.entries
    }

    /// Looks up an entry by name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        self.entries.iter().find(|e| e.name == name)
    }

    /// Serializes in the samtools `.fai` text format.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        for e in &self.entries {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                e.name, e.length, e.offset, e.line_bases, e.line_width
            )?;
        }
        Ok(())
    }

    /// Parses the samtools `.fai` text format.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for lines without five numeric fields.
    pub fn read_from(text: &str) -> Result<Self> {
        let mut entries = Vec::new();
        for (i, line) in text.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let bad =
                || SeqError::malformed(FORMAT, format!("fai line {}: expected 5 fields", i + 1));
            let [name, length, offset, line_bases, line_width] = f.as_slice() else {
                return Err(bad());
            };
            let num = |s: &str| s.parse::<u64>().map_err(|_| bad());
            entries.push(FaiEntry {
                name: (*name).to_string(),
                length: num(length)?,
                offset: num(offset)?,
                line_bases: num(line_bases)?,
                line_width: num(line_width)?,
            });
        }
        Ok(Self { entries })
    }

    /// Fetches bases `[beg, end)` (0-based, half-open) of a sequence from the
    /// FASTA bytes the index was built from. The end is clamped to the length.
    ///
    /// # Errors
    /// Returns [`SeqError::NotFound`] for an unknown name, or [`SeqError::Malformed`]
    /// if the index does not match `data`.
    pub fn fetch(&self, data: &[u8], name: &str, beg: u64, end: u64) -> Result<Vec<u8>> {
        let e = self.get(name).ok_or_else(|| SeqError::NotFound(format!("sequence {name}")))?;
        let end = end.min(e.length);
        if beg >= end {
            return Ok(Vec::new());
        }
        let mut out = Vec::with_capacity((end - beg) as usize);
        let mut base = beg;
        while base < end {
            let line = base / e.line_bases;
            let col = base % e.line_bases;
            let take = (e.line_bases - col).min(end - base);
            let start = e.offset + line * e.line_width + col;
            let slice = data.get(start as usize..(start + take) as usize).ok_or_else(|| {
                SeqError::malformed(FORMAT, "index does not match the FASTA data")
            })?;
            out.extend_from_slice(slice);
            base += take;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FA: &[u8] = b">chr1 test contig\nACGTACGT\nACGTAC\n>chr2\nGGGG\nCC\n";

    #[test]
    fn reads_multiline_records() {
        let mut r = FastaReader::new(FA);
        let a = r.next_record().unwrap().unwrap();
        assert_eq!(a.name, "chr1");
        assert_eq!(a.description.as_deref(), Some("test contig"));
        assert_eq!(a.seq, b"ACGTACGTACGTAC");
        let b = r.next_record().unwrap().unwrap();
        assert_eq!(b.seq, b"GGGGCC");
        assert!(r.next_record().unwrap().is_none());
    }

    #[test]
    fn rejects_input_without_header() {
        assert!(FastaReader::new(&b"ACGT\n"[..]).next_record().is_err());
    }

    #[test]
    fn write_then_read_roundtrip() {
        let mut r = FastaReader::new(FA);
        let mut w = FastaWriter::new(Vec::new(), 5);
        let mut recs = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            w.write_record(&rec).unwrap();
            recs.push(rec);
        }
        let out = w.finish().unwrap();
        let mut back = FastaReader::new(&out[..]);
        for expect in recs {
            assert_eq!(back.next_record().unwrap().unwrap(), expect);
        }
    }

    #[test]
    fn fai_build_fetch_and_text_roundtrip() {
        let fai = Fai::build(FA).unwrap();
        assert_eq!(fai.entries().len(), 2);
        assert_eq!(fai.get("chr1").unwrap().length, 14);
        assert_eq!(fai.get("chr1").unwrap().line_bases, 8);
        assert_eq!(fai.get("chr2").unwrap().length, 6);

        // Fetch spans a line break and matches a plain slice of the sequence.
        assert_eq!(fai.fetch(FA, "chr1", 6, 11).unwrap(), b"GTACG");
        assert_eq!(fai.fetch(FA, "chr1", 0, 14).unwrap(), b"ACGTACGTACGTAC");
        assert_eq!(fai.fetch(FA, "chr2", 3, 100).unwrap(), b"GCC");
        assert!(fai.fetch(FA, "missing", 0, 1).is_err());

        let mut text = Vec::new();
        fai.write_to(&mut text).unwrap();
        let back = Fai::read_from(std::str::from_utf8(&text).unwrap()).unwrap();
        assert_eq!(back, fai);
    }

    #[test]
    fn fai_rejects_irregular_wrapping() {
        assert!(Fai::build(b">x\nACGT\nAC\nACGT\n").is_err());
    }
}
