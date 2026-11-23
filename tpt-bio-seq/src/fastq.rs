//! FASTQ parsing and writing.
//!
//! Four-line records only (no wrapped sequence lines). CRLF line endings are
//! accepted. [`FastqParser`] is zero-copy over a byte slice, which suits
//! memory-mapped input. [`FastqReader`] streams owned records from any
//! [`BufRead`]. Gzip and BGZF input are decoded transparently by [`open_reader`].

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use crate::error::{Result, SeqError};

const FORMAT: &str = "FASTQ";

/// A borrowed FASTQ record. `name` excludes the leading `@`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastqRecord<'a> {
    /// Header line without the leading `@`.
    pub name: &'a [u8],
    /// Bases.
    pub seq: &'a [u8],
    /// Phred+33 qualities, same length as `seq`.
    pub qual: &'a [u8],
}

impl FastqRecord<'_> {
    /// The read identifier: the header up to the first space or tab.
    #[must_use]
    pub fn id(&self) -> &[u8] {
        id_of(self.name)
    }

    /// Copies the record into an owned buffer.
    #[must_use]
    pub fn to_owned_record(&self) -> FastqBuf {
        FastqBuf { name: self.name.to_vec(), seq: self.seq.to_vec(), qual: self.qual.to_vec() }
    }
}

/// An owned FASTQ record.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastqBuf {
    /// Header line without the leading `@`.
    pub name: Vec<u8>,
    /// Bases.
    pub seq: Vec<u8>,
    /// Phred+33 qualities.
    pub qual: Vec<u8>,
}

impl FastqBuf {
    /// Borrows the record.
    #[must_use]
    pub fn as_record(&self) -> FastqRecord<'_> {
        FastqRecord { name: &self.name, seq: &self.seq, qual: &self.qual }
    }
}

/// Zero-copy FASTQ parser over a byte slice.
#[derive(Debug, Clone)]
pub struct FastqParser<'a> {
    data: &'a [u8],
    pos: usize,
    record: u64,
}

impl<'a> FastqParser<'a> {
    /// Creates a parser over the full slice.
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0, record: 0 }
    }

    fn line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.data.len() {
            return None;
        }
        let rest = &self.data[self.pos..];
        let (mut line, advance) = match rest.iter().position(|&b| b == b'\n') {
            Some(i) => (&rest[..i], i + 1),
            None => (rest, rest.len()),
        };
        self.pos += advance;
        if let Some((&b'\r', head)) = line.split_last() {
            line = head;
        }
        Some(line)
    }

    fn fail(&self, reason: impl Into<String>) -> SeqError {
        SeqError::malformed_at(
            FORMAT,
            self.pos as u64,
            format!("record {}: {}", self.record + 1, reason.into()),
        )
    }
}

impl<'a> Iterator for FastqParser<'a> {
    type Item = Result<FastqRecord<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        // Skip blank lines between records.
        let header = loop {
            let line = self.line()?;
            if !line.is_empty() {
                break line;
            }
        };
        let Some(name) = header.strip_prefix(b"@") else {
            return Some(Err(self.fail("header line does not start with '@'")));
        };
        let (Some(seq), Some(plus), Some(qual)) = (self.line(), self.line(), self.line()) else {
            return Some(Err(self.fail("unexpected end of file inside record")));
        };
        if !plus.starts_with(b"+") {
            return Some(Err(self.fail("separator line does not start with '+'")));
        }
        if qual.len() != seq.len() {
            return Some(Err(self.fail(format!(
                "sequence length {} differs from quality length {}",
                seq.len(),
                qual.len()
            ))));
        }
        self.record += 1;
        Some(Ok(FastqRecord { name, seq, qual }))
    }
}

/// Streaming FASTQ reader over any buffered source.
#[derive(Debug)]
pub struct FastqReader<R> {
    inner: R,
    line: Vec<u8>,
    record: u64,
}

impl<R: BufRead> FastqReader<R> {
    /// Wraps a buffered reader.
    pub fn new(inner: R) -> Self {
        Self { inner, line: Vec::new(), record: 0 }
    }

    /// Reads one line into `self.line` without its terminator. `false` at EOF.
    fn read_line(&mut self) -> Result<bool> {
        self.line.clear();
        if self.inner.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }
        if self.line.last() == Some(&b'\n') {
            self.line.pop();
        }
        if self.line.last() == Some(&b'\r') {
            self.line.pop();
        }
        Ok(true)
    }

    /// Reads the next record into `rec`, reusing its allocations. `false` at EOF.
    pub fn read_record(&mut self, rec: &mut FastqBuf) -> Result<bool> {
        loop {
            if !self.read_line()? {
                return Ok(false);
            }
            if !self.line.is_empty() {
                break;
            }
        }
        if self.line.first() != Some(&b'@') {
            return Err(self.fail("header line does not start with '@'"));
        }
        rec.name.clear();
        rec.name.extend_from_slice(&self.line[1..]);

        rec.seq.clear();
        if !self.read_line()? {
            return Err(self.fail("unexpected end of file inside record"));
        }
        rec.seq.extend_from_slice(&self.line);

        if !self.read_line()? || self.line.first() != Some(&b'+') {
            return Err(self.fail("separator line does not start with '+'"));
        }

        rec.qual.clear();
        if !self.read_line()? {
            return Err(self.fail("unexpected end of file inside record"));
        }
        rec.qual.extend_from_slice(&self.line);
        if rec.qual.len() != rec.seq.len() {
            return Err(self.fail(format!(
                "sequence length {} differs from quality length {}",
                rec.seq.len(),
                rec.qual.len()
            )));
        }
        self.record += 1;
        Ok(true)
    }

    /// Reads up to `n` records. The result is shorter only at end of input.
    pub fn next_batch(&mut self, n: usize) -> Result<Vec<FastqBuf>> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let mut rec = FastqBuf::default();
            if !self.read_record(&mut rec)? {
                break;
            }
            out.push(rec);
        }
        Ok(out)
    }

    fn fail(&self, reason: impl Into<String>) -> SeqError {
        SeqError::malformed(FORMAT, format!("record {}: {}", self.record + 1, reason.into()))
    }
}

/// Reads records in pairs from two files (R1 and R2) and checks that the
/// identifiers match.
#[derive(Debug)]
pub struct PairedFastqReader<R1, R2> {
    first: FastqReader<R1>,
    second: FastqReader<R2>,
}

impl<R1: BufRead, R2: BufRead> PairedFastqReader<R1, R2> {
    /// Pairs two readers.
    pub fn new(first: R1, second: R2) -> Self {
        Self { first: FastqReader::new(first), second: FastqReader::new(second) }
    }

    /// The next pair, or `None` when both inputs are exhausted.
    ///
    /// # Errors
    /// Returns an error if one input ends before the other, or if the read IDs differ.
    pub fn next_pair(&mut self) -> Result<Option<(FastqBuf, FastqBuf)>> {
        let mut a = FastqBuf::default();
        let mut b = FastqBuf::default();
        let has_a = self.first.read_record(&mut a)?;
        let has_b = self.second.read_record(&mut b)?;
        match (has_a, has_b) {
            (false, false) => Ok(None),
            (true, true) => {
                if pair_id(&a.name) != pair_id(&b.name) {
                    return Err(SeqError::malformed(
                        FORMAT,
                        format!(
                            "paired reads out of sync: {} vs {}",
                            String::from_utf8_lossy(&a.name),
                            String::from_utf8_lossy(&b.name)
                        ),
                    ));
                }
                Ok(Some((a, b)))
            }
            _ => Err(SeqError::malformed(FORMAT, "paired inputs have different record counts")),
        }
    }
}

/// Reads interleaved pairs (R1, R2, R1, R2, ...) from one source.
#[derive(Debug)]
pub struct InterleavedFastqReader<R> {
    inner: FastqReader<R>,
}

impl<R: BufRead> InterleavedFastqReader<R> {
    /// Wraps an interleaved stream.
    pub fn new(inner: R) -> Self {
        Self { inner: FastqReader::new(inner) }
    }

    /// The next pair, or `None` at end of input.
    ///
    /// # Errors
    /// Returns an error if the stream has an odd number of records or the IDs differ.
    pub fn next_pair(&mut self) -> Result<Option<(FastqBuf, FastqBuf)>> {
        let mut a = FastqBuf::default();
        if !self.inner.read_record(&mut a)? {
            return Ok(None);
        }
        let mut b = FastqBuf::default();
        if !self.inner.read_record(&mut b)? {
            return Err(SeqError::malformed(FORMAT, "interleaved input has an odd record count"));
        }
        if pair_id(&a.name) != pair_id(&b.name) {
            return Err(SeqError::malformed(FORMAT, "interleaved pair IDs do not match"));
        }
        Ok(Some((a, b)))
    }
}

/// The read name with a trailing `/1` or `/2` removed, so mates compare equal.
fn pair_id(name: &[u8]) -> &[u8] {
    match id_of(name) {
        [head @ .., b'/', b'1' | b'2'] => head,
        id => id,
    }
}

/// The header up to the first space or tab.
fn id_of(name: &[u8]) -> &[u8] {
    let end = name.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(name.len());
    &name[..end]
}

/// FASTQ writer.
#[derive(Debug)]
pub struct FastqWriter<W: Write> {
    inner: W,
}

impl<W: Write> FastqWriter<W> {
    /// Wraps a writer.
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes one record.
    ///
    /// # Errors
    /// Returns an error if `seq` and `qual` differ in length, or on I/O failure.
    pub fn write_record(&mut self, rec: &FastqRecord<'_>) -> Result<()> {
        if rec.seq.len() != rec.qual.len() {
            return Err(SeqError::malformed(FORMAT, "sequence and quality lengths differ"));
        }
        self.inner.write_all(b"@")?;
        self.inner.write_all(rec.name)?;
        self.inner.write_all(b"\n")?;
        self.inner.write_all(rec.seq)?;
        self.inner.write_all(b"\n+\n")?;
        self.inner.write_all(rec.qual)?;
        self.inner.write_all(b"\n")?;
        Ok(())
    }

    /// Flushes and returns the inner writer.
    pub fn finish(mut self) -> Result<W> {
        self.inner.flush()?;
        Ok(self.inner)
    }
}

/// Opens a file for buffered reading, transparently decompressing gzip and BGZF.
///
/// Compression is detected from the magic bytes, not the file extension.
///
/// # Errors
/// Returns an I/O error if the file cannot be opened or read.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut file = BufReader::with_capacity(1 << 16, File::open(path)?);
    let magic = file.fill_buf()?;
    let gzipped = magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
    if gzipped {
        Ok(Box::new(BufReader::with_capacity(1 << 16, flate2::bufread::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(file))
    }
}

/// Reads a whole reader into memory as a byte vector, decompressing gzip/BGZF.
pub fn read_all_maybe_compressed(mut reader: impl Read) -> Result<Vec<u8>> {
    let mut raw = Vec::new();
    reader.read_to_end(&mut raw)?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut out)?;
        return Ok(out);
    }
    Ok(raw)
}

/// Streams FASTQ in batches with a reader thread and a bounded channel.
///
/// A producer thread parses batches of `batch_size` records and sends them over
/// a channel that holds at most `capacity` batches, so a slow consumer applies
/// backpressure and memory stays bounded. The consumer calls `process` on each
/// batch, in input order, and passes the result to `sink`. `process` may use
/// rayon internally to parallelize within a batch.
///
/// # Errors
/// Returns the first parse or I/O error. Batches before the error are still
/// delivered to `sink`.
///
/// # Panics
/// Panics if `batch_size` or `capacity` is zero.
pub fn process_batches<R, T, P, S>(
    reader: R,
    batch_size: usize,
    capacity: usize,
    process: P,
    mut sink: S,
) -> Result<()>
where
    R: BufRead + Send,
    T: Send,
    P: Fn(&[FastqBuf]) -> T,
    S: FnMut(T),
{
    assert!(batch_size > 0 && capacity > 0, "batch_size and capacity must be non-zero");
    let (tx, rx) = crossbeam::channel::bounded::<Vec<FastqBuf>>(capacity);
    std::thread::scope(|scope| {
        let producer = scope.spawn(move || -> Result<()> {
            let mut parser = FastqReader::new(reader);
            loop {
                let batch = parser.next_batch(batch_size)?;
                if batch.is_empty() {
                    return Ok(());
                }
                if tx.send(batch).is_err() {
                    // The consumer stopped early; nothing more to do.
                    return Ok(());
                }
            }
        });
        for batch in &rx {
            sink(process(&batch));
        }
        producer
            .join()
            .unwrap_or_else(|_| Err(SeqError::malformed(FORMAT, "reader thread panicked")))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &[u8] = b"@r1 extra\nACGT\n+\nIIII\n@r2\nNNA\n+\n!!#\n";

    #[test]
    fn zero_copy_parser() {
        let recs: Vec<_> = FastqParser::new(SAMPLE).collect::<Result<_>>().unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id(), b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[1].qual, b"!!#");
    }

    #[test]
    fn crlf_and_blank_lines_accepted() {
        let data = b"\r\n@a\r\nAC\r\n+\r\nII\r\n\r\n";
        let recs: Vec<_> = FastqParser::new(data).collect::<Result<_>>().unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].seq, b"AC");
    }

    #[test]
    fn streaming_matches_zero_copy() {
        let mut reader = FastqReader::new(SAMPLE);
        let mut buf = FastqBuf::default();
        let mut got = Vec::new();
        while reader.read_record(&mut buf).unwrap() {
            got.push(buf.clone());
        }
        let expect: Vec<_> =
            FastqParser::new(SAMPLE).map(|r| r.unwrap().to_owned_record()).collect();
        assert_eq!(got, expect);
    }

    #[test]
    fn malformed_inputs_are_rejected() {
        assert!(FastqParser::new(b"r1\nAC\n+\nII\n").next().unwrap().is_err());
        assert!(FastqParser::new(b"@r1\nAC\n+\nI\n").next().unwrap().is_err());
        assert!(FastqParser::new(b"@r1\nAC\nX\nII\n").next().unwrap().is_err());
        assert!(FastqParser::new(b"@r1\nAC\n").next().unwrap().is_err());
    }

    #[test]
    fn write_then_read_roundtrip() {
        let mut w = FastqWriter::new(Vec::new());
        for rec in FastqParser::new(SAMPLE) {
            w.write_record(&rec.unwrap()).unwrap();
        }
        let out = w.finish().unwrap();
        assert_eq!(out, SAMPLE);
    }

    #[test]
    fn paired_and_interleaved() {
        let r1 = b"@p/1\nAC\n+\nII\n@q/1\nGG\n+\nII\n";
        let r2 = b"@p/2\nTT\n+\nII\n@q/2\nCC\n+\nII\n";
        let mut pr = PairedFastqReader::new(&r1[..], &r2[..]);
        assert!(pr.next_pair().unwrap().is_some());
        assert!(pr.next_pair().unwrap().is_some());
        assert!(pr.next_pair().unwrap().is_none());

        let inter = b"@p/1\nAC\n+\nII\n@p/2\nTT\n+\nII\n";
        let mut ir = InterleavedFastqReader::new(&inter[..]);
        let (a, b) = ir.next_pair().unwrap().unwrap();
        assert_eq!(a.seq, b"AC");
        assert_eq!(b.seq, b"TT");

        let bad = b"@x/1\nAC\n+\nII\n";
        let mut bad_pair = PairedFastqReader::new(&bad[..], &r2[..]);
        assert!(bad_pair.next_pair().is_err());
    }

    #[test]
    fn streaming_batches_preserve_order() {
        let mut data = Vec::new();
        for i in 0..25u32 {
            data.extend_from_slice(
                format!(
                    "@r{i}
ACGT
+
IIII
"
                )
                .as_bytes(),
            );
        }
        let mut seen = Vec::new();
        process_batches(
            &data[..],
            4,
            2,
            |batch| batch.iter().map(|r| r.name.clone()).collect::<Vec<_>>(),
            |names| seen.extend(names),
        )
        .unwrap();
        assert_eq!(seen.len(), 25);
        assert_eq!(seen[0], b"r0");
        assert_eq!(seen[24], b"r24");
    }

    #[test]
    fn streaming_batches_report_parse_errors() {
        let data = b"@r0
ACGT
+
IIII
@bad
AC
+
I
";
        let result = process_batches(&data[..], 1, 1, <[FastqBuf]>::len, |_| {});
        assert!(result.is_err());
    }

    #[test]
    fn batch_reading() {
        let mut reader = FastqReader::new(SAMPLE);
        assert_eq!(reader.next_batch(1).unwrap().len(), 1);
        assert_eq!(reader.next_batch(10).unwrap().len(), 1);
        assert!(reader.next_batch(10).unwrap().is_empty());
    }
}
