//! VCF text reading and writing (VCF 4.x, uncompressed).
//!
//! The header is parsed into structured meta lines (`##INFO=<ID=...>` and so on).
//! Records are zero-copy views: each [`VcfRecord`] borrows its line from the
//! reader's buffer, and only the fixed numeric fields are converted.
//!
//! BGZF-compressed VCF (`.vcf.gz`) is supported through [`VcfReader::new_bgzf`] and
//! [`VcfWriter::new_bgzf`]. BCF and tabix indexes are not implemented yet.

use std::io::{BufRead, BufReader, Read, Write};

use crate::bgzf::{BgzfReader, BgzfWriter};
use crate::error::{Result, SeqError};

const FORMAT: &str = "VCF";
const FIXED_COLUMNS: [&str; 8] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

/// One `##key=value` meta line. Structured values (`<ID=x,...>`) are split into fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetaLine {
    /// The key before the first `=`, such as `"INFO"` or `"contig"`.
    pub key: String,
    /// The raw value after the first `=`.
    pub value: String,
    /// The `KEY=value` pairs of a structured value, with quotes removed. Empty when
    /// the value is not wrapped in `<...>`.
    pub fields: Vec<(String, String)>,
}

impl MetaLine {
    /// The value of a structured field, such as `ID` or `Number`.
    #[must_use]
    pub fn get(&self, key: &str) -> Option<&str> {
        self.fields.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }
}

/// A parsed VCF header. `text` is kept verbatim so output round-trips.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfHeader {
    text: String,
    meta: Vec<MetaLine>,
    samples: Vec<String>,
}

/// Splits the inside of a `<...>` value on commas that are outside quotes.
fn split_structured(inner: &str) -> Vec<&str> {
    let mut parts = Vec::new();
    let mut start = 0;
    let mut in_quotes = false;
    for (i, c) in inner.char_indices() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                parts.push(&inner[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    parts.push(&inner[start..]);
    parts
}

fn parse_structured(value: &str) -> Result<Vec<(String, String)>> {
    let Some(inner) = value.strip_prefix('<').and_then(|v| v.strip_suffix('>')) else {
        return Ok(Vec::new());
    };
    split_structured(inner)
        .into_iter()
        .map(|part| {
            let (k, v) = part.split_once('=').ok_or_else(|| {
                SeqError::malformed(FORMAT, format!("structured field without '=': {part}"))
            })?;
            Ok((k.to_string(), v.trim_matches('"').to_string()))
        })
        .collect()
}

impl VcfHeader {
    /// Parses the header lines (everything before the first record).
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] if the text lacks a `##fileformat` line or a
    /// `#CHROM` line, if the `#CHROM` columns are wrong, if a meta line has no `=`,
    /// or if a data line appears before `#CHROM`.
    pub fn parse(text: &str) -> Result<Self> {
        let mut meta = Vec::new();
        let mut samples = Vec::new();
        let mut saw_format = false;
        let mut saw_chrom = false;
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("##") {
                let (key, value) = rest.split_once('=').ok_or_else(|| {
                    SeqError::malformed(FORMAT, format!("meta line without '=': {line}"))
                })?;
                if key == "fileformat" {
                    if !value.starts_with("VCFv") {
                        return Err(SeqError::malformed(
                            FORMAT,
                            format!("unsupported fileformat: {value}"),
                        ));
                    }
                    saw_format = true;
                }
                meta.push(MetaLine {
                    key: key.to_string(),
                    value: value.to_string(),
                    fields: parse_structured(value)?,
                });
            } else if let Some(rest) = line.strip_prefix('#') {
                let cols: Vec<&str> = rest.split('\t').collect();
                if cols.len() < FIXED_COLUMNS.len() || cols[..FIXED_COLUMNS.len()] != FIXED_COLUMNS
                {
                    return Err(SeqError::malformed(FORMAT, "bad #CHROM column line"));
                }
                samples = cols[FIXED_COLUMNS.len()..]
                    .iter()
                    .skip(1) // the FORMAT column is not a sample
                    .map(|s| (*s).to_string())
                    .collect();
                saw_chrom = true;
            } else if !line.is_empty() {
                return Err(SeqError::malformed(FORMAT, "data line before #CHROM"));
            }
        }
        if !saw_format {
            return Err(SeqError::malformed(FORMAT, "missing ##fileformat line"));
        }
        if !saw_chrom {
            return Err(SeqError::malformed(FORMAT, "missing #CHROM line"));
        }
        let mut text = text.to_string();
        if !text.ends_with('\n') {
            text.push('\n');
        }
        Ok(Self { text, meta, samples })
    }

    /// The raw header text, as stored.
    #[must_use]
    pub fn text(&self) -> &str {
        &self.text
    }

    /// All meta lines, in file order.
    #[must_use]
    pub fn meta(&self) -> &[MetaLine] {
        &self.meta
    }

    /// The sample names from the `#CHROM` line.
    #[must_use]
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// The `fileformat` value, such as `"VCFv4.3"`.
    #[must_use]
    pub fn file_format(&self) -> Option<&str> {
        self.meta.iter().find(|m| m.key == "fileformat").map(|m| m.value.as_str())
    }

    /// Contigs as `(ID, length)`. The length is `None` when absent or not a number.
    pub fn contigs(&self) -> impl Iterator<Item = (&str, Option<u64>)> + '_ {
        self.lines_with_key("contig").filter_map(|m| {
            m.get("ID").map(|id| (id, m.get("length").and_then(|l| l.parse().ok())))
        })
    }

    /// `##INFO` definitions.
    pub fn info_defs(&self) -> impl Iterator<Item = &MetaLine> + '_ {
        self.lines_with_key("INFO")
    }

    /// `##FORMAT` definitions.
    pub fn format_defs(&self) -> impl Iterator<Item = &MetaLine> + '_ {
        self.lines_with_key("FORMAT")
    }

    /// `##FILTER` definitions.
    pub fn filter_defs(&self) -> impl Iterator<Item = &MetaLine> + '_ {
        self.lines_with_key("FILTER")
    }

    fn lines_with_key<'a>(&'a self, key: &'a str) -> impl Iterator<Item = &'a MetaLine> + 'a {
        self.meta.iter().filter(move |m| m.key == key)
    }
}

/// A zero-copy view of one VCF data line.
#[derive(Debug, Clone)]
pub struct VcfRecord<'a> {
    line: &'a str,
    cols: Vec<&'a str>,
    pos: u64,
    qual: Option<f32>,
}

impl<'a> VcfRecord<'a> {
    /// Parses one data line (no trailing newline) for a file with `n_samples` samples.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] with `line_no` for fewer than eight columns,
    /// a sample-column count that does not match `n_samples`, a POS that is not a
    /// positive integer, or a QUAL that is neither `.` nor a number.
    pub fn parse(line: &'a str, n_samples: usize, line_no: u64) -> Result<Self> {
        let fail =
            |reason: String| SeqError::malformed(FORMAT, format!("line {line_no}: {reason}"));
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < FIXED_COLUMNS.len() {
            return Err(fail(format!("expected at least 8 columns, found {}", cols.len())));
        }
        // With no FORMAT column there are no samples, so a file with samples needs one.
        let present = if cols.len() > FIXED_COLUMNS.len() { cols.len() - 9 } else { 0 };
        if cols.len() == FIXED_COLUMNS.len() && n_samples > 0 {
            return Err(fail(format!("header declares {n_samples} samples but record has none")));
        }
        if cols.len() > FIXED_COLUMNS.len() && present != n_samples {
            return Err(fail(format!(
                "header declares {n_samples} samples but record has {present}"
            )));
        }
        let pos: u64 = cols[1]
            .parse()
            .ok()
            .filter(|p| *p >= 1)
            .ok_or_else(|| fail(format!("POS is not a positive integer: {}", cols[1])))?;
        let qual = match cols[5] {
            "." => None,
            q => Some(q.parse::<f32>().map_err(|_| fail(format!("QUAL is not a number: {q}")))?),
        };
        Ok(Self { line, cols, pos, qual })
    }

    /// The whole line, without a newline.
    #[must_use]
    pub fn line(&self) -> &str {
        self.line
    }

    /// The chromosome or contig name.
    #[must_use]
    pub fn chrom(&self) -> &str {
        self.cols[0]
    }

    /// The 1-based position.
    #[must_use]
    pub fn pos(&self) -> u64 {
        self.pos
    }

    /// The identifier, or `None` for `.`.
    #[must_use]
    pub fn id(&self) -> Option<&str> {
        dot_to_none(self.cols[2])
    }

    /// The reference allele.
    #[must_use]
    pub fn reference(&self) -> &str {
        self.cols[3]
    }

    /// The alternate alleles. Empty for `.`.
    pub fn alts(&self) -> impl Iterator<Item = &str> + '_ {
        split_list(self.cols[4], ',')
    }

    /// The quality, or `None` for `.`.
    #[must_use]
    pub fn qual(&self) -> Option<f32> {
        self.qual
    }

    /// The filter names. `PASS` is returned as a name. Empty for `.`.
    pub fn filters(&self) -> impl Iterator<Item = &str> + '_ {
        split_list(self.cols[6], ';')
    }

    /// The INFO entries as `(key, value)`. Flags have no value.
    pub fn info(&self) -> impl Iterator<Item = (&str, Option<&str>)> + '_ {
        split_list(self.cols[7], ';').map(|kv| match kv.split_once('=') {
            Some((k, v)) => (k, Some(v)),
            None => (kv, None),
        })
    }

    /// The FORMAT keys, such as `GT` and `DP`. Empty when there are no samples.
    pub fn format_keys(&self) -> impl Iterator<Item = &str> + '_ {
        self.cols.get(8).map_or_else(|| split_list("", ':'), |f| split_list(f, ':'))
    }

    /// The raw sample columns, in header order. Each is a colon-separated list
    /// matching [`VcfRecord::format_keys`].
    pub fn samples(&self) -> impl Iterator<Item = &str> + '_ {
        self.cols.iter().skip(9).copied()
    }
}

fn dot_to_none(s: &str) -> Option<&str> {
    (s != ".").then_some(s)
}

/// Splits a list field on `sep`, treating `.` as an empty list.
fn split_list(s: &str, sep: char) -> impl Iterator<Item = &str> {
    let items = if s == "." || s.is_empty() { None } else { Some(s.split(sep)) };
    items.into_iter().flatten()
}

/// Reads VCF text. The header is parsed on construction.
#[derive(Debug)]
pub struct VcfReader<R> {
    inner: R,
    header: VcfHeader,
    pending: Option<String>,
    buf: String,
    line_no: u64,
}

impl<R: BufRead> VcfReader<R> {
    /// Reads the header and prepares to read records.
    ///
    /// # Errors
    /// Returns an I/O error, or [`SeqError::Malformed`] for an invalid header.
    pub fn new(mut inner: R) -> Result<Self> {
        let mut text = String::new();
        let mut pending = None;
        let mut line = String::new();
        let mut line_no = 0;
        loop {
            line.clear();
            if inner.read_line(&mut line)? == 0 {
                break;
            }
            line_no += 1;
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if trimmed.starts_with('#') {
                text.push_str(trimmed);
                text.push('\n');
            } else {
                if !trimmed.is_empty() {
                    pending = Some(trimmed.to_string());
                }
                break;
            }
        }
        let header = VcfHeader::parse(&text)?;
        Ok(Self { inner, header, pending, buf: String::new(), line_no })
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// Reads the next record, or `None` at end of input. Blank lines are skipped.
    /// The record borrows the reader's line buffer, so it must be dropped before
    /// the next call.
    ///
    /// # Errors
    /// Returns an I/O error, or [`SeqError::Malformed`] with the line number.
    pub fn next_record(&mut self) -> Result<Option<VcfRecord<'_>>> {
        if let Some(line) = self.pending.take() {
            self.buf = line;
        } else {
            loop {
                self.buf.clear();
                if self.inner.read_line(&mut self.buf)? == 0 {
                    return Ok(None);
                }
                self.line_no += 1;
                let len = self.buf.trim_end_matches(['\n', '\r']).len();
                self.buf.truncate(len);
                if !self.buf.is_empty() {
                    break;
                }
            }
        }
        let n_samples = self.header.samples().len();
        VcfRecord::parse(&self.buf, n_samples, self.line_no).map(Some)
    }
}

impl<S: Read> VcfReader<BufReader<BgzfReader<S>>> {
    /// Reads BGZF-compressed VCF (`.vcf.gz`).
    ///
    /// # Errors
    /// Same as [`VcfReader::new`], plus BGZF decoding errors.
    pub fn new_bgzf(inner: S) -> Result<VcfReader<BufReader<BgzfReader<S>>>> {
        VcfReader::new(BufReader::new(BgzfReader::new(inner)))
    }
}

/// Writes VCF text. The header is written on construction.
#[derive(Debug)]
pub struct VcfWriter<W> {
    inner: W,
    n_samples: usize,
}

impl<W: Write> VcfWriter<W> {
    /// Writes the header and returns a writer for records.
    ///
    /// # Errors
    /// Returns an I/O error.
    pub fn new(mut inner: W, header: &VcfHeader) -> Result<Self> {
        inner.write_all(header.text().as_bytes())?;
        Ok(Self { inner, n_samples: header.samples().len() })
    }

    /// Writes a record unchanged.
    ///
    /// # Errors
    /// Returns an I/O error.
    pub fn write_record(&mut self, record: &VcfRecord<'_>) -> Result<()> {
        writeln!(self.inner, "{}", record.line())?;
        Ok(())
    }

    /// Writes BGZF-compressed VCF. Call [`VcfWriter::finish`] to write the EOF block.
    ///
    /// # Errors
    /// Returns an I/O error.
    pub fn new_bgzf(inner: W, header: &VcfHeader) -> Result<VcfWriter<BgzfWriter<W>>> {
        VcfWriter::new(BgzfWriter::new(inner), header)
    }

    /// Validates the columns against the header and writes them as one line.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] if the columns would not parse, or an I/O error.
    pub fn write_columns(&mut self, columns: &[&str]) -> Result<()> {
        let line = columns.join("\t");
        // Parsing validates the same rules a reader enforces.
        VcfRecord::parse(&line, self.n_samples, 0)?;
        writeln!(self.inner, "{line}")?;
        Ok(())
    }

    /// Flushes and returns the underlying writer.
    ///
    /// # Errors
    /// Returns an I/O error from the flush.
    pub fn into_inner(mut self) -> Result<W> {
        self.inner.flush()?;
        Ok(self.inner)
    }
}

impl<W: Write> VcfWriter<BgzfWriter<W>> {
    /// Writes the BGZF EOF block, flushes, and returns the underlying writer.
    ///
    /// # Errors
    /// Returns an I/O error from the final block.
    pub fn finish(self) -> Result<W> {
        self.inner.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=248956422>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth, all samples\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

    #[test]
    fn parses_header_structure() {
        let h = VcfHeader::parse(HEADER).unwrap();
        assert_eq!(h.file_format(), Some("VCFv4.3"));
        assert_eq!(h.contigs().collect::<Vec<_>>(), [("chr1", Some(248_956_422))]);
        let dp = h.info_defs().next().unwrap();
        assert_eq!(dp.get("Description"), Some("Total depth, all samples"));
        assert_eq!(h.info_defs().count(), 2);
        assert_eq!(h.format_defs().count(), 1);
        assert_eq!(h.filter_defs().next().and_then(|f| f.get("ID")), Some("q10"));
        assert_eq!(h.samples(), ["s1", "s2"]);
    }

    #[test]
    fn rejects_bad_headers() {
        assert!(VcfHeader::parse("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").is_err());
        assert!(VcfHeader::parse("##fileformat=VCFv4.3\n").is_err());
        assert!(
            VcfHeader::parse("##fileformat=BCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                .is_err()
        );
        assert!(VcfHeader::parse("##fileformat=VCFv4.3\n#CHROM\tPOS\tREF\n").is_err());
    }

    #[test]
    fn parses_record_fields() {
        let line = "chr1\t100\trs1\tA\tG,T\t30.5\tPASS\tDP=12;DB\tGT:DP\t0/1:9\t.";
        let r = VcfRecord::parse(line, 2, 1).unwrap();
        assert_eq!(r.chrom(), "chr1");
        assert_eq!(r.pos(), 100);
        assert_eq!(r.id(), Some("rs1"));
        assert_eq!(r.alts().collect::<Vec<_>>(), ["G", "T"]);
        assert_eq!(r.qual(), Some(30.5));
        assert_eq!(r.filters().collect::<Vec<_>>(), ["PASS"]);
        assert_eq!(r.info().collect::<Vec<_>>(), [("DP", Some("12")), ("DB", None)]);
        assert_eq!(r.format_keys().collect::<Vec<_>>(), ["GT", "DP"]);
        assert_eq!(r.samples().collect::<Vec<_>>(), ["0/1:9", "."]);
    }

    #[test]
    fn dot_fields_are_empty() {
        let r = VcfRecord::parse("chr1\t5\t.\tC\t.\t.\t.\t.", 0, 1).unwrap();
        assert_eq!(r.id(), None);
        assert_eq!(r.alts().count(), 0);
        assert_eq!(r.qual(), None);
        assert_eq!(r.filters().count(), 0);
        assert_eq!(r.info().count(), 0);
        assert_eq!(r.format_keys().count(), 0);
    }

    #[test]
    fn rejects_bad_records_with_line_number() {
        let short = VcfRecord::parse("chr1\t5\t.\tC", 0, 7).unwrap_err().to_string();
        assert!(short.contains("line 7"), "{short}");
        assert!(VcfRecord::parse("chr1\t0\t.\tC\tG\t.\t.\t.", 0, 1).is_err());
        assert!(VcfRecord::parse("chr1\t5\t.\tC\tG\tabc\t.\t.", 0, 1).is_err());
        // Header declares two samples, record has one.
        assert!(VcfRecord::parse("chr1\t5\t.\tC\tG\t.\t.\t.\tGT\t0/1", 2, 1).is_err());
        // Header declares samples, record has no FORMAT column.
        assert!(VcfRecord::parse("chr1\t5\t.\tC\tG\t.\t.\t.", 2, 1).is_err());
    }

    #[test]
    fn reader_round_trips_through_writer() {
        // A blank line is skipped on read, so the written output omits it.
        let body = "\nchr1\t100\trs1\tA\tG\t30\tPASS\tDP=12\tGT\t0/1\t1/1\nchr1\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1\t0/0\n";
        let input = format!("{HEADER}{body}");
        let mut reader = VcfReader::new(Cursor::new(input.as_bytes())).unwrap();
        let header = reader.header().clone();
        let mut out = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut out, &header).unwrap();
            let mut count = 0;
            while let Some(rec) = reader.next_record().unwrap() {
                writer.write_record(&rec).unwrap();
                count += 1;
            }
            assert_eq!(count, 2);
        }
        let written = String::from_utf8(out).unwrap();
        assert_eq!(written, format!("{HEADER}{}", &body[1..]));
    }

    #[test]
    fn writer_validates_columns() {
        let header = VcfHeader::parse(HEADER).unwrap();
        let mut out = Vec::new();
        let mut writer = VcfWriter::new(&mut out, &header).unwrap();
        assert!(writer.write_columns(&["chr1", "5", ".", "C", "G", ".", ".", "."]).is_err());
        writer
            .write_columns(&["chr1", "5", ".", "C", "G", ".", ".", ".", "GT", "0/1", "1/1"])
            .unwrap();
    }

    #[test]
    fn bgzf_round_trip() {
        let body = "chr1	100	rs1	A	G	30	PASS	DP=12	GT	0/1	1/1
";
        let input = format!("{HEADER}{body}");
        let mut reader = VcfReader::new(Cursor::new(input.as_bytes())).unwrap();
        let header = reader.header().clone();
        let mut compressed = Vec::new();
        {
            let mut writer = VcfWriter::new_bgzf(&mut compressed, &header).unwrap();
            while let Some(rec) = reader.next_record().unwrap() {
                writer.write_record(&rec).unwrap();
            }
            writer.finish().unwrap();
        }
        // The output is compressed, so it must not contain the plain header.
        assert!(!compressed.windows(6).any(|w| w == b"#CHROM"));
        let mut reader = VcfReader::new_bgzf(Cursor::new(compressed)).unwrap();
        assert_eq!(reader.header().samples(), ["s1", "s2"]);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.pos(), 100);
        assert_eq!(rec.samples().collect::<Vec<_>>(), ["0/1", "1/1"]);
        assert!(reader.next_record().unwrap().is_none());
    }
}
