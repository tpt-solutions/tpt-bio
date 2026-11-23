//! SAM text reading and writing.

use std::fmt::Write as _;
use std::io::{BufRead, Write};

use crate::error::{Result, SeqError};
use crate::header::Header;
use crate::record::{CigarKind, CigarOp, Flags, Record, Tag, TagValue};

const FORMAT: &str = "SAM";

fn malformed(line_no: u64, reason: impl Into<String>) -> SeqError {
    SeqError::malformed(FORMAT, format!("line {line_no}: {}", reason.into()))
}

fn parse_int<T: std::str::FromStr>(field: &str, what: &str, line_no: u64) -> Result<T> {
    field.parse::<T>().map_err(|_| malformed(line_no, format!("{what} is not an integer: {field}")))
}

/// Parses one SAM alignment line (not a header line).
///
/// The field count is checked first, so the fixed-position indexing is in bounds.
///
/// # Errors
/// Returns [`SeqError::Malformed`] with the line number for any syntax error,
/// or for a reference name that is missing from `header`.
pub fn parse_record(line: &str, header: &Header, line_no: u64) -> Result<Record> {
    let fields: Vec<&str> = line.split('\t').collect();
    let [
        qname,
        flag,
        rname,
        pos_text,
        mapq_text,
        cigar_text,
        rnext,
        pnext,
        tlen_text,
        seq_text,
        qual_text,
        tag_fields @ ..,
    ] = fields.as_slice()
    else {
        return Err(malformed(
            line_no,
            format!("expected at least 11 fields, found {}", fields.len()),
        ));
    };
    let name = match *qname {
        "*" => Vec::new(),
        n => n.as_bytes().to_vec(),
    };
    let flags = Flags(parse_int::<u16>(flag, "FLAG", line_no)?);
    let ref_id = resolve_ref(rname, header, line_no)?;
    let pos = parse_int::<i64>(pos_text, "POS", line_no)? - 1;
    let mapq = parse_int::<u8>(mapq_text, "MAPQ", line_no)?;
    let cigar = parse_cigar(cigar_text, line_no)?;
    let next_ref_id = match *rnext {
        "=" => ref_id,
        other => resolve_ref(other, header, line_no)?,
    };
    let next_pos = parse_int::<i64>(pnext, "PNEXT", line_no)? - 1;
    let tlen = parse_int::<i64>(tlen_text, "TLEN", line_no)?;
    let seq = match *seq_text {
        "*" => Vec::new(),
        s => s.as_bytes().to_vec(),
    };
    let qual = match *qual_text {
        "*" => Vec::new(),
        q => q
            .bytes()
            .map(|b| {
                if (33..=126).contains(&b) {
                    Ok(b - 33)
                } else {
                    Err(malformed(line_no, "quality byte outside Phred+33 range"))
                }
            })
            .collect::<Result<Vec<u8>>>()?,
    };
    if !qual.is_empty() && qual.len() != seq.len() {
        return Err(malformed(line_no, "SEQ and QUAL lengths differ"));
    }
    let tags = tag_fields.iter().map(|f| parse_tag(f, line_no)).collect::<Result<Vec<_>>>()?;
    Ok(Record {
        name,
        flags,
        ref_id,
        pos,
        mapq,
        cigar,
        next_ref_id,
        next_pos,
        tlen,
        seq,
        qual,
        tags,
    })
}

fn resolve_ref(name: &str, header: &Header, line_no: u64) -> Result<i32> {
    if name == "*" {
        return Ok(-1);
    }
    header
        .ref_id(name)
        .ok_or_else(|| malformed(line_no, format!("reference {name} is not in the header")))
}

fn parse_cigar(text: &str, line_no: u64) -> Result<Vec<CigarOp>> {
    if text == "*" {
        return Ok(Vec::new());
    }
    let mut ops = Vec::new();
    let mut len: u64 = 0;
    let mut digits = 0;
    for b in text.bytes() {
        if b.is_ascii_digit() {
            len = len * 10 + u64::from(b - b'0');
            digits += 1;
            if len > u64::from(u32::MAX) {
                return Err(malformed(line_no, "CIGAR length overflows"));
            }
        } else {
            let kind = CigarKind::from_char(b)
                .ok_or_else(|| malformed(line_no, format!("bad CIGAR operation {}", b as char)))?;
            if digits == 0 {
                return Err(malformed(line_no, "CIGAR operation without length"));
            }
            ops.push(CigarOp { len: len as u32, kind });
            len = 0;
            digits = 0;
        }
    }
    if digits != 0 {
        return Err(malformed(line_no, "CIGAR ends with a length"));
    }
    Ok(ops)
}

fn parse_tag(field: &str, line_no: u64) -> Result<Tag> {
    let b = field.as_bytes();
    if b.len() < 5 || b[2] != b':' || b[4] != b':' {
        return Err(malformed(line_no, format!("bad tag {field}")));
    }
    let key = [b[0], b[1]];
    let body = &field[5..];
    let value = match b[3] {
        b'A' => {
            let c = body.as_bytes();
            if c.len() != 1 {
                return Err(malformed(
                    line_no,
                    format!("tag {field}: A value must be one character"),
                ));
            }
            TagValue::Char(c[0])
        }
        b'i' => TagValue::Int(parse_int::<i64>(body, "integer tag", line_no)?),
        b'f' => TagValue::Float(parse_int_float(body, line_no)?),
        b'Z' => TagValue::String(body.as_bytes().to_vec()),
        b'H' => TagValue::Hex(body.as_bytes().to_vec()),
        b'B' => parse_array(body, line_no)?,
        other => {
            return Err(malformed(line_no, format!("unknown tag type {}", other as char)));
        }
    };
    Ok(Tag { key, value })
}

fn parse_int_float(text: &str, line_no: u64) -> Result<f32> {
    text.parse::<f32>().map_err(|_| malformed(line_no, format!("bad float {text}")))
}

fn parse_array(body: &str, line_no: u64) -> Result<TagValue> {
    let mut parts = body.split(',');
    let subtype = parts.next().unwrap_or("");
    match subtype {
        "f" => {
            let values = parts
                .filter(|s| !s.is_empty())
                .map(|s| parse_int_float(s, line_no))
                .collect::<Result<Vec<f32>>>()?;
            Ok(TagValue::FloatArray(values))
        }
        "c" | "C" | "s" | "S" | "i" | "I" => {
            let values = parts
                .filter(|s| !s.is_empty())
                .map(|s| parse_int::<i64>(s, "array element", line_no))
                .collect::<Result<Vec<i64>>>()?;
            Ok(TagValue::IntArray(values))
        }
        other => Err(malformed(line_no, format!("unknown array subtype {other}"))),
    }
}

/// Formats one record as a SAM line (without a trailing newline).
///
/// # Errors
/// Returns [`SeqError::Malformed`] if a reference ID is outside the header.
pub fn format_record(rec: &Record, header: &Header) -> Result<String> {
    let mut out = String::with_capacity(128 + rec.seq.len() * 2);
    let name = |id: i32| -> Result<&str> {
        if id < 0 {
            Ok("*")
        } else {
            header.ref_name(id).ok_or_else(|| {
                SeqError::malformed(FORMAT, format!("reference id {id} is not in the header"))
            })
        }
    };
    let qname =
        if rec.name.is_empty() { "*" } else { std::str::from_utf8(&rec.name).unwrap_or("*") };
    let rname = name(rec.ref_id)?;
    let rnext = if rec.next_ref_id < 0 {
        "*".to_string()
    } else if rec.next_ref_id == rec.ref_id {
        "=".to_string()
    } else {
        name(rec.next_ref_id)?.to_string()
    };
    let seq = if rec.seq.is_empty() { "*" } else { std::str::from_utf8(&rec.seq).unwrap_or("*") };
    let qual = if rec.qual.is_empty() {
        "*".to_string()
    } else {
        rec.qual.iter().map(|&q| char::from(q + 33)).collect()
    };
    let _ = write!(
        out,
        "{qname}\t{}\t{rname}\t{}\t{}\t{}\t{rnext}\t{}\t{}\t{seq}\t{qual}",
        rec.flags.0,
        rec.pos + 1,
        rec.mapq,
        crate::record::format_cigar(&rec.cigar),
        rec.next_pos + 1,
        rec.tlen,
    );
    for tag in &rec.tags {
        out.push('\t');
        out.push(char::from(tag.key[0]));
        out.push(char::from(tag.key[1]));
        out.push(':');
        out.push_str(&tag.value.to_string());
    }
    Ok(out)
}

/// Streaming SAM reader. The header is read eagerly by [`SamReader::new`].
#[derive(Debug)]
pub struct SamReader<R> {
    inner: R,
    header: Header,
    pending: Option<String>,
    line_no: u64,
    buf: String,
}

impl<R: BufRead> SamReader<R> {
    /// Reads the header and prepares to read records.
    ///
    /// # Errors
    /// Returns an error on I/O failure or an invalid `@SQ` line.
    pub fn new(mut inner: R) -> Result<Self> {
        let mut header_text = String::new();
        let mut line_no = 0u64;
        let mut pending = None;
        let mut buf = String::new();
        loop {
            buf.clear();
            if inner.read_line(&mut buf)? == 0 {
                break;
            }
            line_no += 1;
            let line = buf.trim_end_matches(['\n', '\r']);
            if line.starts_with('@') {
                header_text.push_str(line);
                header_text.push('\n');
            } else {
                pending = Some(line.to_string());
                break;
            }
        }
        let header = Header::parse(&header_text)?;
        Ok(Self { inner, header, pending, line_no, buf: String::new() })
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The next record, or `None` at end of input.
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        if let Some(line) = self.pending.take() {
            return parse_record(&line, &self.header, self.line_no).map(Some);
        }
        loop {
            self.buf.clear();
            if self.inner.read_line(&mut self.buf)? == 0 {
                return Ok(None);
            }
            self.line_no += 1;
            let line = self.buf.trim_end_matches(['\n', '\r']);
            if line.is_empty() {
                continue;
            }
            return parse_record(line, &self.header, self.line_no).map(Some);
        }
    }
}

/// SAM writer. Writes the header on construction.
#[derive(Debug)]
pub struct SamWriter<W: Write> {
    inner: W,
    header: Header,
}

impl<W: Write> SamWriter<W> {
    /// Writes the header and returns a writer for records.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn new(mut inner: W, header: Header) -> Result<Self> {
        inner.write_all(header.text().as_bytes())?;
        Ok(Self { inner, header })
    }

    /// Writes one record.
    ///
    /// # Errors
    /// Returns an error if a reference ID is outside the header, or on I/O failure.
    pub fn write_record(&mut self, rec: &Record) -> Result<()> {
        let line = format_record(rec, &self.header)?;
        self.inner.write_all(line.as_bytes())?;
        self.inner.write_all(b"\n")?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::CigarOp;

    const SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n\
r1\t99\tchr1\t101\t60\t4M2I3M\t=\t201\t150\tACGTACGTAG\tIIIIIIIIII\tNM:i:2\tRG:Z:rg1\tXB:B:i,1,2,3\n\
r2\t4\t*\t0\t0\t*\t*\t0\t0\tACG\t*\tAS:i:-5\n\
r3\t16\tchr2\t10\t30\t3M\tchr1\t50\t0\tTTT\tFFF\tZZ:A:x\tXF:f:1.5\n";

    #[test]
    fn parse_fields() {
        let mut r = SamReader::new(SAM.as_bytes()).unwrap();
        assert_eq!(r.header().references().len(), 2);
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.name, b"r1");
        assert_eq!(rec.pos, 100);
        assert_eq!(rec.next_ref_id, rec.ref_id);
        assert_eq!(rec.qual[0], b'I' - 33);
        assert_eq!(rec.cigar[1], CigarOp { len: 2, kind: CigarKind::Insertion });
        assert_eq!(rec.tags.len(), 3);
        let unmapped = r.next_record().unwrap().unwrap();
        assert_eq!(unmapped.ref_id, -1);
        assert!(!unmapped.is_mapped());
        assert!(unmapped.qual.is_empty());
    }

    #[test]
    fn text_roundtrip() {
        let mut r = SamReader::new(SAM.as_bytes()).unwrap();
        let header = r.header().clone();
        let mut out = Vec::new();
        let mut w = SamWriter::new(&mut out, header).unwrap();
        while let Some(rec) = r.next_record().unwrap() {
            w.write_record(&rec).unwrap();
        }
        w.finish().unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), SAM);
    }

    #[test]
    fn errors_have_line_numbers() {
        let bad = "@SQ\tSN:chr1\tLN:10\nr1\t0\tchrX\t1\t0\t*\t*\t0\t0\tA\t*\n";
        let mut r = SamReader::new(bad.as_bytes()).unwrap();
        let err = r.next_record().unwrap_err().to_string();
        assert!(err.contains("line 2") && err.contains("chrX"), "{err}");
    }
}
