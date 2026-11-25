//! BAM binary reading and writing, with a zero-copy record view.
//!
//! BAM is BGZF-compressed. Each record is a little-endian `block_size`
//! followed by fixed fields, the read name, CIGAR, 4-bit packed sequence,
//! qualities and auxiliary tags. Field layout follows the `SAMv1` specification.

use std::io::{Read, Seek, Write};

use tpt_bio_core::Iupac;

use crate::bgzf::{BgzfReader, BgzfWriter, VirtualOffset, read_bounded};
use crate::error::{Result, SeqError};
use crate::header::{Header, Reference};
use crate::record::{CigarKind, CigarOp, Flags, Record, Tag, TagValue};

const FORMAT: &str = "BAM";
const MAGIC: &[u8; 4] = b"BAM\x01";
const FIXED_LEN: usize = 32;
/// Bin used for records with no reference position (`reg2bin(-1, 0)`).
const UNMAPPED_BIN: u16 = 4680;

fn bad(reason: impl Into<String>) -> SeqError {
    SeqError::malformed(FORMAT, reason)
}

/// Computes the UCSC binning bin for a 0-based half-open interval.
#[must_use]
pub fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let bin = if beg >> 14 == end >> 14 {
        4681 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        585 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        73 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        9 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        1 + (beg >> 26)
    } else {
        0
    };
    bin as u16
}

/// Zero-copy view over one BAM record body (the bytes after `block_size`).
#[derive(Debug, Clone, Copy)]
pub struct RecordView<'a> {
    data: &'a [u8],
    name_len: usize,
    n_cigar: usize,
    l_seq: usize,
}

impl<'a> RecordView<'a> {
    /// Validates the fixed layout and wraps the record body.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] if any length field overruns the body.
    pub fn new(data: &'a [u8]) -> Result<Self> {
        if data.len() < FIXED_LEN {
            return Err(bad("record shorter than fixed fields"));
        }
        let name_len = usize::from(data[8]);
        let n_cigar = usize::from(u16::from_le_bytes([data[12], data[13]]));
        let l_seq = le_i32(data, 16)?;
        let l_seq = usize::try_from(l_seq).map_err(|_| bad("negative l_seq"))?;
        let need = FIXED_LEN + name_len + 4 * n_cigar + l_seq.div_ceil(2) + l_seq;
        if need > data.len() {
            return Err(bad("record fields overrun the record block"));
        }
        Ok(Self { data, name_len, n_cigar, l_seq })
    }

    /// Reference index (-1 if unset).
    #[must_use]
    pub fn ref_id(&self) -> i32 {
        le_i32_unchecked(self.data, 0)
    }

    /// 0-based position (-1 if unset).
    #[must_use]
    pub fn pos(&self) -> i64 {
        i64::from(le_i32_unchecked(self.data, 4))
    }

    /// Mapping quality.
    #[must_use]
    pub fn mapq(&self) -> u8 {
        self.data[9]
    }

    /// Bin from the index.
    #[must_use]
    pub fn bin(&self) -> u16 {
        u16::from_le_bytes([self.data[10], self.data[11]])
    }

    /// FLAG bits.
    #[must_use]
    pub fn flags(&self) -> Flags {
        Flags(u16::from_le_bytes([self.data[14], self.data[15]]))
    }

    /// Mate reference index (-1 if unset).
    #[must_use]
    pub fn next_ref_id(&self) -> i32 {
        le_i32_unchecked(self.data, 20)
    }

    /// Mate 0-based position (-1 if unset).
    #[must_use]
    pub fn next_pos(&self) -> i64 {
        i64::from(le_i32_unchecked(self.data, 24))
    }

    /// Template length.
    #[must_use]
    pub fn tlen(&self) -> i64 {
        i64::from(le_i32_unchecked(self.data, 28))
    }

    /// Read name without the NUL terminator.
    #[must_use]
    pub fn name(&self) -> &'a [u8] {
        let raw = &self.data[FIXED_LEN..FIXED_LEN + self.name_len];
        raw.strip_suffix(&[0]).unwrap_or(raw)
    }

    fn cigar_start(&self) -> usize {
        FIXED_LEN + self.name_len
    }

    /// CIGAR operations.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for an unknown operation code.
    pub fn cigar(&self) -> Result<Vec<CigarOp>> {
        let mut ops = Vec::with_capacity(self.n_cigar);
        for i in 0..self.n_cigar {
            let v = le_u32_unchecked(self.data, self.cigar_start() + 4 * i);
            let kind = CigarKind::from_code((v & 0xf) as u8)
                .ok_or_else(|| bad("unknown CIGAR operation code"))?;
            ops.push(CigarOp { len: v >> 4, kind });
        }
        Ok(ops)
    }

    /// Reference bases spanned by the CIGAR, summed straight from the packed
    /// operations. Unlike [`RecordView::cigar`], this does not allocate, which
    /// matters when a region query tests many records it will not return.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for an unknown operation code.
    pub fn reference_span(&self) -> Result<u64> {
        let mut span = 0u64;
        for i in 0..self.n_cigar {
            let v = le_u32_unchecked(self.data, self.cigar_start() + 4 * i);
            let kind = CigarKind::from_code((v & 0xf) as u8)
                .ok_or_else(|| bad("unknown CIGAR operation code"))?;
            if kind.consumes_reference() {
                span += u64::from(v >> 4);
            }
        }
        Ok(span)
    }

    fn seq_start(&self) -> usize {
        self.cigar_start() + 4 * self.n_cigar
    }

    /// Bases as ASCII, uppercase.
    #[must_use]
    pub fn seq(&self) -> Vec<u8> {
        let packed = &self.data[self.seq_start()..self.seq_start() + self.l_seq.div_ceil(2)];
        let mut out = Vec::with_capacity(self.l_seq);
        for (i, &byte) in packed.iter().enumerate() {
            out.push(SEQ_ASCII[usize::from(byte >> 4)]);
            if 2 * i + 1 < self.l_seq {
                out.push(SEQ_ASCII[usize::from(byte & 0xf)]);
            }
        }
        out
    }

    fn qual_start(&self) -> usize {
        self.seq_start() + self.l_seq.div_ceil(2)
    }

    /// Raw Phred scores. Empty when the record stores `0xFF` (missing) qualities.
    #[must_use]
    pub fn qual(&self) -> Vec<u8> {
        let raw = &self.data[self.qual_start()..self.qual_start() + self.l_seq];
        if self.l_seq > 0 && raw[0] == 0xff { Vec::new() } else { raw.to_vec() }
    }

    /// The raw auxiliary-tag bytes.
    #[must_use]
    pub fn aux_bytes(&self) -> &'a [u8] {
        &self.data[self.qual_start() + self.l_seq..]
    }

    /// Decodes into an owned [`Record`].
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for invalid CIGAR codes or auxiliary data.
    pub fn to_record(&self) -> Result<Record> {
        let name = match self.name() {
            b"*" => Vec::new(),
            n => n.to_vec(),
        };
        Ok(Record {
            name,
            flags: self.flags(),
            ref_id: self.ref_id(),
            pos: self.pos(),
            mapq: self.mapq(),
            cigar: self.cigar()?,
            next_ref_id: self.next_ref_id(),
            next_pos: self.next_pos(),
            tlen: self.tlen(),
            seq: self.seq(),
            qual: self.qual(),
            tags: decode_tags(self.aux_bytes())?,
        })
    }
}

/// ASCII bases for the 4-bit BAM sequence codes, as `Iupac::to_ascii` gives them.
const SEQ_ASCII: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

fn le_i32(data: &[u8], at: usize) -> Result<i32> {
    data.get(at..at + 4)
        .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .ok_or_else(|| bad("truncated fixed field"))
}

fn le_i32_unchecked(data: &[u8], at: usize) -> i32 {
    i32::from_le_bytes([data[at], data[at + 1], data[at + 2], data[at + 3]])
}

fn le_u32_unchecked(data: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([data[at], data[at + 1], data[at + 2], data[at + 3]])
}

/// Smallest BAM integer type code that holds `v`, as used by the spec.
fn int_fits(code: u8, v: i64) -> bool {
    match code {
        b'c' => i8::try_from(v).is_ok(),
        b'C' => u8::try_from(v).is_ok(),
        b's' => i16::try_from(v).is_ok(),
        b'S' => u16::try_from(v).is_ok(),
        b'i' => i32::try_from(v).is_ok(),
        b'I' => u32::try_from(v).is_ok(),
        _ => false,
    }
}

/// The smallest integer type that holds every value in `min..=max`. Each BAM
/// integer type covers a contiguous range, so checking the extremes is enough.
fn int_code_for(min: i64, max: i64) -> Result<u8> {
    for code in [b'c', b'C', b's', b'S', b'i', b'I'] {
        if int_fits(code, min) && int_fits(code, max) {
            return Ok(code);
        }
    }
    Err(bad("integer tag value does not fit in 32 bits"))
}

fn int_size(code: u8) -> usize {
    match code {
        b'c' | b'C' => 1,
        b's' | b'S' => 2,
        _ => 4,
    }
}

fn push_int(out: &mut Vec<u8>, code: u8, v: i64) {
    match code {
        b'c' => out.push(v as i8 as u8),
        b'C' => out.push(v as u8),
        b's' => out.extend_from_slice(&(v as i16).to_le_bytes()),
        b'S' => out.extend_from_slice(&(v as u16).to_le_bytes()),
        b'i' => out.extend_from_slice(&(v as i32).to_le_bytes()),
        _ => out.extend_from_slice(&(v as u32).to_le_bytes()),
    }
}

fn encode_tags(tags: &[Tag], out: &mut Vec<u8>) -> Result<()> {
    for tag in tags {
        out.extend_from_slice(&tag.key);
        match &tag.value {
            TagValue::Char(c) => {
                out.push(b'A');
                out.push(*c);
            }
            TagValue::Int(v) => {
                let code = int_code_for(*v, *v)?;
                out.push(code);
                push_int(out, code, *v);
            }
            TagValue::Float(f) => {
                out.push(b'f');
                out.extend_from_slice(&f.to_le_bytes());
            }
            TagValue::String(s) => {
                out.push(b'Z');
                out.extend_from_slice(s);
                out.push(0);
            }
            TagValue::Hex(s) => {
                out.push(b'H');
                out.extend_from_slice(s);
                out.push(0);
            }
            TagValue::IntArray(values) => {
                let code = int_code_for(
                    values.iter().copied().min().unwrap_or(0),
                    values.iter().copied().max().unwrap_or(0),
                )?;
                out.push(b'B');
                out.push(code);
                out.extend_from_slice(&(values.len() as u32).to_le_bytes());
                for v in values {
                    push_int(out, code, *v);
                }
            }
            TagValue::FloatArray(values) => {
                out.push(b'B');
                out.push(b'f');
                out.extend_from_slice(&(values.len() as u32).to_le_bytes());
                for v in values {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
    }
    Ok(())
}

fn take<'a>(data: &'a [u8], at: &mut usize, n: usize) -> Result<&'a [u8]> {
    let end = at.checked_add(n).ok_or_else(|| bad("tag length overflow"))?;
    let slice = data.get(*at..end).ok_or_else(|| bad("truncated auxiliary data"))?;
    *at = end;
    Ok(slice)
}

fn take_nul(data: &[u8], at: &mut usize) -> Result<Vec<u8>> {
    let rest = &data[*at..];
    let end = rest.iter().position(|&b| b == 0).ok_or_else(|| bad("unterminated string tag"))?;
    *at += end + 1;
    Ok(rest[..end].to_vec())
}

fn read_int(code: u8, bytes: &[u8]) -> Result<i64> {
    Ok(match code {
        b'c' => i64::from(bytes[0] as i8),
        b'C' => i64::from(bytes[0]),
        b's' => i64::from(i16::from_le_bytes([bytes[0], bytes[1]])),
        b'S' => i64::from(u16::from_le_bytes([bytes[0], bytes[1]])),
        b'i' => i64::from(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        b'I' => i64::from(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        _ => return Err(bad("unknown integer tag type")),
    })
}

/// Decodes the auxiliary block of a record.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for truncated or unknown tag types.
pub fn decode_tags(data: &[u8]) -> Result<Vec<Tag>> {
    let mut tags = Vec::new();
    let mut at = 0usize;
    while at < data.len() {
        let key_bytes = take(data, &mut at, 2)?;
        let key = [key_bytes[0], key_bytes[1]];
        let ty = take(data, &mut at, 1)?[0];
        let value = match ty {
            b'A' => TagValue::Char(take(data, &mut at, 1)?[0]),
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
                let v = read_int(ty, take(data, &mut at, int_size(ty))?)?;
                TagValue::Int(v)
            }
            b'f' => {
                let b = take(data, &mut at, 4)?;
                TagValue::Float(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            }
            b'Z' => TagValue::String(take_nul(data, &mut at)?),
            b'H' => TagValue::Hex(take_nul(data, &mut at)?),
            b'B' => {
                let sub = take(data, &mut at, 1)?[0];
                let count_bytes = take(data, &mut at, 4)?;
                let count = u32::from_le_bytes([
                    count_bytes[0],
                    count_bytes[1],
                    count_bytes[2],
                    count_bytes[3],
                ]) as usize;
                if sub == b'f' {
                    let raw = take(data, &mut at, 4 * count)?;
                    TagValue::FloatArray(
                        raw.chunks_exact(4)
                            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect(),
                    )
                } else {
                    let size = int_size(sub);
                    let raw = take(data, &mut at, size * count)?;
                    let values =
                        raw.chunks_exact(size).map(|c| read_int(sub, c)).collect::<Result<_>>()?;
                    TagValue::IntArray(values)
                }
            }
            other => return Err(bad(format!("unknown tag type {}", other as char))),
        };
        tags.push(Tag { key, value });
    }
    Ok(tags)
}

/// Encodes one record into `out` (appending), including its `block_size` prefix.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for bases outside IUPAC, a quality length that
/// does not match the sequence, or values that do not fit the BAM field widths.
pub fn encode_record(rec: &Record, out: &mut Vec<u8>) -> Result<()> {
    if !rec.qual.is_empty() && rec.qual.len() != rec.seq.len() {
        return Err(bad("sequence and quality lengths differ"));
    }
    let start = out.len();
    out.extend_from_slice(&[0u8; 4]); // block_size, patched below
    let end = rec.alignment_end().unwrap_or(rec.pos + 1).max(rec.pos + 1);
    let bin = if rec.pos < 0 { UNMAPPED_BIN } else { reg2bin(rec.pos, end) };

    let name: &[u8] = if rec.name.is_empty() { b"*" } else { &rec.name };
    let name_len = name.len() + 1;
    if name_len > 255 {
        return Err(bad("read name longer than 254 bytes"));
    }
    if rec.cigar.len() > usize::from(u16::MAX) {
        return Err(bad("too many CIGAR operations"));
    }
    let to_i32 = |v: i64, what: &str| -> Result<i32> {
        i32::try_from(v).map_err(|_| bad(format!("{what} does not fit in 32 bits")))
    };
    let pos = to_i32(rec.pos, "position")?;
    let next_pos = to_i32(rec.next_pos, "mate position")?;
    let tlen = to_i32(rec.tlen, "template length")?;

    out.extend_from_slice(&rec.ref_id.to_le_bytes());
    out.extend_from_slice(&pos.to_le_bytes());
    out.push(name_len as u8);
    out.push(rec.mapq);
    out.extend_from_slice(&bin.to_le_bytes());
    out.extend_from_slice(&(rec.cigar.len() as u16).to_le_bytes());
    out.extend_from_slice(&rec.flags.0.to_le_bytes());
    out.extend_from_slice(&(rec.seq.len() as u32).to_le_bytes());
    out.extend_from_slice(&rec.next_ref_id.to_le_bytes());
    out.extend_from_slice(&next_pos.to_le_bytes());
    out.extend_from_slice(&tlen.to_le_bytes());

    out.extend_from_slice(name);
    out.push(0);
    for op in &rec.cigar {
        out.extend_from_slice(&((op.len << 4) | u32::from(op.kind.to_code())).to_le_bytes());
    }
    for pair in rec.seq.chunks(2) {
        let hi = iupac_nibble(pair[0])?;
        let lo = pair.get(1).map_or(Ok(0), |&b| iupac_nibble(b))?;
        out.push((hi << 4) | lo);
    }
    if rec.qual.is_empty() {
        out.extend(std::iter::repeat_n(0xffu8, rec.seq.len()));
    } else {
        out.extend_from_slice(&rec.qual);
    }
    encode_tags(&rec.tags, out)?;

    let block_size = (out.len() - start - 4) as u32;
    out[start..start + 4].copy_from_slice(&block_size.to_le_bytes());
    Ok(())
}

fn iupac_nibble(b: u8) -> Result<u8> {
    Iupac::from_ascii(b)
        .map(Iupac::to_nibble)
        .ok_or_else(|| bad(format!("invalid base 0x{b:02X} in sequence")))
}

/// Reads a BAM stream: header first, then records.
#[derive(Debug)]
pub struct BamReader<R: Read> {
    inner: BgzfReader<R>,
    header: Header,
    body: Vec<u8>,
}

impl<R: Read> BamReader<R> {
    /// Reads and validates the BAM header.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a bad magic number or truncated header.
    pub fn new(inner: R) -> Result<Self> {
        let mut inner = BgzfReader::new(inner);
        let mut magic = [0u8; 4];
        read_exact(&mut inner, &mut magic)?;
        if &magic != MAGIC {
            return Err(bad("missing BAM magic"));
        }
        let text_len = read_u32(&mut inner)? as usize;
        let text = read_vec(&mut inner, text_len)?;
        let text = String::from_utf8(strip_nul(&text).to_vec())
            .map_err(|_| bad("header text is not UTF-8"))?;
        let n_ref = read_u32(&mut inner)?;
        let mut references = Vec::new();
        for _ in 0..n_ref {
            let l_name = read_u32(&mut inner)? as usize;
            let name = read_vec(&mut inner, l_name)?;
            let name = String::from_utf8(strip_nul(&name).to_vec())
                .map_err(|_| bad("reference name is not UTF-8"))?;
            let len = u64::from(read_u32(&mut inner)?);
            references.push(Reference { name, len });
        }
        let header = Header::from_text_and_references(&text, references)?;
        Ok(Self { inner, header, body: Vec::new() })
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The virtual offset of the next record's first byte.
    #[must_use]
    pub fn virtual_position(&self) -> VirtualOffset {
        self.inner.virtual_position()
    }

    /// Reads the next record into `self.body`, returning a view over it.
    ///
    /// The view is valid until the next call. Use [`RecordView::to_record`] to keep it.
    pub fn next_view(&mut self) -> Result<Option<RecordView<'_>>> {
        Ok(self.next_view_with_offsets()?.map(|(_, _, view)| view))
    }

    /// Like [`BamReader::next_view`], also returning the virtual offsets before
    /// and after the record. Index building needs both.
    pub fn next_view_with_offsets(
        &mut self,
    ) -> Result<Option<(VirtualOffset, VirtualOffset, RecordView<'_>)>> {
        let begin = self.inner.virtual_position();
        let mut size_bytes = [0u8; 4];
        let n = read_up_to(&mut self.inner, &mut size_bytes)?;
        if n == 0 {
            return Ok(None);
        }
        if n < 4 {
            return Err(bad("truncated record size"));
        }
        let size = u32::from_le_bytes(size_bytes) as usize;
        if size <= INLINE_RECORD_BYTES {
            self.body.resize(size, 0);
            read_exact(&mut self.inner, &mut self.body)?;
        } else {
            self.body = read_vec(&mut self.inner, size)?;
        }
        let end = self.inner.virtual_position();
        Ok(Some((begin, end, RecordView::new(&self.body)?)))
    }

    /// Reads the next record as an owned [`Record`].
    /// Reads up to `max` raw record bodies (the bytes after each block size) and appends
    /// them to `out`. Returns how many were read, which is 0 at the end of the file.
    ///
    /// Each body can be parsed with [`RecordView::new`] on another thread, which is how
    /// [`crate::sam::SamWriter::write_bam_batch`] formats records in parallel.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a truncated record, and I/O errors.
    pub fn read_raw_batch(&mut self, max: usize, out: &mut Vec<Vec<u8>>) -> Result<usize> {
        let mut read = 0;
        while read < max {
            let mut size_bytes = [0u8; 4];
            let n = read_up_to(&mut self.inner, &mut size_bytes)?;
            if n == 0 {
                break;
            }
            if n < 4 {
                return Err(bad("truncated record size"));
            }
            let size = u32::from_le_bytes(size_bytes) as usize;
            out.push(read_vec(&mut self.inner, size)?);
            read += 1;
        }
        Ok(read)
    }

    pub fn next_record(&mut self) -> Result<Option<Record>> {
        match self.next_view()? {
            Some(view) => view.to_record().map(Some),
            None => Ok(None),
        }
    }
}

impl<R: Read + Seek> BamReader<R> {
    /// Repositions to a virtual offset, such as a chunk start from a BAI index.
    pub fn seek_virtual(&mut self, voffset: VirtualOffset) -> Result<()> {
        self.inner.seek_virtual(voffset)
    }
}

/// Writes a BAM stream. The header is written on construction.
#[derive(Debug)]
pub struct BamWriter<W: Write> {
    inner: BgzfWriter<W>,
    header: Header,
    scratch: Vec<u8>,
}

impl<W: Write> BamWriter<W> {
    /// Writes the header and returns a writer for records.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn new(inner: W, header: Header) -> Result<Self> {
        let mut inner = BgzfWriter::new(inner);
        let mut head = Vec::new();
        head.extend_from_slice(MAGIC);
        head.extend_from_slice(&(header.text().len() as u32).to_le_bytes());
        head.extend_from_slice(header.text().as_bytes());
        head.extend_from_slice(&(header.references().len() as u32).to_le_bytes());
        for r in header.references() {
            head.extend_from_slice(&((r.name.len() + 1) as u32).to_le_bytes());
            head.extend_from_slice(r.name.as_bytes());
            head.push(0);
            head.extend_from_slice(&(r.len as u32).to_le_bytes());
        }
        inner.write_all(&head)?;
        Ok(Self { inner, header, scratch: Vec::new() })
    }

    /// Encodes and writes one record.
    ///
    /// # Errors
    /// Returns an error if the record cannot be encoded or on I/O failure.
    pub fn write_record(&mut self, rec: &Record) -> Result<()> {
        if rec.ref_id >= 0 && self.header.ref_name(rec.ref_id).is_none() {
            return Err(bad(format!("reference id {} is not in the header", rec.ref_id)));
        }
        self.scratch.clear();
        encode_record(rec, &mut self.scratch)?;
        self.inner.write_all(&self.scratch)?;
        Ok(())
    }

    /// Writes the EOF marker and returns the inner writer.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn finish(self) -> Result<W> {
        self.inner.finish()
    }
}

fn strip_nul(bytes: &[u8]) -> &[u8] {
    bytes.strip_suffix(&[0]).unwrap_or(bytes)
}

/// Records up to this size reuse the record buffer; larger ones are read with
/// [`read_vec`] so a corrupt size field cannot force a large zero-filled allocation.
const INLINE_RECORD_BYTES: usize = 1 << 16;

/// Reads `len` bytes that were length-prefixed in the file. See [`read_bounded`].
fn read_vec<R: Read>(r: &mut R, len: usize) -> Result<Vec<u8>> {
    read_bounded(r, len).map_err(|e| match e.kind() {
        std::io::ErrorKind::UnexpectedEof => bad("truncated file"),
        _ => SeqError::Io(e),
    })
}

fn read_exact<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<()> {
    r.read_exact(buf).map_err(|e| match e.kind() {
        std::io::ErrorKind::UnexpectedEof => bad("truncated file"),
        _ => SeqError::Io(e),
    })
}

fn read_up_to<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    read_exact(r, &mut b)?;
    Ok(u32::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    #[test]
    fn seq_table_matches_iupac_for_every_nibble() {
        for v in 0..16u8 {
            let iupac = Iupac::from_nibble(v).expect("every nibble is a code");
            assert_eq!(SEQ_ASCII[usize::from(v)], iupac.to_ascii(), "nibble {v}");
        }
    }

    use super::*;
    use crate::record::{CigarKind, CigarOp};
    use crate::sam::SamReader;
    use std::io::Cursor;

    const SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n\
r1\t99\tchr1\t101\t60\t4M2I3M\t=\t201\t150\tACGTACGTAG\tIIIIIIIIII\tNM:i:2\tRG:Z:rg1\tXB:B:i,1,2,3\n\
r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGNN\t*\tAS:i:-5\tXH:H:1AE3\n\
r3\t16\tchr2\t10\t30\t3M\tchr1\t50\t0\tTTT\tFFF\tZZ:A:x\tXF:f:1.5\tXG:B:f,0.5,2\tXL:i:70000\tXM:i:-70000\n";

    #[test]
    fn reg2bin_known_values() {
        assert_eq!(reg2bin(-1, 0), UNMAPPED_BIN);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 16384), 4681);
        assert_eq!(reg2bin(0, 16385), 585);
        assert_eq!(reg2bin(0, 16400), 585);
    }

    #[test]
    fn sam_bam_sam_roundtrip() {
        let mut sam = SamReader::new(SAM.as_bytes()).unwrap();
        let header = sam.header().clone();
        let mut bam_bytes = Vec::new();
        let mut w = BamWriter::new(&mut bam_bytes, header.clone()).unwrap();
        let mut records = Vec::new();
        while let Some(rec) = sam.next_record().unwrap() {
            w.write_record(&rec).unwrap();
            records.push(rec);
        }
        w.finish().unwrap();

        let mut r = BamReader::new(Cursor::new(&bam_bytes)).unwrap();
        assert_eq!(r.header().text(), header.text());
        assert_eq!(r.header().references(), header.references());
        let mut back = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            back.push(rec);
        }
        assert_eq!(back, records);
    }

    #[test]
    fn view_exposes_fields_without_copying() {
        let mut sam = SamReader::new(SAM.as_bytes()).unwrap();
        let header = sam.header().clone();
        let rec = sam.next_record().unwrap().unwrap();
        let mut buf = Vec::new();
        encode_record(&rec, &mut buf).unwrap();
        let view = RecordView::new(&buf[4..]).unwrap();
        assert_eq!(view.name(), b"r1");
        assert_eq!(view.pos(), 100);
        assert_eq!(view.seq(), b"ACGTACGTAG");
        assert_eq!(view.qual(), vec![40u8; 10]);
        assert_eq!(view.cigar().unwrap()[1], CigarOp { len: 2, kind: CigarKind::Insertion });
        assert_eq!(header.ref_name(view.ref_id()), Some("chr1"));
    }

    #[test]
    fn integer_tags_choose_smallest_type() {
        let mut out = Vec::new();
        encode_tags(&[Tag { key: *b"XX", value: TagValue::Int(-70000) }], &mut out).unwrap();
        assert_eq!(out[2], b'i');
        let back = decode_tags(&out).unwrap();
        assert_eq!(back[0].value, TagValue::Int(-70000));
    }

    #[test]
    fn truncated_records_are_errors_not_panics() {
        let mut sam = SamReader::new(SAM.as_bytes()).unwrap();
        let rec = sam.next_record().unwrap().unwrap();
        let mut buf = Vec::new();
        encode_record(&rec, &mut buf).unwrap();
        for cut in 4..buf.len() {
            let _ = RecordView::new(&buf[4..cut]).and_then(|v| v.to_record());
        }
    }
}
