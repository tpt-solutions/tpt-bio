//! CRAM 3.0 writer.
//!
//! [`CramWriter`] writes the file definition, a SAM header container, data
//! containers and an end-of-file container. Each container holds one reference
//! id (a container is closed when the reference changes) and up to
//! `slices_per_container` slices of `records_per_slice` records.
//!
//! Layout choices, all within the spec's data model:
//! - Every data series is an EXTERNAL block of its own; byte arrays use a
//!   length/bytes pair, read names use a 0-terminated stop array. The core block
//!   is empty, so no bit codes are written.
//! - Alignment starts are absolute (AP delta off). The CRAM flags are always
//!   0x3: qualities stored as an array, and mate fields stored detached.
//! - Mapped records are encoded against the reference: mismatches become `X`
//!   features and CIGAR insertions, deletions, skips and clips become their
//!   features. Without a reference, each aligned run is a `b` feature carrying its
//!   bases, and the reference-required flag is cleared.
//! - Slices carry the MD5 of their reference region when a reference is supplied.
//! - Gzip compresses data blocks. Header blocks are stored raw.
//!
//! The end-of-file marker is a self-consistent empty container (reference -1,
//! alignment start 4542278, no records) with a minimal compression header. The
//! 38-byte marker that samtools writes was not reproduced from the specification
//! text here, so the bytes are not claimed to match it.

use std::collections::HashMap;
use std::io::Write;

use crate::cram::block::{self, Method, content};
use crate::cram::bytes::put_itf8;
use crate::cram::compression::{self, CompressionHeader, Encoding, Preservation};
use crate::cram::container::{self, ContainerHeader};
use crate::cram::decode::substitution_code;
use crate::cram::reference::CramReference;
use crate::cram::slice::{self, SliceHeader};
use crate::error::{Result, SeqError};
use crate::header::Header;
use crate::record::{CigarKind, Flags, Record, TagValue, query_length, reference_length};

const FORMAT: &str = "CRAM writer";
/// Substitution matrix written to the preservation map: alternative `i` has code `i`.
const SUBSTITUTION_MATRIX: [u8; 5] = [0x1B; 5];
/// The CRAM 3.0 end-of-file container, byte for byte as samtools and htslib write
/// it. Readers check for it, so it is written verbatim rather than built.
const EOF_CONTAINER_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, // container length 15
    0xff, 0xff, 0xff, 0xff, // reference id -1
    0x0f, // alignment start 4542278 ("EOF") as ITF8
    0xe0, 0x45, 0x4f, 0x46, //
    0x00, 0x00, 0x00, 0x00, // span, record count, record counter
    0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, // bases, block count, landmarks, crc32
    0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, // empty compression header
    0xee, 0x63, 0x01, 0x4b, // block crc32
];

/// How a data series is stored.
#[derive(Clone, Copy)]
enum Kind {
    /// Integers or single bytes in one external block.
    Ext,
    /// Byte array: length and bytes in two external blocks.
    Arr,
    /// Byte array terminated by a stop byte in one external block.
    Stop(u8),
}

/// Assigns external block ids and encodings for one container.
struct Plan {
    next_id: i32,
    series_ids: HashMap<[u8; 2], (i32, i32)>,
    tag_ids: HashMap<[u8; 3], (i32, i32)>,
    series: Vec<([u8; 2], Encoding)>,
    tags: Vec<([u8; 3], Encoding)>,
    ext_ids: Vec<i32>,
    tag_lines: Vec<Vec<[u8; 3]>>,
    tag_line_index: HashMap<Vec<[u8; 3]>, usize>,
}

impl Plan {
    fn new() -> Self {
        let mut p = Self {
            next_id: 1,
            series_ids: HashMap::new(),
            tag_ids: HashMap::new(),
            series: Vec::new(),
            tags: Vec::new(),
            ext_ids: Vec::new(),
            tag_lines: Vec::new(),
            tag_line_index: HashMap::new(),
        };
        let fixed: &[(&[u8; 2], Kind)] = &[
            (b"BF", Kind::Ext),
            (b"CF", Kind::Ext),
            (b"RL", Kind::Ext),
            (b"AP", Kind::Ext),
            (b"RG", Kind::Ext),
            (b"RN", Kind::Stop(0)),
            (b"MF", Kind::Ext),
            (b"NS", Kind::Ext),
            (b"NP", Kind::Ext),
            (b"TS", Kind::Ext),
            (b"TL", Kind::Ext),
            (b"FN", Kind::Ext),
            (b"FC", Kind::Ext),
            (b"FP", Kind::Ext),
            (b"BS", Kind::Ext),
            (b"IN", Kind::Arr),
            (b"SC", Kind::Arr),
            (b"DL", Kind::Ext),
            (b"RS", Kind::Ext),
            (b"PD", Kind::Ext),
            (b"HC", Kind::Ext),
            (b"BA", Kind::Ext),
            (b"QS", Kind::Ext),
            (b"MQ", Kind::Ext),
            (b"BB", Kind::Arr),
        ];
        for (key, kind) in fixed {
            p.series_key(**key, *kind);
        }
        p
    }

    fn alloc(&mut self) -> i32 {
        let id = self.next_id;
        self.next_id += 1;
        self.ext_ids.push(id);
        id
    }

    fn series_key(&mut self, key: [u8; 2], kind: Kind) {
        let (a, b, enc) = match kind {
            Kind::Ext => {
                let a = self.alloc();
                (a, -1, Encoding::External { content_id: a })
            }
            Kind::Arr => {
                let a = self.alloc();
                let b = self.alloc();
                (
                    a,
                    b,
                    Encoding::ByteArrayLen {
                        len: Box::new(Encoding::External { content_id: a }),
                        bytes: Box::new(Encoding::External { content_id: b }),
                    },
                )
            }
            Kind::Stop(stop) => {
                let a = self.alloc();
                (a, -1, Encoding::ByteArrayStop { stop, content_id: a })
            }
        };
        self.series_ids.insert(key, (a, b));
        self.series.push((key, enc));
    }

    /// Allocates a tag's byte-array blocks on first use.
    fn tag_key(&mut self, key: [u8; 3]) {
        if self.tag_ids.contains_key(&key) {
            return;
        }
        let a = self.alloc();
        let b = self.alloc();
        self.tag_ids.insert(key, (a, b));
        self.tags.push((
            key,
            Encoding::ByteArrayLen {
                len: Box::new(Encoding::External { content_id: a }),
                bytes: Box::new(Encoding::External { content_id: b }),
            },
        ));
    }

    /// Returns the tag line index for a list of tag keys, adding it if new.
    fn tag_line(&mut self, keys: Vec<[u8; 3]>) -> usize {
        if let Some(&i) = self.tag_line_index.get(&keys) {
            return i;
        }
        let i = self.tag_lines.len();
        self.tag_line_index.insert(keys.clone(), i);
        self.tag_lines.push(keys);
        i
    }

    fn compression_header(&self, reference: bool) -> CompressionHeader {
        CompressionHeader {
            preservation: Preservation {
                read_names: true,
                ap_delta: false,
                reference_required: reference,
                substitution_matrix: SUBSTITUTION_MATRIX,
                tag_dictionary: self.tag_lines.clone(),
            },
            data_series: self.series.clone(),
            tags: self.tags.clone(),
        }
    }
}

/// Byte buffers for each external block, plus the plan that names them.
struct Encoder {
    plan: Plan,
    bufs: HashMap<i32, Vec<u8>>,
}

impl Encoder {
    fn buf(&mut self, id: i32) -> &mut Vec<u8> {
        self.bufs.entry(id).or_default()
    }

    fn int(&mut self, key: [u8; 2], v: i32) -> Result<()> {
        let (id, _) = self.ids(key)?;
        put_itf8(self.buf(id), v);
        Ok(())
    }

    fn byte(&mut self, key: [u8; 2], v: u8) -> Result<()> {
        let (id, _) = self.ids(key)?;
        self.buf(id).push(v);
        Ok(())
    }

    fn array(&mut self, key: [u8; 2], data: &[u8]) -> Result<()> {
        let (len_id, bytes_id) = self.ids(key)?;
        put_itf8(self.buf(len_id), data.len() as i32);
        self.buf(bytes_id).extend_from_slice(data);
        Ok(())
    }

    fn stop_array(&mut self, key: [u8; 2], data: &[u8], stop: u8) -> Result<()> {
        let (id, _) = self.ids(key)?;
        let b = self.buf(id);
        b.extend_from_slice(data);
        b.push(stop);
        Ok(())
    }

    fn ids(&self, key: [u8; 2]) -> Result<(i32, i32)> {
        self.plan
            .series_ids
            .get(&key)
            .copied()
            .ok_or_else(|| SeqError::malformed(FORMAT, "unknown data series"))
    }

    fn tag(&mut self, key: [u8; 3], value: &[u8]) {
        let (len_id, bytes_id) = self.plan.tag_ids[&key];
        put_itf8(self.buf(len_id), value.len() as i32);
        self.buf(bytes_id).extend_from_slice(value);
    }
}

/// A read feature to write, with its 1-based read position.
enum FeatureData {
    Substitution(u8),
    Bases(Vec<u8>),
    Insertion(Vec<u8>),
    SoftClip(Vec<u8>),
    Deletion(usize),
    RefSkip(usize),
    HardClip(usize),
    Padding(usize),
}

fn checked_i32(v: i64, what: &str) -> Result<i32> {
    i32::try_from(v).map_err(|_| {
        SeqError::Unsupported(format!("{what} {v} does not fit in CRAM's 32-bit field"))
    })
}

/// Tag value bytes in BAM binary form (without the key), and the CRAM type code.
fn tag_value(v: &TagValue) -> Result<(u8, Vec<u8>)> {
    Ok(match v {
        TagValue::Char(c) => (b'A', vec![*c]),
        TagValue::Int(i) => (b'i', checked_i32(*i, "integer tag")?.to_le_bytes().to_vec()),
        TagValue::Float(f) => (b'f', f.to_le_bytes().to_vec()),
        TagValue::String(s) => {
            let mut out = s.clone();
            out.push(0);
            (b'Z', out)
        }
        TagValue::Hex(s) => {
            let mut out = s.clone();
            out.push(0);
            (b'H', out)
        }
        TagValue::IntArray(values) => {
            let mut out = vec![b'i'];
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for &x in values {
                out.extend_from_slice(&checked_i32(x, "integer array element")?.to_le_bytes());
            }
            (b'B', out)
        }
        TagValue::FloatArray(values) => {
            let mut out = vec![b'f'];
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for x in values {
                out.extend_from_slice(&x.to_le_bytes());
            }
            (b'B', out)
        }
    })
}

/// Builds the read features of a mapped record from its CIGAR, bases and reference.
fn features_of(rec: &Record, window: Option<&[u8]>) -> Result<Vec<(u8, usize, FeatureData)>> {
    let seq = &rec.seq;
    let mut out = Vec::new();
    let mut rp = 1usize; // next 1-based read position
    let mut ro = 0usize; // next reference index in `window`
    for op in &rec.cigar {
        let len = op.len as usize;
        match op.kind {
            CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch => {
                if let Some(w) = window {
                    for _ in 0..len {
                        let read_b = seq[rp - 1].to_ascii_uppercase();
                        let ref_b = *w
                            .get(ro)
                            .ok_or_else(|| SeqError::malformed(FORMAT, "reference too short"))?;
                        if read_b != ref_b {
                            let code = substitution_code(&SUBSTITUTION_MATRIX, ref_b, read_b)
                                .ok_or_else(|| {
                                    SeqError::Unsupported(format!(
                                        "base {} cannot be written as a substitution",
                                        read_b as char
                                    ))
                                })?;
                            out.push((b'X', rp, FeatureData::Substitution(code)));
                        }
                        rp += 1;
                        ro += 1;
                    }
                } else {
                    out.push((b'b', rp, FeatureData::Bases(seq[rp - 1..rp - 1 + len].to_vec())));
                    rp += len;
                }
            }
            CigarKind::Insertion => {
                out.push((b'I', rp, FeatureData::Insertion(seq[rp - 1..rp - 1 + len].to_vec())));
                rp += len;
            }
            CigarKind::SoftClip => {
                out.push((b'S', rp, FeatureData::SoftClip(seq[rp - 1..rp - 1 + len].to_vec())));
                rp += len;
            }
            CigarKind::Deletion => {
                out.push((b'D', rp, FeatureData::Deletion(len)));
                ro += len;
            }
            CigarKind::Skip => {
                out.push((b'N', rp, FeatureData::RefSkip(len)));
                ro += len;
            }
            CigarKind::HardClip => out.push((b'H', rp, FeatureData::HardClip(len))),
            CigarKind::Pad => out.push((b'P', rp, FeatureData::Padding(len))),
        }
    }
    Ok(out)
}

/// Writes one record into the encoder's buffers.
fn encode_record(
    enc: &mut Encoder,
    rec: &Record,
    tag_line: usize,
    window: Option<&[u8]>,
) -> Result<()> {
    let mapped = rec.flags.0 & Flags::UNMAPPED == 0;
    if rec.name.contains(&0) {
        return Err(SeqError::malformed(FORMAT, "read name contains NUL"));
    }
    if !rec.qual.is_empty() && rec.qual.len() != rec.seq.len() {
        return Err(SeqError::malformed(FORMAT, "sequence and quality lengths differ"));
    }
    let flags = i32::from(rec.flags.0);
    enc.int(*b"BF", flags)?;
    enc.int(*b"CF", 0x3)?;
    enc.int(*b"RL", rec.seq.len() as i32)?;
    enc.int(*b"AP", checked_i32(rec.pos + 1, "position")?)?;
    enc.int(*b"RG", -1)?;
    enc.stop_array(*b"RN", &rec.name, 0)?;
    let mut mf = 0;
    if flags & i32::from(Flags::MATE_REVERSE) != 0 {
        mf |= 1;
    }
    if flags & i32::from(Flags::MATE_UNMAPPED) != 0 {
        mf |= 2;
    }
    enc.int(*b"MF", mf)?;
    enc.int(*b"NS", rec.next_ref_id)?;
    enc.int(*b"NP", checked_i32(rec.next_pos + 1, "mate position")?)?;
    enc.int(*b"TS", checked_i32(rec.tlen, "template length")?)?;
    enc.int(*b"TL", tag_line as i32)?;
    for tag in &rec.tags {
        let (ty, value) = tag_value(&tag.value)?;
        let mut body = value;
        let key = [tag.key[0], tag.key[1], ty];
        let mut framed = Vec::with_capacity(body.len());
        framed.append(&mut body);
        enc.tag(key, &framed);
    }

    let qual_bytes: Vec<u8> =
        if rec.qual.is_empty() { vec![0xFF; rec.seq.len()] } else { rec.qual.clone() };
    if mapped {
        if rec.cigar.is_empty() {
            return Err(SeqError::Unsupported("mapped record without a CIGAR".into()));
        }
        if query_length(&rec.cigar) != rec.seq.len() as u64 {
            return Err(SeqError::malformed(FORMAT, "CIGAR does not match sequence length"));
        }
        let feats = features_of(rec, window)?;
        enc.int(*b"FN", feats.len() as i32)?;
        let mut prev = 0usize;
        for (code, pos, data) in &feats {
            enc.byte(*b"FC", *code)?;
            enc.int(*b"FP", (pos - prev) as i32)?;
            prev = *pos;
            match data {
                FeatureData::Substitution(c) => enc.byte(*b"BS", *c)?,
                FeatureData::Bases(b) => enc.array(*b"BB", b)?,
                FeatureData::Insertion(b) => enc.array(*b"IN", b)?,
                FeatureData::SoftClip(b) => enc.array(*b"SC", b)?,
                FeatureData::Deletion(n) => enc.int(*b"DL", *n as i32)?,
                FeatureData::RefSkip(n) => enc.int(*b"RS", *n as i32)?,
                FeatureData::HardClip(n) => enc.int(*b"HC", *n as i32)?,
                FeatureData::Padding(n) => enc.int(*b"PD", *n as i32)?,
            }
        }
        enc.int(*b"MQ", i32::from(rec.mapq))?;
    } else {
        for &b in &rec.seq {
            enc.byte(*b"BA", b)?;
        }
    }
    for &q in &qual_bytes {
        enc.byte(*b"QS", q)?;
    }
    Ok(())
}

/// Position span of a set of records: 1-based start and span on the reference.
fn span_of(records: &[Record]) -> (i32, i32) {
    let mapped: Vec<&Record> = records
        .iter()
        .filter(|r| r.flags.0 & Flags::UNMAPPED == 0 && r.pos >= 0 && !r.cigar.is_empty())
        .collect();
    if mapped.is_empty() {
        let start = records.iter().filter(|r| r.pos >= 0).map(|r| r.pos + 1).min().unwrap_or(0);
        return (start as i32, 0);
    }
    let min = mapped.iter().map(|r| r.pos).min().unwrap_or(0);
    let max_end =
        mapped.iter().map(|r| r.pos + reference_length(&r.cigar) as i64).max().unwrap_or(min);
    ((min + 1) as i32, (max_end - min) as i32)
}

/// Writes CRAM 3.0 to any [`Write`] sink.
#[derive(Debug)]
pub struct CramWriter<W: Write> {
    inner: W,
    header: Header,
    reference: Option<CramReference>,
    records_per_slice: usize,
    slices_per_container: usize,
    pending: Vec<Record>,
    counter: i64,
}

impl<W: Write> CramWriter<W> {
    /// Starts a CRAM file with default layout (10,000 records per slice, one slice
    /// per container). The file definition and header container are written now.
    ///
    /// # Errors
    /// Returns an error if the header cannot be written.
    pub fn new(inner: W, header: Header, reference: Option<CramReference>) -> Result<Self> {
        Self::with_layout(inner, header, reference, 10_000, 1)
    }

    /// As [`CramWriter::new`], with explicit slice and container sizes.
    ///
    /// # Errors
    /// Returns [`SeqError::Unsupported`] for zero sizes, or an I/O error.
    pub fn with_layout(
        mut inner: W,
        header: Header,
        reference: Option<CramReference>,
        records_per_slice: usize,
        slices_per_container: usize,
    ) -> Result<Self> {
        if records_per_slice == 0 || slices_per_container == 0 {
            return Err(SeqError::Unsupported(
                "zero records per slice or slices per container".into(),
            ));
        }
        inner.write_all(&container::write_file_definition(&[0u8; 20]))?;
        let text = header.text().as_bytes();
        let mut payload = Vec::with_capacity(text.len() + 4);
        payload.extend_from_slice(&(text.len() as i32).to_le_bytes());
        payload.extend_from_slice(text);
        let mut blocks = Vec::new();
        block::write_block(Method::Raw, content::FILE_HEADER, 0, &payload, &mut blocks)?;
        write_container(&mut inner, &ContainerFields::bare(1), &blocks)?;
        Ok(Self {
            inner,
            header,
            reference,
            records_per_slice,
            slices_per_container,
            pending: Vec::new(),
            counter: 0,
        })
    }

    /// Adds a record. Records are buffered until their container is full or the
    /// reference id changes.
    ///
    /// # Errors
    /// Returns [`SeqError::Unsupported`] for records CRAM cannot express here, and
    /// [`SeqError::Malformed`] for inconsistent records.
    pub fn write_record(&mut self, rec: &Record) -> Result<()> {
        if self.pending.first().is_some_and(|first| first.ref_id != rec.ref_id) {
            self.flush()?;
        }
        self.pending.push(rec.clone());
        if self.pending.len() >= self.records_per_slice * self.slices_per_container {
            self.flush()?;
        }
        Ok(())
    }

    /// Flushes pending records, writes the end-of-file container and returns the sink.
    ///
    /// # Errors
    /// Returns an error if a write fails.
    pub fn finish(mut self) -> Result<W> {
        self.flush()?;
        self.inner.write_all(&EOF_CONTAINER_V3)?;
        self.inner.flush()?;
        Ok(self.inner)
    }

    /// The reference sequence name for a reference id, if the header has it.
    fn ref_name(&self, ref_id: i32) -> Option<&str> {
        usize::try_from(ref_id)
            .ok()
            .and_then(|i| self.header.references().get(i))
            .map(|r| r.name.as_str())
    }

    /// Reference bases spanned by a mapped record, when a reference is available.
    fn record_window(&self, rec: &Record) -> Result<Option<Vec<u8>>> {
        let Some(reference) = &self.reference else { return Ok(None) };
        if rec.flags.0 & Flags::UNMAPPED != 0 || rec.cigar.is_empty() {
            return Ok(None);
        }
        let Some(name) = self.ref_name(rec.ref_id) else { return Ok(None) };
        if rec.pos < 0 {
            return Err(SeqError::malformed(FORMAT, "mapped record without a position"));
        }
        let beg = rec.pos as u64;
        let end = beg + reference_length(&rec.cigar);
        Ok(Some(reference.fetch(name, beg, end)?))
    }

    /// MD5 of the reference region a slice covers, when a reference is available.
    fn slice_md5(&self, ref_id: i32, start: i32, span: i32) -> Result<Option<[u8; 16]>> {
        let (Some(reference), Some(name)) = (&self.reference, self.ref_name(ref_id)) else {
            return Ok(None);
        };
        if start <= 0 || span <= 0 {
            return Ok(None);
        }
        let beg = (start - 1) as u64;
        Ok(Some(reference.region_md5(name, beg, beg + span as u64)?))
    }

    /// Writes buffered records as one container.
    fn flush(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let records = std::mem::take(&mut self.pending);
        let ref_id = records[0].ref_id;

        // Pre-scan: allocate every block id and the tag dictionary, so the
        // compression header is known before any slice is encoded.
        let mut plan = Plan::new();
        let mut lines = Vec::with_capacity(records.len());
        for rec in &records {
            let mut keys = Vec::with_capacity(rec.tags.len());
            for tag in &rec.tags {
                let (ty, _) = tag_value(&tag.value)?;
                let key = [tag.key[0], tag.key[1], ty];
                plan.tag_key(key);
                keys.push(key);
            }
            lines.push(plan.tag_line(keys));
        }
        let ext_ids = plan.ext_ids.clone();
        let comp = compression::write(&plan.compression_header(self.reference.is_some()));
        let mut enc = Encoder { plan, bufs: HashMap::new() };

        let mut blocks = Vec::new();
        block::write_block(Method::Raw, content::COMPRESSION_HEADER, 0, &comp, &mut blocks)?;
        let mut landmarks = Vec::new();
        let mut n_blocks = 1i32;
        let mut n_bases = 0i64;
        for (chunk_no, chunk) in records.chunks(self.records_per_slice).enumerate() {
            let first = chunk_no * self.records_per_slice;
            for (i, rec) in chunk.iter().enumerate() {
                let window = self.record_window(rec)?;
                encode_record(&mut enc, rec, lines[first + i], window.as_deref())?;
                n_bases += rec.seq.len() as i64;
            }
            let (start, span) = span_of(chunk);
            let sh = SliceHeader {
                ref_id,
                start,
                span,
                n_records: chunk.len() as i32,
                record_counter: self.counter + first as i64,
                n_blocks: 1 + ext_ids.len() as i32,
                content_ids: ext_ids.clone(),
                embedded_ref_id: None,
                md5: self.slice_md5(ref_id, start, span)?,
                tags: Vec::new(),
            };
            landmarks.push(blocks.len() as i32);
            block::write_block(
                Method::Raw,
                content::SLICE_HEADER,
                0,
                &slice::write(&sh),
                &mut blocks,
            )?;
            block::write_block(Method::Raw, content::CORE_DATA, 0, &[], &mut blocks)?;
            for &id in &ext_ids {
                let data = enc.bufs.remove(&id).unwrap_or_default();
                block::write_block(Method::Gzip, content::EXTERNAL_DATA, id, &data, &mut blocks)?;
            }
            enc.bufs = ext_ids.iter().map(|&id| (id, Vec::new())).collect();
            n_blocks += 2 + ext_ids.len() as i32;
        }

        let (start, span) = span_of(&records);
        let fields = ContainerFields {
            ref_id,
            start,
            span,
            n_records: records.len() as i32,
            record_counter: self.counter,
            n_bases,
            landmarks,
            n_blocks,
        };
        write_container(&mut self.inner, &fields, &blocks)?;
        self.counter += records.len() as i64;
        Ok(())
    }
}

/// Fields of a container header that the writer computes.
struct ContainerFields {
    ref_id: i32,
    start: i32,
    span: i32,
    n_records: i32,
    record_counter: i64,
    n_bases: i64,
    landmarks: Vec<i32>,
    n_blocks: i32,
}

impl ContainerFields {
    /// Fields for a container with no records (file header container).
    fn bare(n_blocks: i32) -> Self {
        Self {
            ref_id: -1,
            start: 0,
            span: 0,
            n_records: 0,
            record_counter: 0,
            n_bases: 0,
            landmarks: Vec::new(),
            n_blocks,
        }
    }
}

/// Writes a container header followed by its blocks.
fn write_container<W: Write>(out: &mut W, f: &ContainerFields, blocks: &[u8]) -> Result<()> {
    let header = ContainerHeader {
        length: i32::try_from(blocks.len())
            .map_err(|_| SeqError::Unsupported("container larger than 2 GiB".into()))?,
        ref_id: f.ref_id,
        start: f.start,
        span: f.span,
        n_records: f.n_records,
        record_counter: f.record_counter,
        n_bases: f.n_bases,
        n_blocks: f.n_blocks,
        landmarks: f.landmarks.clone(),
    };
    out.write_all(&container::write_container_header(&header))?;
    out.write_all(blocks)?;
    Ok(())
}
