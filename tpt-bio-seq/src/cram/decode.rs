//! Record decoding for one CRAM slice.
//!
//! Per record the fields are read in the order used by the CRAM 3.0 data
//! series: flags, CRAM flags, reference id (multi-reference slices only), read
//! length, alignment start, read group, read name, detached mate fields or the
//! downstream-mate distance, tag line and tags. Mapped reads then carry read
//! features (FN and the feature list) and a mapping quality; unmapped reads carry
//! their bases. A CRAM flags bit can add a quality array last.
//!
//! Mapped bases are rebuilt by copying from the reference and applying the
//! features. Records that need reference bases fail with
//! [`SeqError::Unsupported`] when no reference was supplied.
//!
//! Substitution codes use the spec's substitution matrix: for reference base `r`,
//! the four alternatives are the bases of `ACGTN` other than `r` (in that order,
//! with `N` replaced by nothing for reference `N`), and the two-bit field `i` of
//! the matrix byte gives the code of alternative `i`.

use std::collections::HashMap;

use crate::bam::decode_tags;
use crate::cram::block::{Block, content};
use crate::cram::codec::{Codec, DataReader};
use crate::cram::compression::CompressionHeader;
use crate::cram::reference::CramReference;
use crate::cram::slice::SliceHeader;
use crate::error::{Result, SeqError};
use crate::header::Reference;
use crate::record::{CigarKind, CigarOp, Flags, Record, Tag};

const FORMAT: &str = "CRAM record";
const NULL: Codec = Codec::Null;
/// The bases that substitution codes index into, for each reference base.
const BASES: &[u8; 5] = b"ACGTN";

/// Data-series and tag codecs for one compression header.
#[derive(Debug)]
pub struct Series {
    map: HashMap<[u8; 2], Codec>,
    tags: HashMap<[u8; 3], Codec>,
    /// Tag dictionary: each tag line is a list of 3-byte keys.
    pub tag_lines: Vec<Vec<[u8; 3]>>,
    /// Substitution matrix (five bytes, one per reference base).
    pub substitution_matrix: [u8; 5],
    /// Whether read names are stored.
    pub read_names: bool,
    /// Whether alignment starts are deltas from the previous record.
    pub ap_delta: bool,
}

impl Series {
    /// Builds codecs for every data series and tag in `h`.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when an encoding has invalid parameters.
    pub fn new(h: &CompressionHeader) -> Result<Self> {
        let mut map = HashMap::new();
        for (key, enc) in &h.data_series {
            map.insert(*key, Codec::from_encoding(enc)?);
        }
        let mut tags = HashMap::new();
        for (key, enc) in &h.tags {
            tags.insert(*key, Codec::from_encoding(enc)?);
        }
        Ok(Self {
            map,
            tags,
            tag_lines: h.preservation.tag_dictionary.clone(),
            substitution_matrix: h.preservation.substitution_matrix,
            read_names: h.preservation.read_names,
            ap_delta: h.preservation.ap_delta,
        })
    }

    /// The codec for a data series. A missing series reads as the NULL encoding.
    #[must_use]
    pub fn series(&self, key: &[u8; 2]) -> &Codec {
        self.map.get(key).unwrap_or(&NULL)
    }
}

/// One read feature, with its 1-based read position.
#[derive(Debug)]
enum Feature {
    Substitution { pos: usize, code: u8 },
    Base { pos: usize, base: u8, qual: u8 },
    Insertion { pos: usize, bases: Vec<u8> },
    SingleBase { pos: usize, base: u8 },
    Deletion { pos: usize, len: usize },
    RefSkip { pos: usize, len: usize },
    SoftClip { pos: usize, bases: Vec<u8> },
    HardClip { pos: usize, len: usize },
    Padding { pos: usize, len: usize },
    Quality { pos: usize, qual: u8 },
    Bases { pos: usize, bases: Vec<u8> },
    Qualities { pos: usize, quals: Vec<u8> },
}

impl Feature {
    fn pos(&self) -> usize {
        match self {
            Self::Substitution { pos, .. }
            | Self::Base { pos, .. }
            | Self::Insertion { pos, .. }
            | Self::SingleBase { pos, .. }
            | Self::Deletion { pos, .. }
            | Self::RefSkip { pos, .. }
            | Self::SoftClip { pos, .. }
            | Self::HardClip { pos, .. }
            | Self::Padding { pos, .. }
            | Self::Quality { pos, .. }
            | Self::Bases { pos, .. }
            | Self::Qualities { pos, .. } => *pos,
        }
    }

    /// Reference bases this feature skips over (deletions and skips).
    fn ref_skip(&self) -> usize {
        match self {
            Self::Deletion { len, .. } | Self::RefSkip { len, .. } => *len,
            _ => 0,
        }
    }
}

/// Substitute the base for code `code` against reference base `r`.
///
/// # Errors
/// Returns [`SeqError::Malformed`] when the matrix has no alternative with that code.
pub fn substitute(matrix: &[u8; 5], r: u8, code: u8) -> Result<u8> {
    let ri = ref_index(r);
    let alts = alternatives(ri);
    for (i, &base) in alts.iter().enumerate() {
        if (matrix[ri] >> (6 - 2 * i)) & 3 == code {
            return Ok(base);
        }
    }
    Err(SeqError::malformed(FORMAT, format!("substitution code {code} for base {}", r as char)))
}

/// The substitution code of `base` against reference base `r`, if representable.
#[must_use]
pub fn substitution_code(matrix: &[u8; 5], r: u8, base: u8) -> Option<u8> {
    let ri = ref_index(r);
    let alts = alternatives(ri);
    let i = alts.iter().position(|&b| b == base)?;
    Some((matrix[ri] >> (6 - 2 * i)) & 3)
}

fn ref_index(r: u8) -> usize {
    BASES.iter().position(|&b| b == r.to_ascii_uppercase()).unwrap_or(4)
}

fn alternatives(ri: usize) -> Vec<u8> {
    BASES.iter().enumerate().filter(|(j, _)| *j != ri).map(|(_, &b)| b).collect()
}

fn usize_of(v: i32, what: &str) -> Result<usize> {
    usize::try_from(v).map_err(|_| SeqError::malformed(FORMAT, format!("negative {what}")))
}

fn byte_of(v: i32, what: &str) -> Result<u8> {
    u8::try_from(v).map_err(|_| SeqError::malformed(FORMAT, format!("{what} out of byte range")))
}

/// Builds the record's bases, CIGAR and per-base quality overrides from features.
struct Rebuild<'a> {
    read_len: usize,
    window: Option<&'a [u8]>,
    ref_off: usize,
    seq: Vec<u8>,
    cigar: Vec<CigarOp>,
    quals: Vec<(usize, u8)>,
    matrix: [u8; 5],
}

impl Rebuild<'_> {
    fn ref_base(&self) -> Result<u8> {
        let w = self.window.ok_or_else(|| {
            SeqError::Unsupported(
                "CRAM record needs bases from a reference; open the file with a reference".into(),
            )
        })?;
        w.get(self.ref_off)
            .copied()
            .ok_or_else(|| SeqError::malformed(FORMAT, "reference region too short for record"))
    }

    fn op(&mut self, kind: CigarKind, len: usize) {
        let len = len as u32;
        if let Some(last) = self.cigar.last_mut()
            && last.kind == kind
        {
            last.len += len;
            return;
        }
        self.cigar.push(CigarOp { len, kind });
    }

    /// Copies reference bases until `upto - 1` bases have been emitted.
    fn copy_to(&mut self, upto: usize) -> Result<()> {
        if upto < self.seq.len() + 1 {
            return Err(SeqError::malformed(FORMAT, "read features out of order"));
        }
        while self.seq.len() + 1 < upto {
            let r = self.ref_base()?;
            self.seq.push(r.to_ascii_uppercase());
            self.ref_off += 1;
            self.op(CigarKind::Match, 1);
        }
        Ok(())
    }

    fn apply(&mut self, f: &Feature) -> Result<()> {
        match f {
            Feature::Substitution { pos, code } => {
                self.copy_to(*pos)?;
                let r = self.ref_base()?;
                self.seq.push(substitute(&self.matrix, r, *code)?);
                self.ref_off += 1;
                self.op(CigarKind::Match, 1);
            }
            Feature::Base { pos, base, qual } => {
                self.copy_to(*pos)?;
                self.seq.push(*base);
                self.quals.push((self.seq.len() - 1, *qual));
                self.ref_off += 1;
                self.op(CigarKind::Match, 1);
            }
            Feature::Insertion { pos, bases } => {
                self.copy_to(*pos)?;
                self.seq.extend_from_slice(bases);
                self.op(CigarKind::Insertion, bases.len());
            }
            Feature::SingleBase { pos, base } => {
                self.copy_to(*pos)?;
                self.seq.push(*base);
                self.op(CigarKind::Insertion, 1);
            }
            Feature::Deletion { pos, len } => {
                self.copy_to(*pos)?;
                self.ref_off += len;
                self.op(CigarKind::Deletion, *len);
            }
            Feature::RefSkip { pos, len } => {
                self.copy_to(*pos)?;
                self.ref_off += len;
                self.op(CigarKind::Skip, *len);
            }
            Feature::SoftClip { pos, bases } => {
                self.copy_to(*pos)?;
                self.seq.extend_from_slice(bases);
                self.op(CigarKind::SoftClip, bases.len());
            }
            Feature::HardClip { pos, len } => {
                self.copy_to(*pos)?;
                self.op(CigarKind::HardClip, *len);
            }
            Feature::Padding { pos, len } => {
                self.copy_to(*pos)?;
                self.op(CigarKind::Pad, *len);
            }
            Feature::Quality { pos, qual } => {
                self.quals.push((pos - 1, *qual));
            }
            Feature::Bases { pos, bases } => {
                self.copy_to(*pos)?;
                self.seq.extend_from_slice(bases);
                self.ref_off += bases.len();
                self.op(CigarKind::Match, bases.len());
            }
            Feature::Qualities { pos, quals } => {
                for (i, &q) in quals.iter().enumerate() {
                    self.quals.push((pos - 1 + i, q));
                }
            }
        }
        Ok(())
    }
}

/// State while decoding the records of one slice.
#[derive(Debug)]
pub struct SliceDecoder<'a> {
    series: &'a Series,
    reader: DataReader,
    refs: &'a [Reference],
    reference: Option<&'a CramReference>,
    slice_ref_id: i32,
    prev_ap: i64,
}

impl SliceDecoder<'_> {
    fn read_int(&mut self, key: [u8; 2]) -> Result<i32> {
        let series = self.series;
        self.reader.int(series.series(&key))
    }

    fn read_byte(&mut self, key: [u8; 2]) -> Result<u8> {
        let series = self.series;
        self.reader.byte(series.series(&key))
    }

    fn read_array(&mut self, key: [u8; 2]) -> Result<Vec<u8>> {
        let series = self.series;
        self.reader.byte_array(series.series(&key))
    }

    fn read_feature(&mut self, code: u8, pos: usize) -> Result<Feature> {
        Ok(match code {
            b'X' => Feature::Substitution { pos, code: self.read_byte(*b"BS")? },
            b'B' => {
                Feature::Base { pos, base: self.read_byte(*b"BA")?, qual: self.read_byte(*b"QS")? }
            }
            b'I' => Feature::Insertion { pos, bases: self.read_array(*b"IN")? },
            b'i' => Feature::SingleBase { pos, base: self.read_byte(*b"BA")? },
            b'D' => Feature::Deletion { pos, len: usize_of(self.read_int(*b"DL")?, "deletion")? },
            b'N' => Feature::RefSkip { pos, len: usize_of(self.read_int(*b"RS")?, "skip")? },
            b'S' => Feature::SoftClip { pos, bases: self.read_array(*b"SC")? },
            b'H' => Feature::HardClip { pos, len: usize_of(self.read_int(*b"HC")?, "clip")? },
            b'P' => Feature::Padding { pos, len: usize_of(self.read_int(*b"PD")?, "padding")? },
            b'Q' => Feature::Quality { pos, qual: self.read_byte(*b"QS")? },
            b'b' => Feature::Bases { pos, bases: self.read_array(*b"BB")? },
            b'q' => Feature::Qualities { pos, quals: self.read_array(*b"QQ")? },
            other => {
                return Err(SeqError::malformed(
                    FORMAT,
                    format!("unknown read feature code {}", other as char),
                ));
            }
        })
    }

    /// Decodes one record. Returns the record and its downstream-mate distance (NF).
    fn record(&mut self) -> Result<(Record, Option<usize>)> {
        let series = self.series;
        let bf = self.read_int(*b"BF")?;
        let cf = self.read_int(*b"CF")?;
        let ref_id =
            if self.slice_ref_id == -2 { self.read_int(*b"RI")? } else { self.slice_ref_id };
        let read_len = usize_of(self.read_int(*b"RL")?, "read length")?;
        let ap = i64::from(self.read_int(*b"AP")?);
        let ap = if series.ap_delta {
            let v = self.prev_ap + ap;
            self.prev_ap = v;
            v
        } else {
            ap
        };
        let pos = ap - 1;
        let _read_group = self.read_int(*b"RG")?;
        let name = if series.read_names { self.read_array(*b"RN")? } else { Vec::new() };

        let mut flags = u16::try_from(bf)
            .map_err(|_| SeqError::malformed(FORMAT, format!("flags {bf} out of range")))?;
        let mut next_ref_id = -1;
        let mut next_pos = -1;
        let mut tlen = 0;
        if cf & 0x2 != 0 {
            // Detached mate: MF holds mate-reverse (bit 0) and mate-unmapped (bit 1).
            let mf = self.read_int(*b"MF")?;
            if mf & 1 != 0 {
                flags |= Flags::MATE_REVERSE;
            }
            if mf & 2 != 0 {
                flags |= Flags::MATE_UNMAPPED;
            }
            next_ref_id = self.read_int(*b"NS")?;
            next_pos = i64::from(self.read_int(*b"NP")?) - 1;
            tlen = i64::from(self.read_int(*b"TS")?);
        }
        let downstream = if cf & 0x4 != 0 {
            Some(usize_of(self.read_int(*b"NF")?, "mate distance")?)
        } else {
            None
        };

        let tag_line = self.read_int(*b"TL")?;
        let keys = usize::try_from(tag_line)
            .ok()
            .and_then(|i| series.tag_lines.get(i))
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("tag line {tag_line} not found")))?;
        let mut tags: Vec<Tag> = Vec::new();
        for key in keys {
            let codec = series.tags.get(key).unwrap_or(&NULL);
            let value = self.reader.byte_array(codec)?;
            let mut buf = Vec::with_capacity(3 + value.len());
            buf.extend_from_slice(&key[..3]);
            buf.extend_from_slice(&value);
            tags.extend(decode_tags(&buf)?);
        }

        let mapped = flags & Flags::UNMAPPED == 0;
        let mut mapq = 0u8;
        let (seq, cigar, mut quals) = if mapped {
            let n_features = usize_of(self.read_int(*b"FN")?, "feature count")?;
            let mut features = Vec::with_capacity(n_features.min(1024));
            let mut prev = 0usize;
            for _ in 0..n_features {
                let code = self.read_byte(*b"FC")?;
                let delta = usize_of(self.read_int(*b"FP")?, "feature position")?;
                let fp = prev
                    .checked_add(delta)
                    .ok_or_else(|| SeqError::malformed(FORMAT, "feature position overflow"))?;
                prev = fp;
                features.push(self.read_feature(code, fp)?);
            }
            mapq = byte_of(self.read_int(*b"MQ")?, "mapping quality")?;
            self.rebuild(ref_id, pos, read_len, &features)?
        } else {
            let mut seq = Vec::with_capacity(read_len.min(1 << 16));
            for _ in 0..read_len {
                seq.push(self.read_byte(*b"BA")?);
            }
            (seq, Vec::new(), vec![0xFF; read_len])
        };
        if cf & 0x1 != 0 {
            // Quality array: one byte per base, overriding feature qualities.
            quals.clear();
            for _ in 0..read_len {
                quals.push(self.read_byte(*b"QS")?);
            }
        }
        // 0xFF in every position is the BAM convention for "no qualities".
        let qual = if quals.iter().all(|&q| q == 0xFF) { Vec::new() } else { quals };

        let rec = Record {
            name,
            flags: Flags(flags),
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
        };
        Ok((rec, downstream))
    }

    /// Rebuilds bases, CIGAR and qualities for a mapped read.
    fn rebuild(
        &self,
        ref_id: i32,
        pos: i64,
        read_len: usize,
        features: &[Feature],
    ) -> Result<(Vec<u8>, Vec<CigarOp>, Vec<u8>)> {
        let window = match self.reference {
            Some(reference) if ref_id >= 0 && pos >= 0 => {
                let name = self
                    .refs
                    .get(ref_id as usize)
                    .map(|r| r.name.as_str())
                    .ok_or_else(|| SeqError::malformed(FORMAT, "reference id not in header"))?;
                let skip: usize = features.iter().map(Feature::ref_skip).sum();
                let beg = pos as u64;
                let end = beg + (read_len + skip) as u64 + 1;
                Some(reference.fetch(name, beg, end)?)
            }
            _ => None,
        };
        let mut b = Rebuild {
            read_len,
            window: window.as_deref(),
            ref_off: 0,
            seq: Vec::with_capacity(read_len),
            cigar: Vec::new(),
            quals: Vec::new(),
            matrix: self.series.substitution_matrix,
        };
        for f in features {
            if f.pos() < 1 {
                return Err(SeqError::malformed(FORMAT, "feature position below 1"));
            }
            b.apply(f)?;
        }
        b.copy_to(b.read_len + 1)?;
        if b.seq.len() != read_len {
            return Err(SeqError::malformed(FORMAT, "read features do not match read length"));
        }
        let mut quals = vec![0xFFu8; read_len];
        for (i, q) in b.quals {
            if let Some(slot) = quals.get_mut(i) {
                *slot = q;
            }
        }
        Ok((b.seq, b.cigar, quals))
    }
}

/// Decodes every record of a slice.
///
/// `blocks` are the slice's data blocks (core and external), in file order.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for inconsistent data or a reference MD5
/// mismatch, and [`SeqError::Unsupported`] for records that need a reference
/// when none was supplied, or for embedded reference blocks.
pub fn decode_slice(
    series: &Series,
    sh: &SliceHeader,
    blocks: Vec<Block>,
    refs: &[Reference],
    reference: Option<&CramReference>,
) -> Result<Vec<Record>> {
    if sh.embedded_ref_id.is_some() {
        return Err(SeqError::Unsupported("CRAM slices with embedded reference blocks".into()));
    }
    let n_records = usize_of(sh.n_records, "record count")?;
    let mut core = Vec::new();
    let mut externals = Vec::new();
    for b in blocks {
        match b.content_type {
            content::CORE_DATA => core = b.data,
            content::EXTERNAL_DATA => externals.push((b.content_id, b.data)),
            _ => {}
        }
    }
    if let (Some(expected), Some(reference)) = (sh.md5, reference)
        && sh.ref_id >= 0
        && sh.start > 0
        && sh.span > 0
        && let Some(r) = refs.get(sh.ref_id as usize)
    {
        let beg = (sh.start - 1) as u64;
        let end = beg + sh.span as u64;
        if reference.region_md5(&r.name, beg, end)? != expected {
            return Err(SeqError::malformed(FORMAT, "reference MD5 does not match slice header"));
        }
    }

    let mut dec = SliceDecoder {
        series,
        reader: DataReader::new(core, externals),
        refs,
        reference,
        slice_ref_id: sh.ref_id,
        prev_ap: i64::from(sh.start),
    };
    let mut records = Vec::with_capacity(n_records.min(1 << 16));
    let mut downstream = Vec::with_capacity(n_records.min(1 << 16));
    for _ in 0..n_records {
        let (rec, nf) = dec.record()?;
        records.push(rec);
        downstream.push(nf);
    }
    // A downstream mate is found by its record distance. Only the upstream record
    // gets the link; detached records already carry their mate fields.
    for (i, nf) in downstream.iter().enumerate() {
        if let Some(d) = *nf {
            let j = i + d + 1;
            if j < records.len() {
                let (ref_id, pos) = (records[j].ref_id, records[j].pos);
                records[i].next_ref_id = ref_id;
                records[i].next_pos = pos;
            }
        }
    }
    Ok(records)
}
