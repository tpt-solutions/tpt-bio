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

use crate::bam::decode_tags;
use crate::cram::block::{Block, content};
use crate::cram::codec::{Codec, DataReader, FastMap};
use crate::cram::compression::CompressionHeader;
use crate::cram::reference::CramReference;
use crate::cram::slice::SliceHeader;
use crate::error::{Result, SeqError};
use crate::header::Reference;
use crate::record::{CigarKind, CigarOp, Flags, Record, Tag, TagValue};
use std::ops::Range;

/// Longest read length accepted from a CRAM file. Real reads are far shorter.
const MAX_READ_LEN: usize = 1 << 24;

const FORMAT: &str = "CRAM record";
const NULL: Codec = Codec::Null;
/// The bases that substitution codes index into, for each reference base.
const BASES: &[u8; 5] = b"ACGTN";

/// Data-series codecs resolved once per compression header, so the record loop
/// reads each series without a map lookup.
#[derive(Debug)]
struct Fields {
    bam_flags: Codec,
    cram_flags: Codec,
    ref_id: Codec,
    read_len: Codec,
    align_pos: Codec,
    read_group: Codec,
    name: Codec,
    mate_flags: Codec,
    mate_ref: Codec,
    mate_pos: Codec,
    tlen: Codec,
    mate_distance: Codec,
    tag_line: Codec,
    feature_count: Codec,
    feature_code: Codec,
    feature_pos: Codec,
    mapq: Codec,
    base: Codec,
    qual: Codec,
    subst: Codec,
    insertion: Codec,
    deletion_len: Codec,
    skip_len: Codec,
    soft_clip: Codec,
    hard_clip_len: Codec,
    padding_len: Codec,
    bases: Codec,
    quals: Codec,
}

/// Data-series and tag codecs for one compression header.
#[derive(Debug)]
pub struct Series {
    map: FastMap<[u8; 2], Codec>,
    fields: Fields,
    tags: FastMap<[u8; 3], Codec>,
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
        let mut map = FastMap::default();
        for (key, enc) in &h.data_series {
            map.insert(*key, Codec::from_encoding(enc)?);
        }
        let mut tags = FastMap::default();
        for (key, enc) in &h.tags {
            tags.insert(*key, Codec::from_encoding(enc)?);
        }
        let get = |key: &[u8; 2]| map.get(key).cloned().unwrap_or(Codec::Null);
        let fields = Fields {
            bam_flags: get(b"BF"),
            cram_flags: get(b"CF"),
            ref_id: get(b"RI"),
            read_len: get(b"RL"),
            align_pos: get(b"AP"),
            read_group: get(b"RG"),
            name: get(b"RN"),
            mate_flags: get(b"MF"),
            mate_ref: get(b"NS"),
            mate_pos: get(b"NP"),
            tlen: get(b"TS"),
            mate_distance: get(b"NF"),
            tag_line: get(b"TL"),
            feature_count: get(b"FN"),
            feature_code: get(b"FC"),
            feature_pos: get(b"FP"),
            mapq: get(b"MQ"),
            base: get(b"BA"),
            qual: get(b"QS"),
            subst: get(b"BS"),
            insertion: get(b"IN"),
            deletion_len: get(b"DL"),
            skip_len: get(b"RS"),
            soft_clip: get(b"SC"),
            hard_clip_len: get(b"HC"),
            padding_len: get(b"PD"),
            bases: get(b"BB"),
            quals: get(b"QQ"),
        };
        Ok(Self {
            map,
            fields,
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
    Substitution {
        pos: usize,
        code: u8,
    },
    Base {
        pos: usize,
        base: u8,
        qual: u8,
    },
    /// Byte-array features hold their bytes as a range of the slice's scratch buffer.
    Insertion {
        pos: usize,
        bases: Range<usize>,
    },
    SingleBase {
        pos: usize,
        base: u8,
    },
    Deletion {
        pos: usize,
        len: usize,
    },
    RefSkip {
        pos: usize,
        len: usize,
    },
    SoftClip {
        pos: usize,
        bases: Range<usize>,
    },
    HardClip {
        pos: usize,
        len: usize,
    },
    Padding {
        pos: usize,
        len: usize,
    },
    Quality {
        pos: usize,
        qual: u8,
    },
    Bases {
        pos: usize,
        bases: Range<usize>,
    },
    Qualities {
        pos: usize,
        quals: Range<usize>,
    },
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
    /// The bytes that byte-array features point into.
    scratch: &'a [u8],
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
                let src = &self.scratch[bases.clone()];
                self.seq.extend_from_slice(src);
                self.op(CigarKind::Insertion, src.len());
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
                let src = &self.scratch[bases.clone()];
                self.seq.extend_from_slice(src);
                self.op(CigarKind::SoftClip, src.len());
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
                let src = &self.scratch[bases.clone()];
                self.seq.extend_from_slice(src);
                self.ref_off += src.len();
                self.op(CigarKind::Match, src.len());
            }
            Feature::Qualities { pos, quals } => {
                for (i, &q) in self.scratch[quals.clone()].iter().enumerate() {
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
    /// Reused per record: the read's features and the bytes they point into.
    features: Vec<Feature>,
    scratch: Vec<u8>,
    refs: &'a [Reference],
    read_groups: &'a [String],
    reference: Option<&'a CramReference>,
    slice_ref_id: i32,
    prev_ap: i64,
}

impl SliceDecoder<'_> {
    /// Decodes the tag values for one tag line, in dictionary order.
    fn read_tag_values(&mut self, keys: &[[u8; 3]]) -> Result<Vec<Tag>> {
        let series = self.series;
        let mut tags: Vec<Tag> = Vec::new();
        for key in keys {
            let codec = series.tags.get(key).unwrap_or(&NULL);
            let value = self.reader.byte_array(codec)?;
            let mut buf = Vec::with_capacity(3 + value.len());
            buf.extend_from_slice(key);
            buf.extend_from_slice(&value);
            tags.extend(decode_tags(&buf)?);
        }
        Ok(tags)
    }

    fn read_int(&mut self, codec: &Codec) -> Result<i32> {
        self.reader.int(codec)
    }

    fn read_byte(&mut self, codec: &Codec) -> Result<u8> {
        self.reader.byte(codec)
    }

    fn read_array(&mut self, codec: &Codec) -> Result<Vec<u8>> {
        self.reader.byte_array(codec)
    }

    /// Appends one byte array to `scratch` and returns its range there.
    fn read_array_into(&mut self, codec: &Codec, scratch: &mut Vec<u8>) -> Result<Range<usize>> {
        let start = scratch.len();
        self.reader.byte_array_into(codec, scratch)?;
        Ok(start..scratch.len())
    }

    /// Reads one read feature. Byte-array payloads are appended to `scratch`.
    fn read_feature(&mut self, code: u8, pos: usize, scratch: &mut Vec<u8>) -> Result<Feature> {
        let series = self.series;
        Ok(match code {
            b'X' => Feature::Substitution { pos, code: self.read_byte(&series.fields.subst)? },
            b'B' => Feature::Base {
                pos,
                base: self.read_byte(&series.fields.base)?,
                qual: self.read_byte(&series.fields.qual)?,
            },
            b'I' => Feature::Insertion {
                pos,
                bases: self.read_array_into(&series.fields.insertion, scratch)?,
            },
            b'i' => Feature::SingleBase { pos, base: self.read_byte(&series.fields.base)? },
            b'D' => Feature::Deletion {
                pos,
                len: usize_of(self.read_int(&series.fields.deletion_len)?, "deletion")?,
            },
            b'N' => Feature::RefSkip {
                pos,
                len: usize_of(self.read_int(&series.fields.skip_len)?, "skip")?,
            },
            b'S' => Feature::SoftClip {
                pos,
                bases: self.read_array_into(&series.fields.soft_clip, scratch)?,
            },
            b'H' => Feature::HardClip {
                pos,
                len: usize_of(self.read_int(&series.fields.hard_clip_len)?, "clip")?,
            },
            b'P' => Feature::Padding {
                pos,
                len: usize_of(self.read_int(&series.fields.padding_len)?, "padding")?,
            },
            b'Q' => Feature::Quality { pos, qual: self.read_byte(&series.fields.qual)? },
            b'b' => {
                Feature::Bases { pos, bases: self.read_array_into(&series.fields.bases, scratch)? }
            }
            b'q' => Feature::Qualities {
                pos,
                quals: self.read_array_into(&series.fields.quals, scratch)?,
            },
            other => {
                return Err(SeqError::malformed(
                    FORMAT,
                    format!("unknown read feature code {}", other as char),
                ));
            }
        })
    }

    /// Reads the features and mapping quality of a mapped read, and rebuilds its bases.
    fn read_mapped(
        &mut self,
        ref_id: i32,
        pos: i64,
        read_len: usize,
    ) -> Result<(u8, Vec<u8>, Vec<CigarOp>, Vec<u8>)> {
        let series = self.series;
        let n_features = usize_of(self.read_int(&series.fields.feature_count)?, "feature count")?;
        // The feature list and its byte payloads reuse the slice's scratch buffers.
        let mut features = std::mem::take(&mut self.features);
        let mut scratch = std::mem::take(&mut self.scratch);
        features.clear();
        scratch.clear();
        let mut prev = 0usize;
        for _ in 0..n_features {
            let code = self.read_byte(&series.fields.feature_code)?;
            let delta = usize_of(self.read_int(&series.fields.feature_pos)?, "feature position")?;
            let fp = prev
                .checked_add(delta)
                .ok_or_else(|| SeqError::malformed(FORMAT, "feature position overflow"))?;
            prev = fp;
            features.push(self.read_feature(code, fp, &mut scratch)?);
        }
        let mapq = byte_of(self.read_int(&series.fields.mapq)?, "mapping quality")?;
        let rebuilt = self.rebuild(ref_id, pos, read_len, &features, &scratch);
        self.features = features;
        self.scratch = scratch;
        let (seq, cigar, quals) = rebuilt?;
        Ok((mapq, seq, cigar, quals))
    }

    /// Decodes one record. Returns the record and its downstream-mate distance (NF).
    fn record(&mut self) -> Result<(Record, Option<usize>)> {
        let series = self.series;
        let bf = self.read_int(&series.fields.bam_flags)?;
        let cf = self.read_int(&series.fields.cram_flags)?;
        let ref_id = if self.slice_ref_id == -2 {
            self.read_int(&series.fields.ref_id)?
        } else {
            self.slice_ref_id
        };
        let read_len = usize_of(self.read_int(&series.fields.read_len)?, "read length")?;
        let ap = i64::from(self.read_int(&series.fields.align_pos)?);
        let ap = if series.ap_delta {
            let v = self.prev_ap + ap;
            self.prev_ap = v;
            v
        } else {
            ap
        };
        let pos = ap - 1;
        let read_group = self.read_int(&series.fields.read_group)?;
        let name =
            if series.read_names { self.read_array(&series.fields.name)? } else { Vec::new() };

        let mut flags = u16::try_from(bf)
            .map_err(|_| SeqError::malformed(FORMAT, format!("flags {bf} out of range")))?;
        let mut next_ref_id = -1;
        let mut next_pos = -1;
        let mut tlen = 0;
        if cf & 0x2 != 0 {
            // Detached mate: MF holds mate-reverse (bit 0) and mate-unmapped (bit 1).
            let mf = self.read_int(&series.fields.mate_flags)?;
            if mf & 1 != 0 {
                flags |= Flags::MATE_REVERSE;
            }
            if mf & 2 != 0 {
                flags |= Flags::MATE_UNMAPPED;
            }
            next_ref_id = self.read_int(&series.fields.mate_ref)?;
            next_pos = i64::from(self.read_int(&series.fields.mate_pos)?) - 1;
            tlen = i64::from(self.read_int(&series.fields.tlen)?);
        }
        let downstream = if cf & 0x4 != 0 {
            Some(usize_of(self.read_int(&series.fields.mate_distance)?, "mate distance")?)
        } else {
            None
        };

        let tag_line = self.read_int(&series.fields.tag_line)?;
        let keys = usize::try_from(tag_line)
            .ok()
            .and_then(|i| series.tag_lines.get(i))
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("tag line {tag_line} not found")))?;
        let mut tags = self.read_tag_values(keys)?;
        // The RG series holds an index into the header's @RG lines. samtools writes it
        // as the RG tag, after the tags from the tag list.
        if let Some(id) = usize::try_from(read_group).ok().and_then(|i| self.read_groups.get(i)) {
            tags.push(Tag { key: *b"RG", value: TagValue::String(id.as_bytes().to_vec()) });
        }

        let mapped = flags & Flags::UNMAPPED == 0;
        let (mapq, seq, cigar, mut quals) = if mapped {
            self.read_mapped(ref_id, pos, read_len)?
        } else {
            // Bases are read one per value, with the codec resolved once.
            let seq = self.reader.byte_run(&series.fields.base, read_len)?;
            (0u8, seq, Vec::new(), vec![0xFF; read_len])
        };
        if cf & 0x1 != 0 {
            // Quality array: one byte per base, overriding feature qualities.
            quals = self.reader.byte_run(&series.fields.qual, read_len)?;
        }
        // 0xFF in every position is the BAM convention for "no qualities".
        let qual = if quals.iter().all(|&q| q == 0xFF) { Vec::new() } else { quals };

        if mapped {
            self.add_md_nm(&mut tags, ref_id, pos, &seq, &cigar)?;
        }

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

    /// Adds NM and MD when the file leaves them out, as samtools does on decode.
    /// Both are derived from the reference, so this needs one.
    fn add_md_nm(
        &self,
        tags: &mut Vec<Tag>,
        ref_id: i32,
        pos: i64,
        seq: &[u8],
        cigar: &[CigarOp],
    ) -> Result<()> {
        let has = |tags: &[Tag], key: &[u8; 2]| tags.iter().any(|t| &t.key == key);
        if has(tags, b"NM") && has(tags, b"MD") {
            return Ok(());
        }
        let (Some(reference), true) = (self.reference, ref_id >= 0 && pos >= 0) else {
            return Ok(());
        };
        let name = self
            .refs
            .get(ref_id as usize)
            .map(|r| r.name.as_str())
            .ok_or_else(|| SeqError::malformed(FORMAT, "reference id not in header"))?;
        let span = crate::record::reference_length(cigar);
        let ref_bases = reference.fetch(name, pos as u64, pos as u64 + span)?;
        let (nm, md) = md_and_nm(&ref_bases, seq, cigar);
        if !has(tags, b"MD") {
            tags.push(Tag { key: *b"MD", value: TagValue::String(md) });
        }
        if !has(tags, b"NM") {
            tags.push(Tag { key: *b"NM", value: TagValue::Int(nm) });
        }
        Ok(())
    }

    /// Rebuilds bases, CIGAR and qualities for a mapped read.
    fn rebuild(
        &self,
        ref_id: i32,
        pos: i64,
        read_len: usize,
        features: &[Feature],
        scratch: &[u8],
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
        // The length comes from the file. Refuse absurd reads before allocating for them.
        if read_len > MAX_READ_LEN {
            return Err(SeqError::malformed(
                FORMAT,
                format!("read length {read_len} is too large"),
            ));
        }
        let mut b = Rebuild {
            read_len,
            window: window.as_deref(),
            scratch,
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
        // No overrides leaves no qualities: the caller turns all-0xFF into "absent" anyway.
        if b.quals.is_empty() {
            return Ok((b.seq, b.cigar, Vec::new()));
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
    read_groups: &[String],
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
        features: Vec::new(),
        scratch: Vec::new(),
        refs,
        read_groups,
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

/// Computes the NM and MD tags of an alignment from its reference bases.
///
/// NM counts mismatches and inserted and deleted bases. MD is the run of matching
/// bases, then the reference base at each mismatch, with deleted bases written as
/// `^` followed by the reference bases. Runs of zero are written, as the spec requires.
fn md_and_nm(ref_bases: &[u8], seq: &[u8], cigar: &[CigarOp]) -> (i64, Vec<u8>) {
    let mut nm = 0i64;
    let mut md = Vec::new();
    let mut run = 0u32;
    let (mut rpos, mut qpos) = (0usize, 0usize);
    for op in cigar {
        let n = op.len as usize;
        match op.kind {
            CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch => {
                for _ in 0..n {
                    let r = ref_bases.get(rpos).map_or(b'N', u8::to_ascii_uppercase);
                    let q = seq.get(qpos).map_or(b'N', u8::to_ascii_uppercase);
                    if r == q {
                        run += 1;
                    } else {
                        md.extend_from_slice(run.to_string().as_bytes());
                        md.push(r);
                        run = 0;
                        nm += 1;
                    }
                    rpos += 1;
                    qpos += 1;
                }
            }
            CigarKind::Insertion => {
                nm += i64::from(op.len);
                qpos += n;
            }
            CigarKind::Deletion => {
                nm += i64::from(op.len);
                md.extend_from_slice(run.to_string().as_bytes());
                md.push(b'^');
                for i in 0..n {
                    md.push(ref_bases.get(rpos + i).map_or(b'N', u8::to_ascii_uppercase));
                }
                run = 0;
                rpos += n;
            }
            CigarKind::Skip => rpos += n,
            CigarKind::SoftClip => qpos += n,
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    md.extend_from_slice(run.to_string().as_bytes());
    (nm, md)
}
