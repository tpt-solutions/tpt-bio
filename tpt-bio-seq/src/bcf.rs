//! BCF2 binary VCF: header, typed values, record codec, and BGZF-wrapped reader and writer.
//!
//! The layout follows the BCF 2.2 specification; BCF 2.1 is accepted because its
//! records are identical. Real BCF files are BGZF-compressed, so use
//! [`BcfReader::new_bgzf`] and [`BcfWriter::new_bgzf`]. Records decode to VCF text, so
//! [`BcfReader::next_record`] returns the same [`VcfRecord`] views as
//! [`crate::vcf::VcfReader`].
//!
//! Encoding takes each INFO and FORMAT value's BCF type from the header's `Type`
//! attribute and uses the smallest integer width that fits the values. `GT` uses the
//! BCF allele encoding. Region queries use [`fetch_lines`] together with a
//! [`TabixIndex`] built by [`build_bcf_index`].

use std::collections::HashMap;
use std::fmt::{Display, Write as _};
use std::io::{Read, Seek, Write};

use crate::bgzf::{BgzfReader, BgzfWriter, read_bounded};
use crate::error::{Result, SeqError};
use crate::tabix::{IndexBuilder, IndexKind, TabixIndex};
use crate::vcf::{VcfHeader, VcfRecord};

const FORMAT: &str = "BCF";
const MAGIC: &[u8; 3] = b"BCF";
const MAX_HEADER_TEXT: usize = 1 << 28;
const MAX_RECORD_BYTES: usize = 1 << 31;

const TYPE_MISSING: u8 = 0;
const TYPE_INT8: u8 = 1;
const TYPE_INT16: u8 = 2;
const TYPE_INT32: u8 = 3;
const TYPE_FLOAT: u8 = 5;
const TYPE_CHAR: u8 = 7;
/// A descriptor count of 15 means the real count follows as a typed integer.
const COUNT_ESCAPE: u8 = 15;
const FLOAT_MISSING: u32 = 0x7F80_0001;
const FLOAT_EOV: u32 = 0x7F80_0002;

/// A BCF header: the VCF header plus the string and contig dictionaries that BCF
/// uses in place of names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BcfHeader {
    vcf: VcfHeader,
    strings: Vec<Option<String>>,
    string_ids: HashMap<String, usize>,
    contigs: Vec<Option<String>>,
    contig_ids: HashMap<String, usize>,
    info_types: HashMap<String, String>,
    format_types: HashMap<String, String>,
}

impl BcfHeader {
    /// Parses VCF header text and builds the BCF dictionaries from it.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for an invalid VCF header or dictionary entry.
    pub fn parse(text: &str) -> Result<Self> {
        Self::from_vcf(VcfHeader::parse(text)?)
    }

    /// Builds the dictionaries from a parsed VCF header.
    ///
    /// `PASS` is always dictionary entry 0. Other `FILTER`, `INFO` and `FORMAT` IDs
    /// take their `IDX` value when given, and otherwise the next free slot in header
    /// order.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when an entry lacks `ID`, has a bad `IDX`, or two
    /// entries claim the same index.
    pub fn from_vcf(vcf: VcfHeader) -> Result<Self> {
        let mut string_entries = Vec::new();
        let mut contig_entries = Vec::new();
        let mut info_types = HashMap::new();
        let mut format_types = HashMap::new();
        for meta in vcf.meta() {
            let kind = meta.key.as_str();
            if !matches!(kind, "FILTER" | "INFO" | "FORMAT" | "contig") {
                continue;
            }
            let id = meta
                .get("ID")
                .ok_or_else(|| SeqError::malformed(FORMAT, format!("##{kind} line without ID")))?;
            let idx = meta
                .get("IDX")
                .map(|v| {
                    v.parse::<usize>()
                        .map_err(|_| SeqError::malformed(FORMAT, format!("bad IDX for {id}: {v}")))
                })
                .transpose()?;
            if kind == "contig" {
                contig_entries.push((id.to_string(), idx));
                continue;
            }
            string_entries.push((id.to_string(), idx));
            let ty = meta.get("Type").unwrap_or_default().to_string();
            match kind {
                "INFO" => {
                    info_types.insert(id.to_string(), ty);
                }
                "FORMAT" => {
                    format_types.insert(id.to_string(), ty);
                }
                _ => {}
            }
        }
        let (strings, string_ids) = build_dictionary(&string_entries, true)?;
        let (contigs, contig_ids) = build_dictionary(&contig_entries, false)?;
        Ok(Self { vcf, strings, string_ids, contigs, contig_ids, info_types, format_types })
    }

    /// The underlying VCF header, including the verbatim header text.
    #[must_use]
    pub fn vcf(&self) -> &VcfHeader {
        &self.vcf
    }

    /// Number of samples declared in the header.
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.vcf.samples().len()
    }

    /// Number of contigs in the contig dictionary (highest index plus one).
    #[must_use]
    pub fn n_contigs(&self) -> usize {
        self.contigs.len()
    }

    /// The dictionary index of a contig name, if declared.
    #[must_use]
    pub fn contig_id(&self, name: &str) -> Option<usize> {
        self.contig_ids.get(name).copied()
    }

    /// The contig name at a dictionary index, if defined.
    #[must_use]
    pub fn contig_name(&self, id: usize) -> Option<&str> {
        self.contigs.get(id).and_then(Option::as_deref)
    }

    fn string_name(&self, id: i64) -> Result<&str> {
        usize::try_from(id)
            .ok()
            .and_then(|i| self.strings.get(i))
            .and_then(Option::as_deref)
            .ok_or_else(|| {
                SeqError::malformed(FORMAT, format!("dictionary index {id} is undefined"))
            })
    }

    fn string_id(&self, name: &str) -> Result<i64> {
        self.string_ids
            .get(name)
            .map(|&i| i as i64)
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("{name} is not declared in header")))
    }

    fn info_type(&self, key: &str) -> Result<&str> {
        self.info_types
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("INFO {key} is not declared")))
    }

    fn format_type(&self, key: &str) -> Result<&str> {
        self.format_types
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("FORMAT {key} is not declared")))
    }
}

/// Assigns dictionary slots. Explicit `IDX` values are placed first; other names fill
/// the lowest free slot in header order. With `pass_first`, slot 0 holds `PASS`.
fn build_dictionary(
    entries: &[(String, Option<usize>)],
    pass_first: bool,
) -> Result<(Vec<Option<String>>, HashMap<String, usize>)> {
    let mut slots: Vec<Option<String>> = Vec::new();
    let mut ids: HashMap<String, usize> = HashMap::new();
    if pass_first {
        place(&mut slots, &mut ids, "PASS", 0)?;
    }
    for (name, idx) in entries {
        if let Some(i) = *idx {
            place(&mut slots, &mut ids, name, i)?;
        }
    }
    for (name, idx) in entries {
        if idx.is_none() && !ids.contains_key(name) {
            let free = (0..=slots.len())
                .find(|&i| slots.get(i).is_none_or(Option::is_none))
                .unwrap_or(slots.len());
            place(&mut slots, &mut ids, name, free)?;
        }
    }
    Ok((slots, ids))
}

fn place(
    slots: &mut Vec<Option<String>>,
    ids: &mut HashMap<String, usize>,
    name: &str,
    index: usize,
) -> Result<()> {
    if let Some(&prev) = ids.get(name) {
        if prev == index {
            return Ok(());
        }
        return Err(SeqError::malformed(FORMAT, format!("{name} has conflicting IDX values")));
    }
    if slots.len() <= index {
        slots.resize(index + 1, None);
    }
    if slots[index].is_some() {
        return Err(SeqError::malformed(FORMAT, format!("IDX {index} is used twice")));
    }
    slots[index] = Some(name.to_string());
    ids.insert(name.to_string(), index);
    Ok(())
}

fn int_size(ty: u8) -> Result<usize> {
    match ty {
        TYPE_INT8 => Ok(1),
        TYPE_INT16 => Ok(2),
        TYPE_INT32 => Ok(4),
        _ => Err(SeqError::malformed(FORMAT, format!("expected an integer type, found {ty}"))),
    }
}

/// The `(missing, end-of-vector)` sentinels for an integer type.
fn int_sentinels(ty: u8) -> (i64, i64) {
    match ty {
        TYPE_INT8 => (-128, -127),
        TYPE_INT16 => (-32768, -32767),
        _ => (i64::from(i32::MIN), i64::from(i32::MIN) + 1),
    }
}

fn int_read(raw: &[u8], ty: u8) -> i64 {
    match ty {
        TYPE_INT8 => i64::from(raw[0] as i8),
        TYPE_INT16 => i64::from(i16::from_le_bytes([raw[0], raw[1]])),
        _ => i64::from(i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]])),
    }
}

fn int_write(buf: &mut Vec<u8>, ty: u8, v: i64) {
    match ty {
        TYPE_INT8 => buf.push(v as i8 as u8),
        TYPE_INT16 => buf.extend_from_slice(&(v as i16).to_le_bytes()),
        _ => buf.extend_from_slice(&(v as i32).to_le_bytes()),
    }
}

/// The smallest integer type whose non-sentinel range holds every value.
fn int_type(vals: &[Option<i64>]) -> Result<u8> {
    let (mut lo, mut hi) = (0i64, 0i64);
    for v in vals.iter().flatten() {
        lo = lo.min(*v);
        hi = hi.max(*v);
    }
    if lo >= -120 && hi <= 127 {
        Ok(TYPE_INT8)
    } else if lo >= -32_760 && hi <= 32_767 {
        Ok(TYPE_INT16)
    } else if lo >= -2_147_483_640 && hi <= 2_147_483_647 {
        Ok(TYPE_INT32)
    } else {
        Err(SeqError::malformed(FORMAT, "integer value outside the BCF range"))
    }
}

/// Cursor over a record section. All reads are bounds-checked.
struct Fields<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Fields<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.data.len())
            .ok_or_else(|| SeqError::malformed(FORMAT, "record section is truncated"))?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn take_elems(&mut self, count: usize, size: usize) -> Result<&'a [u8]> {
        let n = count
            .checked_mul(size)
            .ok_or_else(|| SeqError::malformed(FORMAT, "value count overflows"))?;
        self.take(n)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }

    /// A typed integer scalar: a descriptor with count 1, then one value.
    fn typed_scalar_int(&mut self) -> Result<i64> {
        let byte = self.u8()?;
        let (ty, count) = (byte & 0x0f, byte >> 4);
        if count != 1 {
            return Err(SeqError::malformed(FORMAT, "expected a typed integer scalar"));
        }
        let size = int_size(ty)?;
        let v = int_read(self.take(size)?, ty);
        let (missing, eov) = int_sentinels(ty);
        if v == missing || v == eov {
            return Err(SeqError::malformed(FORMAT, "typed scalar is missing"));
        }
        Ok(v)
    }

    /// A type descriptor, returning `(type, count)`.
    fn descriptor(&mut self) -> Result<(u8, usize)> {
        let byte = self.u8()?;
        let ty = byte & 0x0f;
        let count = byte >> 4;
        if count == COUNT_ESCAPE {
            let n = usize::try_from(self.typed_scalar_int()?)
                .map_err(|_| SeqError::malformed(FORMAT, "negative value count"))?;
            return Ok((ty, n));
        }
        Ok((ty, usize::from(count)))
    }

    /// Integer values. Stops at the end-of-vector sentinel; `None` marks a missing value.
    fn ints(&mut self, ty: u8, count: usize) -> Result<Vec<Option<i64>>> {
        if count == 0 {
            return Ok(Vec::new());
        }
        let size = int_size(ty)?;
        let raw = self.take_elems(count, size)?;
        let (missing, eov) = int_sentinels(ty);
        let mut out = Vec::with_capacity(count);
        for chunk in raw.chunks_exact(size) {
            let v = int_read(chunk, ty);
            if v == eov {
                break;
            }
            out.push((v != missing).then_some(v));
        }
        Ok(out)
    }

    fn floats(&mut self, count: usize) -> Result<Vec<Option<f32>>> {
        let raw = self.take_elems(count, 4)?;
        let mut out = Vec::with_capacity(count);
        for chunk in raw.chunks_exact(4) {
            let bits = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            if bits == FLOAT_EOV {
                break;
            }
            out.push(if bits == FLOAT_MISSING { None } else { Some(f32::from_bits(bits)) });
        }
        Ok(out)
    }

    /// Characters, with trailing NUL padding removed.
    fn chars(&mut self, count: usize) -> Result<String> {
        let raw = self.take(count)?;
        let end = raw.iter().rposition(|&b| b != 0).map_or(0, |p| p + 1);
        String::from_utf8(raw[..end].to_vec())
            .map_err(|_| SeqError::malformed(FORMAT, "string value is not UTF-8"))
    }

    /// A typed string such as ID, REF or ALT. An empty string is a missing value.
    fn typed_string(&mut self) -> Result<String> {
        match self.descriptor()? {
            (TYPE_CHAR, count) => self.chars(count),
            (TYPE_MISSING, 0) => Ok(String::new()),
            _ => Err(SeqError::malformed(FORMAT, "expected a typed string")),
        }
    }
}

/// Renders a genotype's allele values, such as `0|1` or `./.`.
fn genotype_text(vals: &[Option<i64>]) -> String {
    let mut s = String::new();
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            s.push(if v.is_some_and(|x| x & 1 == 1) { '|' } else { '/' });
        }
        match v {
            Some(x) if (x >> 1) > 0 => {
                let _ = write!(s, "{}", (x >> 1) - 1);
            }
            _ => s.push('.'),
        }
    }
    if s.is_empty() {
        s.push('.');
    }
    s
}

/// Joins values with commas, rendering missing values as `.`.
fn join<T: Display>(vals: &[Option<T>]) -> String {
    let mut s = String::new();
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        match v {
            Some(x) => {
                let _ = write!(s, "{x}");
            }
            None => s.push('.'),
        }
    }
    s
}

/// Reads one INFO or FORMAT value of `count` elements as VCF text. `is_gt` selects the
/// genotype rendering for integer data. Returns an empty string when there is no value.
fn read_cell(c: &mut Fields<'_>, ty: u8, count: usize, is_gt: bool) -> Result<String> {
    match ty {
        TYPE_CHAR => c.chars(count),
        TYPE_MISSING if count == 0 => Ok(String::new()),
        TYPE_MISSING => Err(SeqError::malformed(FORMAT, "missing type with a nonzero count")),
        TYPE_FLOAT => Ok(join(&c.floats(count)?)),
        _ => {
            let vals = c.ints(ty, count)?;
            Ok(if is_gt { genotype_text(&vals) } else { join(&vals) })
        }
    }
}

fn or_dot(s: String) -> String {
    if s.is_empty() { ".".to_string() } else { s }
}

/// Decodes the shared and individual sections of one record into a VCF line.
fn decode_record(header: &BcfHeader, shared: &[u8], indiv: &[u8]) -> Result<String> {
    let mut c = Fields::new(shared);
    let chrom_id = usize::try_from(c.i32()?)
        .map_err(|_| SeqError::malformed(FORMAT, "record has a negative CHROM index"))?;
    let chrom = header
        .contig_name(chrom_id)
        .ok_or_else(|| SeqError::malformed(FORMAT, format!("contig index {chrom_id} undefined")))?;
    let pos0 = c.i32()?;
    let _rlen = c.i32()?;
    let qual_bits = c.u32()?;
    let allele_info = c.u32()?;
    let fmt_samples = c.u32()?;
    let n_allele = (allele_info >> 16) as usize;
    let n_info = (allele_info & 0xffff) as usize;
    let n_fmt = (fmt_samples >> 24) as usize;
    let n_sample = (fmt_samples & 0x00ff_ffff) as usize;
    if n_allele == 0 {
        return Err(SeqError::malformed(FORMAT, "record has no REF allele"));
    }
    if n_sample != header.n_samples() {
        return Err(SeqError::malformed(
            FORMAT,
            format!("record has {n_sample} samples, header has {}", header.n_samples()),
        ));
    }
    let pos = u64::try_from(pos0)
        .map_err(|_| SeqError::malformed(FORMAT, "negative POS"))?
        .checked_add(1)
        .ok_or_else(|| SeqError::malformed(FORMAT, "POS overflows"))?;

    let id = c.typed_string()?;
    let reference = c.typed_string()?;
    let mut alts = Vec::with_capacity(n_allele - 1);
    for _ in 1..n_allele {
        alts.push(c.typed_string()?);
    }
    let (ty, count) = c.descriptor()?;
    let mut filters = Vec::new();
    for id in c.ints(ty, count)?.into_iter().flatten() {
        filters.push(header.string_name(id)?);
    }

    let mut info_parts = Vec::with_capacity(n_info);
    for _ in 0..n_info {
        let key = header.string_name(c.typed_scalar_int()?)?;
        let (ty, count) = c.descriptor()?;
        let value = read_cell(&mut c, ty, count, false)?;
        info_parts.push(if value.is_empty() { key.to_string() } else { format!("{key}={value}") });
    }

    let mut keys = Vec::with_capacity(n_fmt);
    let mut cells: Vec<Vec<String>> = vec![Vec::with_capacity(n_fmt); n_sample];
    let mut f = Fields::new(indiv);
    for _ in 0..n_fmt {
        let key = header.string_name(f.typed_scalar_int()?)?;
        let (ty, count) = f.descriptor()?;
        for sample in &mut cells {
            let text = read_cell(&mut f, ty, count, key == "GT")?;
            sample.push(or_dot(text));
        }
        keys.push(key);
    }

    let qual = if qual_bits == FLOAT_MISSING {
        ".".to_string()
    } else {
        f32::from_bits(qual_bits).to_string()
    };
    let mut cols = vec![
        chrom.to_string(),
        pos.to_string(),
        or_dot(id),
        reference,
        or_dot(alts.join(",")),
        qual,
        or_dot(filters.join(";")),
        if info_parts.is_empty() { ".".to_string() } else { info_parts.join(";") },
    ];
    if n_fmt > 0 || n_sample > 0 {
        cols.push(if keys.is_empty() { ".".to_string() } else { keys.join(":") });
        for sample in &cells {
            cols.push(if sample.is_empty() { ".".to_string() } else { sample.join(":") });
        }
    }
    Ok(cols.join("\t"))
}

fn push_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_descriptor(buf: &mut Vec<u8>, ty: u8, count: usize) -> Result<()> {
    if count < usize::from(COUNT_ESCAPE) {
        buf.push(((count as u8) << 4) | ty);
        Ok(())
    } else {
        buf.push((COUNT_ESCAPE << 4) | ty);
        put_int_scalar(buf, count as i64)
    }
}

fn put_int_scalar(buf: &mut Vec<u8>, v: i64) -> Result<()> {
    let ty = int_type(&[Some(v)])?;
    buf.push((1 << 4) | ty);
    int_write(buf, ty, v);
    Ok(())
}

fn put_string(buf: &mut Vec<u8>, s: &str) -> Result<()> {
    // An empty string is a zero-length char array, which is what bcftools writes for
    // a missing ID. Type 0 (null) is not valid for these fields.
    put_descriptor(buf, TYPE_CHAR, s.len())?;
    buf.extend_from_slice(s.as_bytes());
    Ok(())
}

/// Writes integers padded to `count` elements with the end-of-vector sentinel.
fn put_ints(buf: &mut Vec<u8>, ty: u8, vals: &[Option<i64>], count: usize) {
    let (missing, eov) = int_sentinels(ty);
    for i in 0..count {
        let v = match vals.get(i) {
            Some(Some(x)) => *x,
            Some(None) => missing,
            None => eov,
        };
        int_write(buf, ty, v);
    }
}

fn put_int_vector(buf: &mut Vec<u8>, vals: &[Option<i64>]) -> Result<()> {
    let ty = int_type(vals)?;
    put_descriptor(buf, ty, vals.len())?;
    put_ints(buf, ty, vals, vals.len());
    Ok(())
}

fn put_floats(buf: &mut Vec<u8>, vals: &[Option<f32>], count: usize) {
    for i in 0..count {
        let bits = match vals.get(i) {
            Some(Some(x)) => x.to_bits(),
            Some(None) => FLOAT_MISSING,
            None => FLOAT_EOV,
        };
        push_u32(buf, bits);
    }
}

fn parse_ints(s: &str) -> Result<Vec<Option<i64>>> {
    s.split(',')
        .map(|t| {
            if t == "." {
                Ok(None)
            } else {
                t.parse::<i64>()
                    .map(Some)
                    .map_err(|_| SeqError::malformed(FORMAT, format!("not an integer: {t}")))
            }
        })
        .collect()
}

fn parse_floats(s: &str) -> Result<Vec<Option<f32>>> {
    s.split(',')
        .map(|t| {
            if t == "." {
                Ok(None)
            } else {
                t.parse::<f32>()
                    .map(Some)
                    .map_err(|_| SeqError::malformed(FORMAT, format!("not a number: {t}")))
            }
        })
        .collect()
}

/// Parses a genotype such as `0|1`, `./.` or `2`. Each allele becomes `(index + 1) << 1`
/// with the phase bit of the separator before it. A missing allele is 0.
fn parse_genotype(s: &str) -> Result<Vec<Option<i64>>> {
    if s == "." || s.is_empty() {
        return Ok(vec![Some(0)]);
    }
    let mut tokens: Vec<(&str, i64)> = Vec::new();
    let mut phase = 0;
    let mut start = 0;
    for (i, c) in s.char_indices() {
        if c == '/' || c == '|' {
            tokens.push((&s[start..i], phase));
            phase = i64::from(c == '|');
            start = i + 1;
        }
    }
    tokens.push((&s[start..], phase));
    tokens
        .into_iter()
        .map(|(tok, bit)| {
            let allele = if tok == "." {
                0
            } else {
                tok.parse::<i64>()
                    .ok()
                    .filter(|a| *a >= 0)
                    .map(|a| a + 1)
                    .ok_or_else(|| SeqError::malformed(FORMAT, format!("bad allele: {tok}")))?
            };
            Ok(Some((allele << 1) | bit))
        })
        .collect()
}

fn encode_info(buf: &mut Vec<u8>, declared: &str, key: &str, value: Option<&str>) -> Result<()> {
    match (declared, value) {
        (_, None) => put_descriptor(buf, TYPE_MISSING, 0),
        ("Flag", Some(_)) => {
            Err(SeqError::malformed(FORMAT, format!("flag INFO {key} has a value")))
        }
        ("Integer", Some(v)) => put_int_vector(buf, &parse_ints(v)?),
        ("Float", Some(v)) => {
            let vals = parse_floats(v)?;
            put_descriptor(buf, TYPE_FLOAT, vals.len())?;
            put_floats(buf, &vals, vals.len());
            Ok(())
        }
        ("String" | "Character", Some(v)) => put_string(buf, v),
        (other, _) => {
            Err(SeqError::malformed(FORMAT, format!("INFO {key} has unsupported Type {other}")))
        }
    }
}

/// Writes one integer field across samples with a shared type and count.
fn put_int_block(buf: &mut Vec<u8>, per: &[Vec<Option<i64>>]) -> Result<()> {
    let count = per.iter().map(Vec::len).max().unwrap_or(0);
    let all: Vec<Option<i64>> = per.iter().flat_map(|v| v.iter().copied()).collect();
    let ty = int_type(&all)?;
    put_descriptor(buf, ty, count)?;
    for v in per {
        put_ints(buf, ty, v, count);
    }
    Ok(())
}

fn put_float_block(buf: &mut Vec<u8>, per: &[Vec<Option<f32>>]) -> Result<()> {
    let count = per.iter().map(Vec::len).max().unwrap_or(0);
    put_descriptor(buf, TYPE_FLOAT, count)?;
    for v in per {
        put_floats(buf, v, count);
    }
    Ok(())
}

fn put_char_block(buf: &mut Vec<u8>, texts: &[&str]) -> Result<()> {
    let bytes: Vec<&[u8]> =
        texts.iter().map(|t| if *t == "." { &b""[..] } else { t.as_bytes() }).collect();
    let count = bytes.iter().map(|b| b.len()).max().unwrap_or(0);
    put_descriptor(buf, TYPE_CHAR, count)?;
    for b in bytes {
        buf.extend_from_slice(b);
        let padded = buf.len() + count - b.len();
        buf.resize(padded, 0);
    }
    Ok(())
}

/// Encodes one FORMAT field (column `field`) across all samples.
fn encode_format(
    buf: &mut Vec<u8>,
    header: &BcfHeader,
    key: &str,
    samples: &[&str],
    field: usize,
) -> Result<()> {
    let texts: Vec<&str> = samples
        .iter()
        .map(|s| s.split(':').nth(field).filter(|t| !t.is_empty()).unwrap_or("."))
        .collect();
    if key == "GT" {
        let per = texts.iter().map(|t| parse_genotype(t)).collect::<Result<Vec<_>>>()?;
        return put_int_block(buf, &per);
    }
    match header.format_type(key)? {
        "Integer" => {
            let per = texts.iter().map(|t| parse_ints(t)).collect::<Result<Vec<_>>>()?;
            put_int_block(buf, &per)
        }
        "Float" => {
            let per = texts.iter().map(|t| parse_floats(t)).collect::<Result<Vec<_>>>()?;
            put_float_block(buf, &per)
        }
        "String" | "Character" => put_char_block(buf, &texts),
        other => {
            Err(SeqError::malformed(FORMAT, format!("FORMAT {key} has unsupported Type {other}")))
        }
    }
}

/// Encodes a record into the shared and individual sections.
fn encode_record(
    header: &BcfHeader,
    rec: &VcfRecord<'_>,
    shared: &mut Vec<u8>,
    indiv: &mut Vec<u8>,
) -> Result<()> {
    let samples: Vec<&str> = rec.samples().collect();
    if samples.len() != header.n_samples() {
        return Err(SeqError::malformed(
            FORMAT,
            format!("record has {} samples, header has {}", samples.len(), header.n_samples()),
        ));
    }
    let chrom = header.contig_id(rec.chrom()).ok_or_else(|| {
        SeqError::malformed(FORMAT, format!("contig {} not in header", rec.chrom()))
    })?;
    let chrom =
        i32::try_from(chrom).map_err(|_| SeqError::malformed(FORMAT, "too many contigs"))?;
    let pos0 =
        i32::try_from(rec.pos() - 1).map_err(|_| SeqError::malformed(FORMAT, "POS too large"))?;
    let info: Vec<(&str, Option<&str>)> = rec.info().collect();
    let alts: Vec<&str> = rec.alts().collect();
    let filters: Vec<&str> = rec.filters().collect();
    let keys: Vec<&str> = rec.format_keys().collect();
    let end = info
        .iter()
        .find(|(k, _)| *k == "END")
        .and_then(|(_, v)| v.and_then(|v| v.parse::<i64>().ok()));
    let rlen = end.map_or(rec.reference().len() as i64, |e| e - i64::from(pos0));
    let rlen = i32::try_from(rlen).map_err(|_| SeqError::malformed(FORMAT, "END out of range"))?;
    let n_info = u32::try_from(info.len())
        .ok()
        .filter(|n| *n <= 0xffff)
        .ok_or_else(|| SeqError::malformed(FORMAT, "too many INFO fields"))?;
    let n_fmt = u32::try_from(keys.len())
        .ok()
        .filter(|n| *n <= 0xff)
        .ok_or_else(|| SeqError::malformed(FORMAT, "too many FORMAT fields"))?;
    let n_allele = (alts.len() + 1) as u32;
    let qual = rec.qual().map_or(FLOAT_MISSING, f32::to_bits);

    shared.clear();
    indiv.clear();
    shared.extend_from_slice(&chrom.to_le_bytes());
    shared.extend_from_slice(&pos0.to_le_bytes());
    shared.extend_from_slice(&rlen.to_le_bytes());
    push_u32(shared, qual);
    push_u32(shared, (n_allele << 16) | n_info);
    push_u32(shared, (n_fmt << 24) | samples.len() as u32);
    put_string(shared, rec.id().unwrap_or(""))?;
    put_string(shared, rec.reference())?;
    for alt in &alts {
        put_string(shared, alt)?;
    }
    let filter_ids =
        filters.iter().map(|f| header.string_id(f).map(Some)).collect::<Result<Vec<_>>>()?;
    put_int_vector(shared, &filter_ids)?;
    for (key, value) in &info {
        put_int_scalar(shared, header.string_id(key)?)?;
        encode_info(shared, header.info_type(key)?, key, *value)?;
    }
    for (field, key) in keys.iter().enumerate() {
        put_int_scalar(indiv, header.string_id(key)?)?;
        encode_format(indiv, header, key, &samples, field)?;
    }
    Ok(())
}

/// Reads the `BCF` magic, version, and header text.
fn read_header<R: Read>(r: &mut R) -> Result<BcfHeader> {
    let mut magic = [0u8; 5];
    read_full(r, &mut magic, "magic")?;
    if &magic[..3] != MAGIC {
        return Err(SeqError::malformed(FORMAT, "missing BCF magic"));
    }
    if magic[3] != 2 || !matches!(magic[4], 1 | 2) {
        return Err(SeqError::Unsupported(format!("BCF version {}.{}", magic[3], magic[4])));
    }
    let mut len = [0u8; 4];
    read_full(r, &mut len, "header length")?;
    let l_text = u32::from_le_bytes(len) as usize;
    if l_text > MAX_HEADER_TEXT {
        return Err(SeqError::malformed(FORMAT, "header text is too large"));
    }
    let text = read_bounded(r, l_text).map_err(|e| {
        if e.kind() == std::io::ErrorKind::UnexpectedEof {
            SeqError::malformed(FORMAT, "truncated header text")
        } else {
            SeqError::Io(e)
        }
    })?;
    let end = text.iter().position(|&b| b == 0).unwrap_or(text.len());
    let text = std::str::from_utf8(&text[..end])
        .map_err(|_| SeqError::malformed(FORMAT, "header text is not UTF-8"))?;
    BcfHeader::parse(text)
}

/// Records up to this size reuse the record buffer. See `read_raw_record`.
const INLINE_RECORD_BYTES: usize = 1 << 16;

fn read_full<R: Read>(r: &mut R, buf: &mut [u8], what: &str) -> Result<()> {
    r.read_exact(buf).map_err(|e| {
        if e.kind() == std::io::ErrorKind::UnexpectedEof {
            SeqError::malformed(FORMAT, format!("truncated {what}"))
        } else {
            SeqError::Io(e)
        }
    })
}

/// Reads one record into `buf` as shared then individual bytes. Returns the length of
/// the shared part, or `None` at a clean end of input.
pub(crate) fn read_raw_record<R: Read>(r: &mut R, buf: &mut Vec<u8>) -> Result<Option<usize>> {
    let mut prefix = [0u8; 8];
    let mut filled = 0;
    while filled < prefix.len() {
        match r.read(&mut prefix[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(SeqError::Io(e)),
        }
    }
    if filled == 0 {
        return Ok(None);
    }
    if filled < prefix.len() {
        return Err(SeqError::malformed(FORMAT, "truncated record length prefix"));
    }
    let l_shared = u32::from_le_bytes([prefix[0], prefix[1], prefix[2], prefix[3]]) as usize;
    let l_indiv = u32::from_le_bytes([prefix[4], prefix[5], prefix[6], prefix[7]]) as usize;
    let total = l_shared
        .checked_add(l_indiv)
        .filter(|&t| t <= MAX_RECORD_BYTES)
        .ok_or_else(|| SeqError::malformed(FORMAT, "record length is too large"))?;
    if total <= INLINE_RECORD_BYTES {
        buf.clear();
        buf.resize(total, 0);
        read_full(r, buf, "record")?;
    } else {
        // A corrupt length must not force a large zero-filled allocation.
        *buf = read_bounded(r, total).map_err(|e| {
            if e.kind() == std::io::ErrorKind::UnexpectedEof {
                SeqError::malformed(FORMAT, "truncated record")
            } else {
                SeqError::Io(e)
            }
        })?;
    }
    Ok(Some(l_shared))
}

/// Reads a BCF file and yields records as VCF views.
#[derive(Debug)]
pub struct BcfReader<R> {
    inner: R,
    header: BcfHeader,
    buf: Vec<u8>,
    line: String,
    n_records: u64,
}

impl<R: Read> BcfReader<R> {
    /// Reads the BCF header from an uncompressed stream. Use [`BcfReader::new_bgzf`] for
    /// standard BGZF-compressed BCF.
    ///
    /// # Errors
    /// Returns an I/O error, or [`SeqError::Malformed`] or [`SeqError::Unsupported`] for
    /// an invalid or unsupported header.
    pub fn new(mut inner: R) -> Result<Self> {
        let header = read_header(&mut inner)?;
        Ok(Self { inner, header, buf: Vec::new(), line: String::new(), n_records: 0 })
    }

    /// The parsed BCF header.
    #[must_use]
    pub fn header(&self) -> &BcfHeader {
        &self.header
    }

    /// Reads the next record, or `None` at end of input. The record borrows the reader's
    /// buffer, so drop it before the next call.
    ///
    /// # Errors
    /// Returns an I/O error, or [`SeqError::Malformed`] for a corrupt record.
    pub fn next_record(&mut self) -> Result<Option<VcfRecord<'_>>> {
        let Some(split) = read_raw_record(&mut self.inner, &mut self.buf)? else {
            return Ok(None);
        };
        self.n_records += 1;
        self.line = decode_record(&self.header, &self.buf[..split], &self.buf[split..])?;
        VcfRecord::parse(&self.line, self.header.n_samples(), self.n_records).map(Some)
    }
}

impl<S: Read> BcfReader<BgzfReader<S>> {
    /// Reads a BGZF-compressed BCF file, the usual form.
    ///
    /// # Errors
    /// Same as [`BcfReader::new`], plus BGZF decoding errors.
    pub fn new_bgzf(inner: S) -> Result<Self> {
        Self::new(BgzfReader::new(inner))
    }
}

/// Writes BCF records for a given header.
#[derive(Debug)]
pub struct BcfWriter<W> {
    inner: W,
    header: BcfHeader,
    shared: Vec<u8>,
    indiv: Vec<u8>,
}

impl<W: Write> BcfWriter<W> {
    /// Writes the uncompressed magic and header.
    ///
    /// # Errors
    /// Returns an I/O error.
    pub fn new(mut inner: W, header: &BcfHeader) -> Result<Self> {
        let text = header.vcf().text().as_bytes();
        inner.write_all(b"BCF\x02\x02")?;
        inner.write_all(&((text.len() + 1) as u32).to_le_bytes())?;
        inner.write_all(text)?;
        inner.write_all(&[0])?;
        Ok(Self { inner, header: header.clone(), shared: Vec::new(), indiv: Vec::new() })
    }

    /// Encodes and writes one record.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] if the record does not fit the header (unknown
    /// contig, filter, INFO or FORMAT key, or sample count), or an I/O error.
    pub fn write_record(&mut self, record: &VcfRecord<'_>) -> Result<()> {
        encode_record(&self.header, record, &mut self.shared, &mut self.indiv)?;
        self.inner.write_all(&(self.shared.len() as u32).to_le_bytes())?;
        self.inner.write_all(&(self.indiv.len() as u32).to_le_bytes())?;
        self.inner.write_all(&self.shared)?;
        self.inner.write_all(&self.indiv)?;
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

impl<W: Write> BcfWriter<BgzfWriter<W>> {
    /// Writes BGZF-compressed BCF. Call [`BcfWriter::new_bgzf`] to create one.
    ///
    /// # Errors
    /// Returns an I/O error.
    pub fn new_bgzf(inner: W, header: &BcfHeader) -> Result<Self> {
        Self::new(BgzfWriter::new(inner), header)
    }

    /// Writes the BGZF EOF block and returns the underlying writer.
    ///
    /// # Errors
    /// Returns an I/O error from the final block.
    pub fn finish(self) -> Result<W> {
        self.inner.finish()
    }
}

/// Builds a tabix-style index (CSI or TBI) for a BGZF-compressed, coordinate-sorted
/// BCF. Records are keyed by contig index, and the `CHROM`, `POS` and `rlen` fields are
/// read from the record.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for a corrupt file or unsorted records, or an error
/// from the underlying reader.
pub fn build_bcf_index<R: Read>(input: R, kind: IndexKind) -> Result<TabixIndex> {
    let mut bgzf = BgzfReader::new(input);
    let header = read_header(&mut bgzf)?;
    let mut builder = IndexBuilder::new(kind, None);
    builder.ensure_refs(header.n_contigs());
    let mut buf = Vec::new();
    loop {
        let begin = bgzf.virtual_position();
        let Some(split) = read_raw_record(&mut bgzf, &mut buf)? else {
            break;
        };
        let end = bgzf.virtual_position();
        let (tid, beg0, end0) = record_extent(&buf[..split])?;
        builder.add(tid, beg0, end0, begin, end)?;
    }
    Ok(builder.finish())
}

/// The contig index, 0-based start, and exclusive end of a record's shared section.
fn record_extent(shared: &[u8]) -> Result<(usize, u64, u64)> {
    let mut c = Fields::new(shared);
    let tid = usize::try_from(c.i32()?)
        .map_err(|_| SeqError::malformed(FORMAT, "record has a negative CHROM index"))?;
    let pos = u64::try_from(c.i32()?).map_err(|_| SeqError::malformed(FORMAT, "negative POS"))?;
    let rlen = c.i32()?.max(1) as u64;
    Ok((tid, pos, pos + rlen))
}

/// Returns VCF text lines for the records of a BCF file that overlap a region.
///
/// `reader` must be a BGZF reader over the same file as `index` (`Read + Seek`).
/// `start` and `end` are 1-based and inclusive. Records come back in file order.
///
/// # Errors
/// Returns [`SeqError::NotFound`] for an unknown contig, [`SeqError::Malformed`] for a
/// bad region or corrupt data, or an I/O error.
pub fn fetch_lines<R: Read + Seek>(
    reader: &mut BgzfReader<R>,
    header: &BcfHeader,
    index: &TabixIndex,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<String>> {
    let tid =
        header.contig_id(chrom).ok_or_else(|| SeqError::NotFound(format!("contig {chrom}")))?;
    if start < 1 || end < start {
        return Err(SeqError::malformed(FORMAT, "region start must be >= 1 and end >= start"));
    }
    let (q_beg, q_end) = (start - 1, end);
    let mut out = Vec::new();
    let mut buf = Vec::new();
    for chunk in index.query_by_id(tid, start, end) {
        reader.seek_virtual(chunk.beg)?;
        while reader.virtual_position() < chunk.end {
            let Some(split) = read_raw_record(reader, &mut buf)? else {
                break;
            };
            let (rtid, beg0, end0) = record_extent(&buf[..split])?;
            if rtid == tid && beg0 < q_end && end0 > q_beg {
                out.push(decode_record(header, &buf[..split], &buf[split..])?);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tabix::{IndexKind, TabixIndex, test_support};
    use std::io::Cursor;

    const HEADER: &str = "##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP\">
##INFO=<ID=RANGE,Number=.,Type=Integer,Description=\"Range probe\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
";

    const RECORDS: [&str; 4] = [
        "chr1\t100\trs1\tA\tG\t30\tPASS\tDP=12;AF=0.5;GENE=abc;DB\tGT:DP:GQ:FT\t0/1:9:44.5:ok\t1|1:.:.:.",
        "chr1\t200\t.\tC\tA,T\t.\tq10\tDP=.;AF=0.1,0.2;RANGE=-120,127,-121,128,32767,32768,-2147483640,2147483647\tGT:DP:GQ:FT\t./.:3:.:q10\t0/0:.:10:.",
        "chr2\t5\t.\tT\t.\t50\t.\t.\tGT\t0/0\t.",
        "chr2\t300\trs9\tGATC\tG\t99.5\tPASS\tDP=40;GENE=xyz;DB;RANGE=5\tGT:DP\t1/0:12\t0|0:40",
    ];

    fn vcf_text(records: &[&str]) -> String {
        let mut s = HEADER.to_string();
        for r in records {
            s.push_str(r);
            s.push('\n');
        }
        s
    }

    #[test]
    fn vcf_to_bcf_round_trip_preserves_text() {
        let text = vcf_text(&RECORDS);
        let mut reader = crate::vcf::VcfReader::new(Cursor::new(text.as_bytes())).unwrap();
        let header = BcfHeader::from_vcf(reader.header().clone()).unwrap();
        let mut bytes = Vec::new();
        {
            let mut writer = BcfWriter::new_bgzf(&mut bytes, &header).unwrap();
            while let Some(rec) = reader.next_record().unwrap() {
                writer.write_record(&rec).unwrap();
            }
            writer.finish().unwrap();
        }
        assert!(!bytes.windows(3).any(|w| w == b"DP="), "payload must be compressed");

        let mut back = BcfReader::new_bgzf(Cursor::new(bytes)).unwrap();
        assert_eq!(back.header().vcf().samples(), ["s1", "s2"]);
        let mut lines = Vec::new();
        while let Some(rec) = back.next_record().unwrap() {
            lines.push(rec.line().to_string());
        }
        assert_eq!(lines, RECORDS);
    }

    #[test]
    fn decoded_records_expose_typed_fields() {
        let text = vcf_text(&RECORDS);
        let mut reader = crate::vcf::VcfReader::new(Cursor::new(text.as_bytes())).unwrap();
        let header = BcfHeader::from_vcf(reader.header().clone()).unwrap();
        let mut bytes = Vec::new();
        let mut writer = BcfWriter::new(&mut bytes, &header).unwrap();
        while let Some(rec) = reader.next_record().unwrap() {
            writer.write_record(&rec).unwrap();
        }
        writer.into_inner().unwrap();

        let mut back = BcfReader::new(Cursor::new(bytes)).unwrap();
        let first = back.next_record().unwrap().unwrap();
        assert_eq!(first.pos(), 100);
        assert_eq!(first.qual(), Some(30.0));
        assert_eq!(first.filters().collect::<Vec<_>>(), ["PASS"]);
        assert!(first.info().any(|(k, v)| k == "DB" && v.is_none()));
        let second = back.next_record().unwrap().unwrap();
        assert_eq!(second.qual(), None);
        assert_eq!(second.alts().collect::<Vec<_>>(), ["A", "T"]);
        assert_eq!(second.filters().collect::<Vec<_>>(), ["q10"]);
    }

    #[test]
    fn dictionary_places_pass_first_and_honours_idx() {
        let text = "##fileformat=VCFv4.3
##contig=<ID=chrX,IDX=3>
##contig=<ID=chrY>
##FILTER=<ID=lowq,Description=\"low\">
##INFO=<ID=ZZ,Number=1,Type=Integer,Description=\"z\",IDX=7>
##INFO=<ID=AA,Number=1,Type=Integer,Description=\"a\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";
        let h = BcfHeader::parse(text).unwrap();
        assert_eq!(h.string_id("PASS").unwrap(), 0);
        assert_eq!(h.string_id("ZZ").unwrap(), 7);
        assert_eq!(h.string_id("lowq").unwrap(), 1);
        assert_eq!(h.string_id("AA").unwrap(), 2);
        assert_eq!(h.contig_id("chrX"), Some(3));
        assert_eq!(h.contig_id("chrY"), Some(0));
        assert!(BcfHeader::parse(
            "##fileformat=VCFv4.3\n##INFO=<ID=A,Number=1,Type=Integer,Description=\"a\",IDX=1>\n\
             ##INFO=<ID=B,Number=1,Type=Integer,Description=\"b\",IDX=1>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        .is_err());
    }

    #[test]
    fn integer_widths_pick_smallest_fitting_type() {
        assert_eq!(int_type(&[Some(-120), Some(127)]).unwrap(), TYPE_INT8);
        assert_eq!(int_type(&[Some(128)]).unwrap(), TYPE_INT16);
        assert_eq!(int_type(&[Some(-121)]).unwrap(), TYPE_INT16);
        assert_eq!(int_type(&[Some(32_768)]).unwrap(), TYPE_INT32);
        assert_eq!(int_type(&[None]).unwrap(), TYPE_INT8);
        assert!(int_type(&[Some(i64::MAX)]).is_err());
    }

    #[test]
    fn genotype_encoding_round_trips() {
        for text in ["0/1", "1|1", "./.", "0", ".", "2|0|1"] {
            let vals = parse_genotype(text).unwrap();
            assert_eq!(genotype_text(&vals), text);
        }
    }

    #[test]
    fn rejects_truncated_and_bad_input() {
        let text = vcf_text(&RECORDS);
        let mut reader = crate::vcf::VcfReader::new(Cursor::new(text.as_bytes())).unwrap();
        let header = BcfHeader::from_vcf(reader.header().clone()).unwrap();
        let mut bytes = Vec::new();
        let mut writer = BcfWriter::new(&mut bytes, &header).unwrap();
        while let Some(rec) = reader.next_record().unwrap() {
            writer.write_record(&rec).unwrap();
        }
        writer.into_inner().unwrap();

        let cut = bytes.len() - 5;
        let mut back = BcfReader::new(Cursor::new(bytes[..cut].to_vec())).unwrap();
        let mut failed = false;
        for _ in 0..RECORDS.len() {
            match back.next_record() {
                Ok(Some(_)) => {}
                Ok(None) => break,
                Err(_) => {
                    failed = true;
                    break;
                }
            }
        }
        assert!(failed, "truncated record must be an error");
        assert!(BcfReader::new(Cursor::new(b"VCF\x02\x02".to_vec())).is_err());
        assert!(BcfReader::new(Cursor::new(b"BCF\x02\x03\x00\x00\x00\x00".to_vec())).is_err());
    }

    #[test]
    fn writer_rejects_undeclared_keys() {
        let header = BcfHeader::parse(HEADER).unwrap();
        let mut out = Vec::new();
        let mut writer = BcfWriter::new(&mut out, &header).unwrap();
        let line = "chr1\t5\t.\tC\tG\t.\t.\tXX=1\tGT\t0/1\t0/1";
        let rec = VcfRecord::parse(line, 2, 1).unwrap();
        assert!(writer.write_record(&rec).is_err());
    }

    #[test]
    fn region_fetch_matches_brute_force_on_bcf() {
        let text = test_support::synthetic_vcf(4000, 7);
        let mut reader = crate::vcf::VcfReader::new(Cursor::new(text.as_bytes())).unwrap();
        let header = BcfHeader::from_vcf(reader.header().clone()).unwrap();
        let mut bytes = Vec::new();
        {
            let mut writer = BcfWriter::new_bgzf(&mut bytes, &header).unwrap();
            while let Some(rec) = reader.next_record().unwrap() {
                writer.write_record(&rec).unwrap();
            }
            writer.finish().unwrap();
        }
        let index: TabixIndex =
            build_bcf_index(&bytes[..], IndexKind::Csi { min_shift: 14, depth: 5 }).unwrap();
        let all = test_support::data_lines(&text);
        let mut rng = test_support::Lcg::new(99);
        for _ in 0..60 {
            let chrom = if rng.next().is_multiple_of(2) { "chr1" } else { "chr2" };
            let start = 1 + rng.next() % 1_500_000;
            let end = start + rng.next() % 40_000;
            let expected = test_support::brute_force(&all, chrom, start, end);
            let mut src = BgzfReader::new(Cursor::new(bytes.clone()));
            let got = fetch_lines(&mut src, &header, &index, chrom, start, end).unwrap();

            assert_eq!(got, expected, "region {chrom}:{start}-{end}");
        }
    }
}
