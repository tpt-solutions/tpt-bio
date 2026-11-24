//! CRAM compression header: preservation map, data-series encodings and tag encodings.
//!
//! The compression header block holds three ITF8-prefixed maps, each with a byte
//! length, an entry count, then entries. Keys are two ASCII bytes. Preservation
//! map values are typed (booleans, or the five-byte substitution matrix and the
//! tag dictionary). Data-series and tag maps map a key to an encoding descriptor:
//! an ITF8 encoding id, an ITF8 parameter length and the parameter bytes.

use crate::cram::bytes::{Cursor, put_itf8};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM compression header";

/// Encoding ids from the CRAM specification.
pub mod encoding_id {
    /// No data.
    pub const NULL: i32 = 0;
    /// Bytes from an external block.
    pub const EXTERNAL: i32 = 1;
    /// Golomb-style Gamma code (bit stream).
    pub const GAMMA: i32 = 2;
    /// Canonical Huffman code (bit stream).
    pub const HUFFMAN: i32 = 3;
    /// Fixed-width bit field (bit stream).
    pub const BYTE_ARRAY_LEN: i32 = 4;
    /// Bytes up to a stop byte.
    pub const BYTE_ARRAY_STOP: i32 = 5;
    /// Beta code (bit stream).
    pub const BETA: i32 = 6;
    /// Sub-exponential code (bit stream).
    pub const SUBEXP: i32 = 7;
}

/// A parsed encoding descriptor. Parameters are kept in their raw form and
/// decoded when the record codec is built.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Encoding {
    /// No data.
    Null,
    /// External block with the given content id.
    External {
        /// Content id of the external block.
        content_id: i32,
    },
    /// Canonical Huffman code.
    Huffman {
        /// Alphabet, in file order.
        symbols: Vec<i32>,
        /// Code length for each symbol.
        lengths: Vec<u32>,
    },
    /// Fixed-width integer read from the core stream.
    Beta {
        /// Offset subtracted from the raw value.
        offset: i32,
        /// Number of bits.
        bits: u32,
    },
    /// Gamma code (`offset`-adjusted).
    Gamma {
        /// Offset subtracted from the decoded value.
        offset: i32,
    },
    /// Byte array with a length encoding and a byte encoding.
    ByteArrayLen {
        /// Encoding of the array length.
        len: Box<Encoding>,
        /// Encoding of each byte.
        bytes: Box<Encoding>,
    },
    /// Byte array terminated by a stop byte.
    ByteArrayStop {
        /// The terminating byte value.
        stop: u8,
        /// Content id of the external block holding the bytes.
        content_id: i32,
    },
    /// Sub-exponential code.
    Subexp {
        /// Offset subtracted from the decoded value.
        offset: i32,
        /// The `k` parameter.
        k: i32,
    },
}

/// Preservation map values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Preservation {
    /// Whether read names are preserved (RN).
    pub read_names: bool,
    /// Whether alignment positions are coordinate-sorted (AP).
    pub ap_delta: bool,
    /// Whether reference sequences are stored in the file (RR).
    pub reference_required: bool,
    /// Substitution matrix (SM), five bytes.
    pub substitution_matrix: [u8; 5],
    /// Tag dictionary (TD): each entry is a list of 3-byte tag keys with types.
    pub tag_dictionary: Vec<Vec<[u8; 3]>>,
}

/// The parsed compression header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompressionHeader {
    /// Preservation map.
    pub preservation: Preservation,
    /// Data-series encodings keyed by their two-byte name (`RL`, `FN`, ...).
    pub data_series: Vec<([u8; 2], Encoding)>,
    /// Tag encodings keyed by the three-byte tag id (name and type), in dictionary order.
    pub tags: Vec<([u8; 3], Encoding)>,
}

impl CompressionHeader {
    /// Looks up a data-series encoding by its two-byte key.
    #[must_use]
    pub fn series(&self, key: &[u8; 2]) -> Option<&Encoding> {
        self.data_series.iter().find(|(k, _)| k == key).map(|(_, e)| e)
    }
}

fn read_key(c: &mut Cursor<'_>) -> Result<[u8; 2]> {
    let b = c.bytes(2)?;
    Ok([b[0], b[1]])
}

fn read_bool(c: &mut Cursor<'_>) -> Result<bool> {
    Ok(c.u8()? != 0)
}

/// Reads one encoding descriptor (id, parameter length, parameters).
///
/// # Errors
/// Returns [`SeqError::Malformed`] for truncated or inconsistent parameters, and
/// [`SeqError::Unsupported`] for an encoding id this reader does not implement.
pub fn read_encoding(c: &mut Cursor<'_>) -> Result<Encoding> {
    let id = c.itf8()?;
    let len = usize::try_from(c.itf8()?)
        .map_err(|_| SeqError::malformed(FORMAT, "negative encoding length"))?;
    let params = c.bytes(len)?;
    let mut p = Cursor::new(params);
    let enc = match id {
        encoding_id::NULL => Encoding::Null,
        encoding_id::EXTERNAL => Encoding::External { content_id: p.itf8()? },
        encoding_id::HUFFMAN => {
            let n = p.itf8()?;
            let symbols = read_itf8_array(&mut p, n)?;
            let m = p.itf8()?;
            let lengths = read_itf8_array(&mut p, m)?;
            let lengths = lengths
                .into_iter()
                .map(|l| {
                    u32::try_from(l)
                        .map_err(|_| SeqError::malformed(FORMAT, "negative code length"))
                })
                .collect::<Result<Vec<u32>>>()?;
            if symbols.len() != lengths.len() {
                return Err(SeqError::malformed(FORMAT, "Huffman symbol and length counts differ"));
            }
            Encoding::Huffman { symbols, lengths }
        }
        encoding_id::BETA => {
            let offset = p.itf8()?;
            let bits = p.itf8()?;
            Encoding::Beta {
                offset,
                bits: u32::try_from(bits)
                    .map_err(|_| SeqError::malformed(FORMAT, "negative BETA bits"))?,
            }
        }
        encoding_id::GAMMA => Encoding::Gamma { offset: p.itf8()? },
        encoding_id::BYTE_ARRAY_LEN => {
            let len = read_encoding(&mut p)?;
            let bytes = read_encoding(&mut p)?;
            Encoding::ByteArrayLen { len: Box::new(len), bytes: Box::new(bytes) }
        }
        encoding_id::BYTE_ARRAY_STOP => {
            let stop = p.u8()?;
            let content_id = p.itf8()?;
            Encoding::ByteArrayStop { stop, content_id }
        }
        encoding_id::SUBEXP => {
            let offset = p.itf8()?;
            let k = p.itf8()?;
            Encoding::Subexp { offset, k }
        }
        other => return Err(SeqError::Unsupported(format!("CRAM encoding id {other}"))),
    };
    Ok(enc)
}

fn read_itf8_array(c: &mut Cursor<'_>, n: i32) -> Result<Vec<i32>> {
    let n = usize::try_from(n).map_err(|_| SeqError::malformed(FORMAT, "negative array count"))?;
    let mut v = Vec::with_capacity(n.min(4096));
    for _ in 0..n {
        v.push(c.itf8()?);
    }
    Ok(v)
}

fn read_preservation(c: &mut Cursor<'_>) -> Result<Preservation> {
    let byte_len =
        usize::try_from(c.itf8()?).map_err(|_| SeqError::malformed(FORMAT, "negative map size"))?;
    let body = c.bytes(byte_len)?;
    let mut m = Cursor::new(body);
    let n = m.itf8()?;
    let mut p = Preservation {
        read_names: true,
        ap_delta: true,
        reference_required: true,
        substitution_matrix: [0; 5],
        tag_dictionary: Vec::new(),
    };
    for _ in 0..n {
        match read_key(&mut m)? {
            [b'R', b'N'] => p.read_names = read_bool(&mut m)?,
            [b'A', b'P'] => p.ap_delta = read_bool(&mut m)?,
            [b'R', b'R'] => p.reference_required = read_bool(&mut m)?,
            [b'S', b'M'] => {
                let sm = m.bytes(5)?;
                p.substitution_matrix.copy_from_slice(sm);
            }
            [b'T', b'D'] => {
                let td_len = usize::try_from(m.itf8()?)
                    .map_err(|_| SeqError::malformed(FORMAT, "negative TD size"))?;
                let td = m.bytes(td_len)?;
                p.tag_dictionary = parse_tag_dictionary(td)?;
            }
            other => {
                return Err(SeqError::Unsupported(format!(
                    "preservation map key {}{}",
                    other[0] as char, other[1] as char
                )));
            }
        }
    }
    Ok(p)
}

/// The tag dictionary is NUL-terminated lists of 3-byte tag ids (two name bytes
/// and a type byte), one list per record layout.
fn parse_tag_dictionary(td: &[u8]) -> Result<Vec<Vec<[u8; 3]>>> {
    let mut out = Vec::new();
    if td.is_empty() {
        return Ok(out);
    }
    // Each entry is NUL-terminated, so drop the final terminator before splitting.
    let body = td.strip_suffix(&[0]).unwrap_or(td);
    for chunk in body.split(|&b| b == 0) {
        if chunk.len() % 3 != 0 {
            return Err(SeqError::malformed(
                FORMAT,
                "tag dictionary entry not a multiple of 3 bytes",
            ));
        }
        out.push(chunk.chunks_exact(3).map(|t| [t[0], t[1], t[2]]).collect());
    }
    Ok(out)
}

fn read_encoding_map(c: &mut Cursor<'_>) -> Result<Vec<([u8; 2], Encoding)>> {
    let byte_len =
        usize::try_from(c.itf8()?).map_err(|_| SeqError::malformed(FORMAT, "negative map size"))?;
    let body = c.bytes(byte_len)?;
    let mut m = Cursor::new(body);
    let n = m.itf8()?;
    let mut out = Vec::new();
    for _ in 0..n {
        let key = read_key(&mut m)?;
        let enc = read_encoding(&mut m)?;
        out.push((key, enc));
    }
    Ok(out)
}

/// Parses a compression header block payload.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for truncated maps and [`SeqError::Unsupported`]
/// for unknown preservation keys or encoding ids.
pub fn parse(data: &[u8]) -> Result<CompressionHeader> {
    let mut c = Cursor::new(data);
    let preservation = read_preservation(&mut c)?;
    let data_series = read_encoding_map(&mut c)?;
    // Tag encodings use a three-byte key (name and type) with the same map layout.
    let byte_len =
        usize::try_from(c.itf8()?).map_err(|_| SeqError::malformed(FORMAT, "negative map size"))?;
    let body = c.bytes(byte_len)?;
    let mut m = Cursor::new(body);
    let n = m.itf8()?;
    let mut tags = Vec::new();
    for _ in 0..n {
        // The key is one ITF8 integer packing the tag name and type: name0 << 16 |
        // name1 << 8 | type. samtools writes it this way.
        let packed = m.itf8()? as u32;
        let key = [(packed >> 16) as u8, (packed >> 8) as u8, packed as u8];
        let enc = read_encoding(&mut m)?;
        tags.push((key, enc));
    }
    Ok(CompressionHeader { preservation, data_series, tags })
}

/// Appends one encoding descriptor (id, parameter length, parameters).
pub fn write_encoding(enc: &Encoding, out: &mut Vec<u8>) {
    let mut p = Vec::new();
    let id = match enc {
        Encoding::Null => encoding_id::NULL,
        Encoding::External { content_id } => {
            put_itf8(&mut p, *content_id);
            encoding_id::EXTERNAL
        }
        Encoding::Huffman { symbols, lengths } => {
            put_itf8(&mut p, symbols.len() as i32);
            for &s in symbols {
                put_itf8(&mut p, s);
            }
            put_itf8(&mut p, lengths.len() as i32);
            for &l in lengths {
                put_itf8(&mut p, l as i32);
            }
            encoding_id::HUFFMAN
        }
        Encoding::Beta { offset, bits } => {
            put_itf8(&mut p, *offset);
            put_itf8(&mut p, *bits as i32);
            encoding_id::BETA
        }
        Encoding::Gamma { offset } => {
            put_itf8(&mut p, *offset);
            encoding_id::GAMMA
        }
        Encoding::ByteArrayLen { len, bytes } => {
            write_encoding(len, &mut p);
            write_encoding(bytes, &mut p);
            encoding_id::BYTE_ARRAY_LEN
        }
        Encoding::ByteArrayStop { stop, content_id } => {
            p.push(*stop);
            put_itf8(&mut p, *content_id);
            encoding_id::BYTE_ARRAY_STOP
        }
        Encoding::Subexp { offset, k } => {
            put_itf8(&mut p, *offset);
            put_itf8(&mut p, *k);
            encoding_id::SUBEXP
        }
    };
    put_itf8(out, id);
    put_itf8(out, p.len() as i32);
    out.extend_from_slice(&p);
}

/// Appends a map: byte length, entry count, then the entries.
fn write_map(entries: &[(Vec<u8>, Vec<u8>)], out: &mut Vec<u8>) {
    let mut body = Vec::new();
    put_itf8(&mut body, entries.len() as i32);
    for (key, value) in entries {
        body.extend_from_slice(key);
        body.extend_from_slice(value);
    }
    put_itf8(out, body.len() as i32);
    out.extend_from_slice(&body);
}

/// Serializes a compression header. [`parse`] reads the result back.
#[must_use]
pub fn write(h: &CompressionHeader) -> Vec<u8> {
    let p = &h.preservation;
    let mut pres: Vec<(Vec<u8>, Vec<u8>)> = vec![
        (b"RN".to_vec(), vec![u8::from(p.read_names)]),
        (b"AP".to_vec(), vec![u8::from(p.ap_delta)]),
        (b"RR".to_vec(), vec![u8::from(p.reference_required)]),
        (b"SM".to_vec(), p.substitution_matrix.to_vec()),
    ];
    let mut td = Vec::new();
    for entry in &p.tag_dictionary {
        for key in entry {
            td.extend_from_slice(key);
        }
        td.push(0);
    }
    let mut td_value = Vec::new();
    put_itf8(&mut td_value, td.len() as i32);
    td_value.extend_from_slice(&td);
    pres.push((b"TD".to_vec(), td_value));

    let mut out = Vec::new();
    write_map(&pres, &mut out);

    let series: Vec<(Vec<u8>, Vec<u8>)> = h
        .data_series
        .iter()
        .map(|(k, enc)| {
            let mut v = Vec::new();
            write_encoding(enc, &mut v);
            (k.to_vec(), v)
        })
        .collect();
    write_map(&series, &mut out);

    let tags: Vec<(Vec<u8>, Vec<u8>)> = h
        .tags
        .iter()
        .map(|(k, enc)| {
            let packed = (i32::from(k[0]) << 16) | (i32::from(k[1]) << 8) | i32::from(k[2]);
            let mut key = Vec::new();
            put_itf8(&mut key, packed);
            let mut v = Vec::new();
            write_encoding(enc, &mut v);
            (key, v)
        })
        .collect();
    write_map(&tags, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compression_header_write_parse_roundtrip() {
        let h = CompressionHeader {
            preservation: Preservation {
                read_names: true,
                ap_delta: false,
                reference_required: true,
                substitution_matrix: [0x1B; 5],
                tag_dictionary: vec![vec![], vec![*b"NMC", *b"MDZ"]],
            },
            data_series: vec![
                (*b"BF", Encoding::External { content_id: 1 }),
                (*b"RN", Encoding::ByteArrayStop { stop: 0, content_id: 2 }),
                (
                    *b"QS",
                    Encoding::ByteArrayLen {
                        len: Box::new(Encoding::Beta { offset: 0, bits: 8 }),
                        bytes: Box::new(Encoding::Gamma { offset: 1 }),
                    },
                ),
            ],
            tags: vec![(*b"NMC", Encoding::Subexp { offset: 0, k: 3 })],
        };
        let bytes = write(&h);
        assert_eq!(parse(&bytes).unwrap(), h);
    }

    fn encoding_bytes(id: i32, params: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        put_itf8(&mut out, id);
        put_itf8(&mut out, params.len() as i32);
        out.extend_from_slice(params);
        out
    }

    #[test]
    fn external_and_byte_array_stop_encodings() {
        let mut ext_params = Vec::new();
        put_itf8(&mut ext_params, 7);
        let bytes = encoding_bytes(encoding_id::EXTERNAL, &ext_params);
        assert_eq!(
            read_encoding(&mut Cursor::new(&bytes)).unwrap(),
            Encoding::External { content_id: 7 }
        );

        let bytes = encoding_bytes(encoding_id::BYTE_ARRAY_STOP, &[0, 9]);
        assert_eq!(
            read_encoding(&mut Cursor::new(&bytes)).unwrap(),
            Encoding::ByteArrayStop { stop: 0, content_id: 9 }
        );
    }

    #[test]
    fn huffman_constant_encoding() {
        // Single-symbol Huffman code: symbol 0x41 with length 0.
        let mut params = Vec::new();
        put_itf8(&mut params, 1);
        put_itf8(&mut params, 0x41);
        put_itf8(&mut params, 1);
        put_itf8(&mut params, 0);
        let bytes = encoding_bytes(encoding_id::HUFFMAN, &params);
        assert_eq!(
            read_encoding(&mut Cursor::new(&bytes)).unwrap(),
            Encoding::Huffman { symbols: vec![0x41], lengths: vec![0] }
        );
    }

    #[test]
    fn unknown_encoding_is_unsupported() {
        let bytes = encoding_bytes(99, &[]);
        assert!(matches!(read_encoding(&mut Cursor::new(&bytes)), Err(SeqError::Unsupported(_))));
    }

    #[test]
    fn tag_dictionary_split() {
        let td = b"NMC\0MDZ\0";
        let parsed = parse_tag_dictionary(td).unwrap();
        assert_eq!(parsed, vec![vec![*b"NMC"], vec![*b"MDZ"]]);
    }
}
