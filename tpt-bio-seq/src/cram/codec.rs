//! Value decoding for CRAM data series.
//!
//! Each data series names one encoding. Integer and byte values come either from
//! an external block (byte-aligned, read in order, shared by every series that
//! names the same content id) or from the core bit stream (Huffman, Beta, Gamma
//! and Subexp codes, read MSB first). Byte arrays are a length encoding followed
//! by per-byte encodings, or a run of external bytes ended by a stop byte.
//!
//! The Gamma and Subexp bit layouts follow the usual Elias-gamma and
//! sub-exponential definitions; the CRAM specification text was not available
//! to check them against, so they are covered only by round-trip tests.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// A hash map with a cheap hasher. CRAM looks up small fixed keys (series, tags and
/// external block ids) once per value, so `SipHash` dominates decoding time.
pub(crate) type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<KeyHasher>>;

/// A multiplicative hasher for small integer and byte keys. It is not DoS-resistant, which is
/// fine here: the keys come from the CRAM header, and their table is bounded.
#[derive(Default)]
pub(crate) struct KeyHasher(u64);

impl Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ u64::from(b)).wrapping_mul(0x0100_0000_01B3);
        }
    }

    fn write_i32(&mut self, i: i32) {
        self.0 = u64::from(i as u32).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

use crate::cram::bytes::Cursor;
use crate::cram::compression::Encoding;
use crate::cram::huffman::{BitReader, Huffman};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM codec";

/// A ready-to-use value decoder built from an [`Encoding`].
#[derive(Debug, Clone)]
pub enum Codec {
    /// No data: integers read as zero, byte arrays as empty.
    Null,
    /// ITF8 integers or single bytes from an external block.
    External(i32),
    /// Canonical Huffman symbols from the core stream.
    Huffman(Box<Huffman>),
    /// `bits`-wide fixed-width value from the core stream, minus `offset`.
    Beta {
        /// Subtracted from the raw value.
        offset: i32,
        /// Width in bits.
        bits: u32,
    },
    /// Elias-gamma value from the core stream, minus `offset`.
    Gamma {
        /// Subtracted from the raw value.
        offset: i32,
    },
    /// Sub-exponential value from the core stream, minus `offset`.
    Subexp {
        /// Subtracted from the raw value.
        offset: i32,
        /// The `k` parameter.
        k: u32,
    },
    /// Length-prefixed byte array.
    ByteArrayLen {
        /// Decodes the array length.
        len: Box<Codec>,
        /// Decodes each byte.
        bytes: Box<Codec>,
    },
    /// Bytes from an external block up to and excluding `stop`.
    ByteArrayStop {
        /// Terminating byte value.
        stop: u8,
        /// External block id.
        content_id: i32,
    },
}

impl Codec {
    /// Builds a decoder for `enc`.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for parameters out of range.
    pub fn from_encoding(enc: &Encoding) -> Result<Self> {
        Ok(match enc {
            Encoding::Null => Self::Null,
            Encoding::External { content_id } => Self::External(*content_id),
            Encoding::Huffman { symbols, lengths } => {
                Self::Huffman(Box::new(Huffman::new(symbols, lengths)?))
            }
            Encoding::Beta { offset, bits } => {
                if *bits > 32 {
                    return Err(SeqError::malformed(FORMAT, "BETA width above 32 bits"));
                }
                Self::Beta { offset: *offset, bits: *bits }
            }
            Encoding::Gamma { offset } => Self::Gamma { offset: *offset },
            Encoding::Subexp { offset, k } => {
                let k = u32::try_from(*k)
                    .ok()
                    .filter(|&k| k <= 31)
                    .ok_or_else(|| SeqError::malformed(FORMAT, "SUBEXP k out of range"))?;
                Self::Subexp { offset: *offset, k }
            }
            Encoding::ByteArrayLen { len, bytes } => Self::ByteArrayLen {
                len: Box::new(Self::from_encoding(len)?),
                bytes: Box::new(Self::from_encoding(bytes)?),
            },
            Encoding::ByteArrayStop { stop, content_id } => {
                Self::ByteArrayStop { stop: *stop, content_id: *content_id }
            }
        })
    }
}

/// Content ids below this bound are looked up in a table rather than scanned.
const DENSE_IDS: usize = 1024;

/// An external block with its read position.
#[derive(Debug)]
struct External {
    data: Vec<u8>,
    pos: usize,
}

/// Reads values for data series from a slice's core block and external blocks.
#[derive(Debug)]
pub struct DataReader {
    core: Vec<u8>,
    core_bit: usize,
    /// External blocks, with their content ids in `ids` at the same index.
    ids: Vec<i32>,
    externals: Vec<External>,
    /// For small non-negative content ids, the index into `externals` plus one (zero if none).
    by_id: Vec<u32>,
}

impl DataReader {
    /// Creates a reader over `core` and `(content id, bytes)` external blocks.
    #[must_use]
    pub fn new(core: Vec<u8>, externals: Vec<(i32, Vec<u8>)>) -> Self {
        let (ids, externals): (Vec<i32>, Vec<External>) =
            externals.into_iter().map(|(id, data)| (id, External { data, pos: 0 })).unzip();
        let mut by_id = vec![0u32; DENSE_IDS];
        for (i, &id) in ids.iter().enumerate() {
            // The first block with an id wins, as a scan would find it.
            if let Some(slot) = usize::try_from(id).ok().and_then(|id| by_id.get_mut(id))
                && *slot == 0
            {
                *slot = i as u32 + 1;
            }
        }
        Self { core, core_bit: 0, ids, externals, by_id }
    }

    fn with_bits<T>(&mut self, f: impl FnOnce(&mut BitReader<'_>) -> Result<T>) -> Result<T> {
        let mut bits = BitReader::at(&self.core, self.core_bit);
        let out = f(&mut bits);
        self.core_bit = bits.position();
        out
    }

    fn external(&mut self, id: i32) -> Result<&mut External> {
        // A slice has a handful of external blocks, so a scan beats hashing per value.
        let dense = usize::try_from(id).ok().and_then(|id| self.by_id.get(id)).copied();
        let found = match dense {
            Some(slot) => slot.checked_sub(1).map(|i| i as usize),
            None => self.ids.iter().position(|&known| known == id),
        };
        match found {
            Some(i) => Ok(&mut self.externals[i]),
            None => Err(SeqError::malformed(FORMAT, format!("missing external block {id}"))),
        }
    }

    fn external_itf8(&mut self, id: i32) -> Result<i32> {
        let ext = self.external(id)?;
        let mut c = Cursor::new(&ext.data[ext.pos.min(ext.data.len())..]);
        let v = c.itf8()?;
        ext.pos += c.position();
        Ok(v)
    }

    fn external_byte(&mut self, id: i32) -> Result<u8> {
        let ext = self.external(id)?;
        let b = *ext
            .data
            .get(ext.pos)
            .ok_or_else(|| SeqError::malformed(FORMAT, format!("external block {id} exhausted")))?;
        ext.pos += 1;
        Ok(b)
    }

    /// Decodes one integer value.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when the data runs out or a code is invalid,
    /// and for byte-array codecs used where an integer is expected.
    pub fn int(&mut self, codec: &Codec) -> Result<i32> {
        match codec {
            Codec::Null => Ok(0),
            Codec::External(id) => self.external_itf8(*id),
            Codec::Huffman(h) => self.with_bits(|b| h.decode(b)),
            Codec::Beta { offset, bits } => {
                let raw = self.with_bits(|b| b.bits(*bits))?;
                to_i32(i64::from(raw) - i64::from(*offset))
            }
            Codec::Gamma { offset } => {
                let raw = self.with_bits(read_gamma)?;
                to_i32(raw - i64::from(*offset))
            }
            Codec::Subexp { offset, k } => {
                let raw = self.with_bits(|b| read_subexp(b, *k))?;
                to_i32(raw - i64::from(*offset))
            }
            Codec::ByteArrayLen { .. } | Codec::ByteArrayStop { .. } => {
                Err(SeqError::malformed(FORMAT, "byte array encoding used for an integer"))
            }
        }
    }

    /// Decodes one byte value.
    ///
    /// # Errors
    /// As [`DataReader::int`], and when the decoded value does not fit in a byte.
    pub fn byte(&mut self, codec: &Codec) -> Result<u8> {
        if let Codec::External(id) = codec {
            return self.external_byte(*id);
        }
        let v = self.int(codec)?;
        u8::try_from(v).map_err(|_| SeqError::malformed(FORMAT, format!("byte value {v}")))
    }

    /// Decodes `n` byte values with one codec. An external codec resolves its block once,
    /// not once per byte.
    ///
    /// # Errors
    /// As [`DataReader::byte`].
    pub(crate) fn byte_run(&mut self, codec: &Codec, n: usize) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.byte_run_into(codec, n, &mut out)?;
        Ok(out)
    }

    /// Appends `n` byte values decoded with one codec to `out`.
    fn byte_run_into(&mut self, codec: &Codec, n: usize, out: &mut Vec<u8>) -> Result<()> {
        if let Codec::External(id) = codec {
            return self.external_run_into(*id, n, out);
        }
        out.reserve(n.min(1 << 16));
        for _ in 0..n {
            out.push(self.byte(codec)?);
        }
        Ok(())
    }

    /// Appends `n` bytes from an external block to `out`.
    fn external_run_into(&mut self, id: i32, n: usize, out: &mut Vec<u8>) -> Result<()> {
        let ext = self.external(id)?;
        let rest = &ext.data[ext.pos.min(ext.data.len())..];
        if rest.len() < n {
            return Err(SeqError::malformed(FORMAT, format!("external block {id} exhausted")));
        }
        out.extend_from_slice(&rest[..n]);
        ext.pos += n;
        Ok(())
    }

    /// Decodes one byte array.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a negative length or exhausted data.
    pub fn byte_array(&mut self, codec: &Codec) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.byte_array_into(codec, &mut out)?;
        Ok(out)
    }

    /// Decodes one byte array and appends it to `out`. Returns the number of bytes appended.
    ///
    /// # Errors
    /// As [`DataReader::byte_array`].
    pub(crate) fn byte_array_into(&mut self, codec: &Codec, out: &mut Vec<u8>) -> Result<usize> {
        match codec {
            Codec::Null => Ok(0),
            Codec::ByteArrayLen { len, bytes } => {
                let n = usize::try_from(self.int(len)?)
                    .map_err(|_| SeqError::malformed(FORMAT, "negative array length"))?;
                self.byte_run_into(bytes, n, out)?;
                Ok(n)
            }
            Codec::ByteArrayStop { stop, content_id } => {
                let stop = *stop;
                let ext = self.external(*content_id)?;
                let rest = &ext.data[ext.pos.min(ext.data.len())..];
                let end = rest.iter().position(|&b| b == stop).ok_or_else(|| {
                    SeqError::malformed(FORMAT, format!("external block {content_id} exhausted"))
                })?;
                out.extend_from_slice(&rest[..end]);
                ext.pos += end + 1;
                Ok(end)
            }
            _ => Err(SeqError::malformed(FORMAT, "integer encoding used for a byte array")),
        }
    }
}

fn to_i32(v: i64) -> Result<i32> {
    i32::try_from(v).map_err(|_| SeqError::malformed(FORMAT, format!("value {v} out of range")))
}

/// Elias-gamma code: `n` zero bits, a one bit, then `n` bits. Value is `2^n + bits`.
fn read_gamma(b: &mut BitReader<'_>) -> Result<i64> {
    let mut n = 0u32;
    while b.bit()? == 0 {
        n += 1;
        if n > 31 {
            return Err(SeqError::malformed(FORMAT, "gamma prefix too long"));
        }
    }
    let rest = b.bits(n)?;
    Ok((1i64 << n) + i64::from(rest))
}

/// Sub-exponential code with parameter `k`: a run of `n` one bits then a zero.
/// `n == 0` reads `k` bits; otherwise reads `k + n - 1` bits above `2^(k+n-1)`.
fn read_subexp(b: &mut BitReader<'_>, k: u32) -> Result<i64> {
    let mut n = 0u32;
    while b.bit()? == 1 {
        n += 1;
        if n > 32 {
            return Err(SeqError::malformed(FORMAT, "sub-exponential prefix too long"));
        }
    }
    if n == 0 {
        return Ok(i64::from(b.bits(k)?));
    }
    let width = k + n - 1;
    if width > 32 {
        return Err(SeqError::malformed(FORMAT, "sub-exponential width too large"));
    }
    Ok((1i64 << width) + i64::from(b.bits(width)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gamma_codec(offset: i32) -> Codec {
        Codec::Gamma { offset }
    }

    #[test]
    fn beta_reads_fixed_width_minus_offset() {
        // 4-bit values: 0b1010 = 10, 0b0001 = 1; offset 1 gives 9 and 0.
        let mut r = DataReader::new(vec![0b1010_0001], vec![]);
        let c = Codec::Beta { offset: 1, bits: 4 };
        assert_eq!(r.int(&c).unwrap(), 9);
        assert_eq!(r.int(&c).unwrap(), 0);
    }

    #[test]
    fn gamma_hand_derived_bits() {
        // Elias gamma: 1 -> "1", 2 -> "010", 3 -> "011", 4 -> "00100".
        // Stream "1 010 011 00100" = 1010 0110 0100 0000 -> 0xA6 0x40.
        let mut r = DataReader::new(vec![0xA6, 0x40], vec![]);
        let c = gamma_codec(0);
        let got: Vec<i32> = (0..4).map(|_| r.int(&c).unwrap()).collect();
        assert_eq!(got, vec![1, 2, 3, 4]);
    }

    #[test]
    fn gamma_offset_is_subtracted() {
        // Gamma of 1 with offset 1 is value 0.
        let mut r = DataReader::new(vec![0b1000_0000], vec![]);
        assert_eq!(r.int(&gamma_codec(1)).unwrap(), 0);
    }

    #[test]
    fn subexp_ranges_are_contiguous() {
        // k = 2: codes "0xx" give 0..3, "10xx" gives 4..7, "110xxx" gives 8..15.
        let mut r = DataReader::new(vec![0b0111_0000, 0b1110_0000], vec![]);
        let c = Codec::Subexp { offset: 0, k: 2 };
        assert_eq!(r.int(&c).unwrap(), 3); // "0" then "11"
        assert_eq!(r.int(&c).unwrap(), 4); // "10" then "00" (n=1, width 2, base 4)
    }

    #[test]
    fn external_itf8_and_bytes_share_a_block() {
        let mut r = DataReader::new(vec![], vec![(3, vec![0x05, 0x81, 0x02, 0x7A])]);
        let ints = Codec::External(3);
        assert_eq!(r.int(&ints).unwrap(), 5);
        assert_eq!(r.int(&ints).unwrap(), 0x102);
        assert_eq!(r.byte(&Codec::External(3)).unwrap(), 0x7A);
        assert!(r.byte(&Codec::External(3)).is_err());
    }

    #[test]
    fn byte_array_stop_and_length() {
        let mut r = DataReader::new(vec![], vec![(1, b"read1\0AB".to_vec())]);
        let stop = Codec::ByteArrayStop { stop: 0, content_id: 1 };
        assert_eq!(r.byte_array(&stop).unwrap(), b"read1");
        let len = Codec::ByteArrayLen {
            len: Box::new(Codec::External(2)),
            bytes: Box::new(Codec::External(1)),
        };
        let mut r = DataReader::new(vec![], vec![(1, b"XYZ".to_vec()), (2, vec![2])]);
        assert_eq!(r.byte_array(&len).unwrap(), b"XY");
    }

    #[test]
    fn huffman_symbols_from_core() {
        let h = Huffman::new(&[10, 20], &[1, 1]).unwrap();
        let mut r = DataReader::new(vec![0b0100_0000], vec![]);
        let c = Codec::Huffman(Box::new(h));
        assert_eq!(r.int(&c).unwrap(), 10);
        assert_eq!(r.int(&c).unwrap(), 20);
    }

    #[test]
    fn null_codec_reads_nothing() {
        let mut r = DataReader::new(vec![], vec![]);
        assert_eq!(r.int(&Codec::Null).unwrap(), 0);
        assert!(r.byte_array(&Codec::Null).unwrap().is_empty());
    }
}
