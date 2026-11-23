//! Low-level CRAM primitives: ITF8/LTF8 variable-length integers and CRC32.
//!
//! ITF8 and LTF8 are the CRAM 3.0 variable-length integer encodings: the number
//! of leading one bits in the first byte gives the count of extra bytes.

use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM";

/// CRC32 (IEEE, reflected, polynomial 0xEDB88320) as used by CRAM block and
/// container checksums. Implemented with a 256-entry table.
#[must_use]
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &b in data {
        crc = CRC_TABLE[((crc ^ u32::from(b)) & 0xFF) as usize] ^ (crc >> 8);
    }
    !crc
}

const CRC_TABLE: [u32; 256] = {
    let mut table = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut c = i as u32;
        let mut k = 0;
        while k < 8 {
            c = if c & 1 != 0 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
            k += 1;
        }
        table[i] = c;
        i += 1;
    }
    table
};

/// Appends an ITF8-encoded value. Values are treated as unsigned 32-bit.
pub fn put_itf8(out: &mut Vec<u8>, v: i32) {
    let u = v as u32;
    if u < 0x80 {
        out.push(u as u8);
    } else if u < 0x4000 {
        out.push(0x80 | (u >> 8) as u8);
        out.push(u as u8);
    } else if u < 0x20_0000 {
        out.push(0xC0 | (u >> 16) as u8);
        out.push((u >> 8) as u8);
        out.push(u as u8);
    } else if u < 0x1000_0000 {
        out.push(0xE0 | (u >> 24) as u8);
        out.push((u >> 16) as u8);
        out.push((u >> 8) as u8);
        out.push(u as u8);
    } else {
        out.push(0xF0 | ((u >> 28) as u8 & 0x0F));
        out.push((u >> 20) as u8);
        out.push((u >> 12) as u8);
        out.push((u >> 4) as u8);
        out.push((u & 0x0F) as u8);
    }
}

/// Appends an LTF8-encoded 64-bit value.
pub fn put_ltf8(out: &mut Vec<u8>, v: i64) {
    let u = v as u64;
    let n = if u < 1 << 7 {
        1
    } else if u < 1 << 14 {
        2
    } else if u < 1 << 21 {
        3
    } else if u < 1 << 28 {
        4
    } else if u < 1 << 35 {
        5
    } else if u < 1 << 42 {
        6
    } else if u < 1 << 49 {
        7
    } else if u < 1 << 56 {
        8
    } else {
        9
    };
    if n == 9 {
        out.push(0xFF);
        out.extend_from_slice(&u.to_be_bytes());
        return;
    }
    // Prefix: n-1 leading one bits followed by a zero bit.
    let prefix: u8 = !(0xFFu8 >> (n - 1));
    let high_bits = (u >> (8 * (n - 1))) as u8;
    out.push(prefix | high_bits);
    for i in (0..n - 1).rev() {
        out.push((u >> (8 * i)) as u8);
    }
}

/// A cursor over CRAM bytes that decodes integers and fixed-width fields.
#[derive(Debug, Clone)]
pub struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    /// Starts reading at the beginning of `data`.
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    /// Current byte offset.
    #[must_use]
    pub fn position(&self) -> usize {
        self.pos
    }

    /// Bytes not yet consumed.
    #[must_use]
    pub fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    /// Whether all bytes have been consumed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.remaining() == 0
    }

    fn truncated() -> SeqError {
        SeqError::malformed(FORMAT, "truncated data")
    }

    /// The bytes between two absolute offsets, for checksum computation.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when the range is out of bounds.
    pub fn slice_between(&self, start: usize, end: usize) -> Result<&'a [u8]> {
        self.data.get(start..end).ok_or_else(Self::truncated)
    }

    /// Reads one byte.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when no bytes remain.
    pub fn u8(&mut self) -> Result<u8> {
        let b = *self.data.get(self.pos).ok_or_else(Self::truncated)?;
        self.pos += 1;
        Ok(b)
    }

    /// Reads `n` bytes as a slice.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when fewer than `n` bytes remain.
    pub fn bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).ok_or_else(Self::truncated)?;
        let s = self.data.get(self.pos..end).ok_or_else(Self::truncated)?;
        self.pos = end;
        Ok(s)
    }

    /// Reads a little-endian `i32`.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on truncation.
    pub fn i32_le(&mut self) -> Result<i32> {
        let b = self.bytes(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// Reads a little-endian `u32`.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on truncation.
    pub fn u32_le(&mut self) -> Result<u32> {
        let b = self.bytes(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// Reads an ITF8 integer, returned as the raw 32-bit pattern.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on truncation.
    pub fn itf8(&mut self) -> Result<i32> {
        let b0 = u32::from(self.u8()?);
        let v = if b0 & 0x80 == 0 {
            b0
        } else if b0 & 0x40 == 0 {
            ((b0 & 0x3F) << 8) | u32::from(self.u8()?)
        } else if b0 & 0x20 == 0 {
            let b1 = u32::from(self.u8()?);
            let b2 = u32::from(self.u8()?);
            ((b0 & 0x1F) << 16) | (b1 << 8) | b2
        } else if b0 & 0x10 == 0 {
            let b1 = u32::from(self.u8()?);
            let b2 = u32::from(self.u8()?);
            let b3 = u32::from(self.u8()?);
            ((b0 & 0x0F) << 24) | (b1 << 16) | (b2 << 8) | b3
        } else {
            let b1 = u32::from(self.u8()?);
            let b2 = u32::from(self.u8()?);
            let b3 = u32::from(self.u8()?);
            let b4 = u32::from(self.u8()?);
            ((b0 & 0x0F) << 28) | (b1 << 20) | (b2 << 12) | (b3 << 4) | (b4 & 0x0F)
        };
        Ok(v as i32)
    }

    /// Reads an LTF8 integer, returned as the raw 64-bit pattern.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on truncation.
    pub fn ltf8(&mut self) -> Result<i64> {
        let b0 = self.u8()?;
        let extra = match b0 {
            x if x & 0x80 == 0 => 0,
            x if x & 0x40 == 0 => 1,
            x if x & 0x20 == 0 => 2,
            x if x & 0x10 == 0 => 3,
            x if x & 0x08 == 0 => 4,
            x if x & 0x04 == 0 => 5,
            x if x & 0x02 == 0 => 6,
            x if x & 0x01 == 0 => 7,
            _ => 8,
        };
        let mut v: u64 = if extra == 8 { 0 } else { u64::from(b0) & (0xFFu64 >> (extra + 1)) };
        for _ in 0..extra {
            v = (v << 8) | u64::from(self.u8()?);
        }
        Ok(v as i64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc32_check_value() {
        // Standard CRC-32 check value for the ASCII digits "123456789".
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc32(b""), 0);
    }

    #[test]
    fn itf8_known_encodings() {
        // Byte patterns taken directly from the CRAM 3.0 specification table.
        let cases: &[(i32, &[u8])] = &[
            (0, &[0x00]),
            (127, &[0x7F]),
            (128, &[0x80, 0x80]),
            (0x3FFF, &[0xBF, 0xFF]),
            (0x4000, &[0xC0, 0x40, 0x00]),
            (0x1F_FFFF, &[0xDF, 0xFF, 0xFF]),
            (0x20_0000, &[0xE0, 0x20, 0x00, 0x00]),
            (0x0FFF_FFFF, &[0xEF, 0xFF, 0xFF, 0xFF]),
            (0x1000_0000, &[0xF1, 0x00, 0x00, 0x00, 0x00]),
        ];
        for &(v, bytes) in cases {
            let mut out = Vec::new();
            put_itf8(&mut out, v);
            assert_eq!(out, bytes, "encode {v}");
            let mut c = Cursor::new(bytes);
            assert_eq!(c.itf8().unwrap(), v, "decode {v}");
            assert!(c.is_empty());
        }
    }

    #[test]
    fn ltf8_roundtrip_boundaries() {
        for &v in &[0i64, 1, 127, 128, 16_383, 16_384, 1 << 35, (1 << 56) - 1, 1 << 56, -1] {
            let mut out = Vec::new();
            put_ltf8(&mut out, v);
            let mut c = Cursor::new(&out);
            assert_eq!(c.ltf8().unwrap(), v, "value {v}");
            assert!(c.is_empty(), "value {v} left bytes");
        }
    }
}
