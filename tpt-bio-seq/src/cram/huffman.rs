//! Canonical Huffman decoding for CRAM HUFFMAN encodings, and an MSB-first bit reader.
//!
//! Codes are assigned canonically: symbols sorted by (code length, symbol value)
//! receive consecutive codes, with a zero-length code meaning a single-symbol
//! alphabet that consumes no bits. Bits are read most-significant first within
//! each byte, as in the CRAM core data stream.

use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM Huffman";
/// Longest code length this decoder accepts.
const MAX_LEN: usize = 32;

/// MSB-first bit reader over a byte slice.
#[derive(Debug, Clone)]
pub struct BitReader<'a> {
    data: &'a [u8],
    bit: usize,
}

impl<'a> BitReader<'a> {
    /// Starts at the first bit of `data`.
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, bit: 0 }
    }

    /// Starts at bit `bit` (counted from the first bit of `data`).
    #[must_use]
    pub fn at(data: &'a [u8], bit: usize) -> Self {
        Self { data, bit }
    }

    /// Index of the next bit to read.
    #[must_use]
    pub fn position(&self) -> usize {
        self.bit
    }

    /// Reads one bit.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when the data is exhausted.
    pub fn bit(&mut self) -> Result<u32> {
        let byte = *self
            .data
            .get(self.bit / 8)
            .ok_or_else(|| SeqError::malformed(FORMAT, "bit stream exhausted"))?;
        let b = (byte >> (7 - (self.bit % 8))) & 1;
        self.bit += 1;
        Ok(u32::from(b))
    }

    /// Reads `n` bits as an unsigned value, most significant first.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] when the data is exhausted.
    pub fn bits(&mut self, n: u32) -> Result<u32> {
        let mut v = 0u32;
        for _ in 0..n {
            v = (v << 1) | self.bit()?;
        }
        Ok(v)
    }
}

/// A canonical Huffman decoder built from symbol and code-length lists.
#[derive(Debug, Clone)]
pub struct Huffman {
    /// Symbols sorted by (length, symbol).
    symbols: Vec<i32>,
    /// For each code length `l`, the first canonical code of that length.
    first_code: [u32; MAX_LEN + 1],
    /// For each code length `l`, the index into `symbols` of its first symbol.
    first_index: [usize; MAX_LEN + 1],
    /// Number of codes of each length.
    count: [usize; MAX_LEN + 1],
    /// Length zero: the alphabet has one symbol and no bits are read.
    single: Option<i32>,
}

impl Huffman {
    /// Builds the decoder. `symbols[i]` has code length `lengths[i]`.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for mismatched lists, lengths above 32, or
    /// a code set that over-subscribes the code space.
    pub fn new(symbols: &[i32], lengths: &[u32]) -> Result<Self> {
        if symbols.len() != lengths.len() || symbols.is_empty() {
            return Err(SeqError::malformed(
                FORMAT,
                "symbol and length lists must match and be non-empty",
            ));
        }
        if lengths.iter().any(|&l| l as usize > MAX_LEN) {
            return Err(SeqError::malformed(FORMAT, "code length above 32"));
        }
        if lengths.iter().all(|&l| l == 0) {
            if symbols.len() != 1 {
                return Err(SeqError::malformed(FORMAT, "zero-length codes need one symbol"));
            }
            return Ok(Self {
                symbols: symbols.to_vec(),
                first_code: [0; MAX_LEN + 1],
                first_index: [0; MAX_LEN + 1],
                count: [0; MAX_LEN + 1],
                single: Some(symbols[0]),
            });
        }
        let mut pairs: Vec<(u32, i32)> =
            lengths.iter().copied().zip(symbols.iter().copied()).collect();
        pairs.sort_unstable();
        let mut count = [0usize; MAX_LEN + 1];
        for &(l, _) in &pairs {
            count[l as usize] += 1;
        }
        count[0] = 0;
        // Kraft check: the code space must not be over-subscribed.
        let mut space: u64 = 0;
        for (l, &n) in count.iter().enumerate().skip(1) {
            space += (n as u64) << (MAX_LEN - l);
        }
        if space > 1u64 << MAX_LEN {
            return Err(SeqError::malformed(FORMAT, "code lengths over-subscribe the code space"));
        }
        let mut first_code = [0u32; MAX_LEN + 1];
        let mut first_index = [0usize; MAX_LEN + 1];
        let mut code = 0u32;
        let mut index = 0usize;
        for l in 1..=MAX_LEN {
            code = (code + count[l - 1] as u32) << 1;
            first_code[l] = code;
            first_index[l] = index;
            index += count[l];
        }
        Ok(Self {
            symbols: pairs.into_iter().map(|(_, s)| s).collect(),
            first_code,
            first_index,
            count,
            single: None,
        })
    }

    /// Decodes one symbol from the bit reader.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a bit pattern that is not a valid code,
    /// or when the bit stream ends early.
    pub fn decode(&self, bits: &mut BitReader<'_>) -> Result<i32> {
        if let Some(s) = self.single {
            return Ok(s);
        }
        let mut code = 0u32;
        for l in 1..=MAX_LEN {
            code = (code << 1) | bits.bit()?;
            let n = self.count[l];
            if n > 0 && code >= self.first_code[l] && code < self.first_code[l] + n as u32 {
                let idx = self.first_index[l] + (code - self.first_code[l]) as usize;
                return Ok(self.symbols[idx]);
            }
        }
        Err(SeqError::malformed(FORMAT, "invalid Huffman code"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_codes_from_lengths() {
        // Lengths 1,2,3,3 give A=0, B=10, C=110, D=111.
        let h = Huffman::new(&[65, 66, 67, 68], &[1, 2, 3, 3]).unwrap();
        // Stream: 0 | 10 | 110 | 111 -> 0101 1011 1000 0000
        let data = [0x5B, 0x80];
        let mut r = BitReader::new(&data);
        let got: Vec<i32> = (0..4).map(|_| h.decode(&mut r).unwrap()).collect();
        assert_eq!(got, vec![65, 66, 67, 68]);
    }

    #[test]
    fn unsorted_symbols_get_canonical_order() {
        // Symbols listed out of order: codes still follow (length, symbol) order.
        let h = Huffman::new(&[68, 65], &[1, 1]).unwrap();
        let data = [0b0100_0000];
        let mut r = BitReader::new(&data);
        assert_eq!(h.decode(&mut r).unwrap(), 65);
        assert_eq!(h.decode(&mut r).unwrap(), 68);
    }

    #[test]
    fn single_symbol_reads_no_bits() {
        let h = Huffman::new(&[42], &[0]).unwrap();
        let mut r = BitReader::new(&[]);
        assert_eq!(h.decode(&mut r).unwrap(), 42);
    }

    #[test]
    fn oversubscribed_code_is_rejected() {
        assert!(Huffman::new(&[1, 2, 3], &[1, 1, 1]).is_err());
    }
}
