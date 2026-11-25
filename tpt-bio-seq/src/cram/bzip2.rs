//! bzip2 decoder, decode only. Used for CRAM block method 2 and for the external
//! (bzip2) form of TOK3 arithmetic streams.
//!
//! Written from the bzip2 file format: stream header `BZh1`-`BZh9`, blocks headed by
//! the magic `0x314159265359`, an end-of-stream marker `0x177245385090` carrying the
//! combined CRC, then Huffman-coded MTF and run-length symbols, the inverse
//! Burrows-Wheeler transform, and the initial run-length stage (four equal bytes are
//! followed by a repeat count). Block and stream CRCs use the MSB-first CRC-32 with
//! polynomial `0x04C11DB7`.
//!
//! Randomised blocks (a format feature retired in bzip2 0.9.5) are rejected. The
//! encoder has not produced them since then.

use crate::error::{Result, SeqError};

const FORMAT: &str = "bzip2";
const BLOCK_MAGIC: u64 = 0x3141_5926_5359;
const EOS_MAGIC: u64 = 0x1772_4538_5090;
const GROUP_SIZE: usize = 50;
const MAX_CODE_LEN: usize = 20;
const MAX_GROUPS: usize = 6;
const MIN_GROUPS: usize = 2;

const CRC_TABLE: [u32; 256] = make_crc_table();

const fn make_crc_table() -> [u32; 256] {
    let mut table = [0_u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut c = (i as u32) << 24;
        let mut bit = 0;
        while bit < 8 {
            c = if c & 0x8000_0000 != 0 { (c << 1) ^ 0x04C1_1DB7 } else { c << 1 };
            bit += 1;
        }
        table[i] = c;
        i += 1;
    }
    table
}

fn block_crc(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFF_u32;
    for &b in data {
        crc = (crc << 8) ^ CRC_TABLE[((crc >> 24) ^ u32::from(b)) as usize];
    }
    !crc
}

fn malformed(reason: impl Into<String>) -> SeqError {
    SeqError::malformed(FORMAT, reason)
}

fn truncated() -> SeqError {
    malformed("stream is truncated")
}

/// Reads bits most significant first.
struct Bits<'a> {
    data: &'a [u8],
    /// Position in bits from the start of `data`.
    pos: usize,
}

impl Bits<'_> {
    fn bit(&mut self) -> Result<u32> {
        let byte = *self.data.get(self.pos >> 3).ok_or_else(truncated)?;
        let bit = (byte >> (7 - (self.pos & 7))) & 1;
        self.pos += 1;
        Ok(u32::from(bit))
    }

    /// Reads `n` bits (at most 32) as an unsigned value.
    fn bits(&mut self, n: u32) -> Result<u32> {
        let mut v = 0_u32;
        for _ in 0..n {
            v = (v << 1) | self.bit()?;
        }
        Ok(v)
    }

    fn bits48(&mut self) -> Result<u64> {
        let hi = u64::from(self.bits(24)?);
        let lo = u64::from(self.bits(24)?);
        Ok((hi << 24) | lo)
    }

    /// Skips to the next byte boundary.
    fn align(&mut self) {
        self.pos = (self.pos + 7) & !7;
    }

    /// Whether the next whole byte is available.
    fn has_byte(&self) -> bool {
        (self.pos >> 3) < self.data.len()
    }
}

/// Canonical Huffman table for one coding group, decoded bit by bit.
struct Huffman {
    /// First code of each length.
    first: [u32; MAX_CODE_LEN + 1],
    /// Number of codes of each length.
    count: [u32; MAX_CODE_LEN + 1],
    /// Index into `symbols` of the first code of each length.
    offset: [u32; MAX_CODE_LEN + 1],
    /// Symbols ordered by code length, then by symbol value.
    symbols: Vec<u16>,
    max_len: usize,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Self {
        let mut count = [0_u32; MAX_CODE_LEN + 1];
        for &len in lengths {
            count[usize::from(len)] += 1;
        }
        let max_len = (1..=MAX_CODE_LEN).rev().find(|&l| count[l] > 0).unwrap_or(1);
        let mut first = [0_u32; MAX_CODE_LEN + 1];
        let mut offset = [0_u32; MAX_CODE_LEN + 1];
        let mut code = 0_u32;
        let mut next = 0_u32;
        for len in 1..=max_len {
            first[len] = code;
            offset[len] = next;
            code = (code + count[len]) << 1;
            next += count[len];
        }
        let mut symbols = Vec::with_capacity(lengths.len());
        for len in 1..=max_len {
            for (sym, &l) in lengths.iter().enumerate() {
                if usize::from(l) == len {
                    symbols.push(sym as u16);
                }
            }
        }
        Self { first, count, offset, symbols, max_len }
    }

    fn decode(&self, bits: &mut Bits<'_>) -> Result<u16> {
        let mut code = 0_u32;
        for len in 1..=self.max_len {
            code = (code << 1) | bits.bit()?;
            let c = self.count[len];
            if c > 0 && code >= self.first[len] && code - self.first[len] < c {
                let index = self.offset[len] + (code - self.first[len]);
                return Ok(self.symbols[index as usize]);
            }
        }
        Err(malformed("invalid Huffman code"))
    }
}

/// Decompresses a bzip2 stream (or several concatenated streams).
///
/// `max_out` bounds the output size, so corrupt input cannot allocate without limit.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for corrupt or truncated data, a CRC mismatch, or
/// output longer than `max_out`, and [`SeqError::Unsupported`] for randomised blocks.
pub(super) fn decompress(input: &[u8], max_out: usize) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut bits = Bits { data: input, pos: 0 };
    loop {
        // Stream header: "BZh" then the level digit, byte aligned.
        if bits.bits(8)? != u32::from(b'B')
            || bits.bits(8)? != u32::from(b'Z')
            || bits.bits(8)? != u32::from(b'h')
        {
            return Err(malformed("missing BZh stream header"));
        }
        let level = bits.bits(8)?;
        if !(u32::from(b'1')..=u32::from(b'9')).contains(&level) {
            return Err(malformed("invalid block size level"));
        }
        let block_max = (level - u32::from(b'0')) as usize * 100_000;

        let mut combined = 0_u32;
        loop {
            let magic = bits.bits48()?;
            let stored_crc = bits.bits(32)?;
            if magic == BLOCK_MAGIC {
                let block = decode_block(&mut bits, block_max)?;
                if block_crc(&block) != stored_crc {
                    return Err(malformed("block CRC mismatch"));
                }
                combined = combined.rotate_left(1) ^ stored_crc;
                if out.len() + block.len() > max_out {
                    return Err(malformed("output exceeds the declared size"));
                }
                out.extend_from_slice(&block);
            } else if magic == EOS_MAGIC {
                if stored_crc != combined {
                    return Err(malformed("stream CRC mismatch"));
                }
                break;
            } else {
                return Err(malformed("bad block magic"));
            }
        }

        // Another stream may follow, starting at the next byte.
        bits.align();
        if !bits.has_byte() {
            break;
        }
    }
    Ok(out)
}

/// Decodes one block after its magic and CRC. Returns the block's bytes.
fn decode_block(bits: &mut Bits<'_>, block_max: usize) -> Result<Vec<u8>> {
    if bits.bit()? != 0 {
        return Err(SeqError::Unsupported("randomised bzip2 blocks".into()));
    }
    let orig_ptr = bits.bits(24)? as usize;
    let seq_to_unseq = read_symbol_map(bits)?;
    let alpha = seq_to_unseq.len() + 2;

    let n_groups = bits.bits(3)? as usize;
    if !(MIN_GROUPS..=MAX_GROUPS).contains(&n_groups) {
        return Err(malformed(format!("{n_groups} Huffman groups is out of range")));
    }
    let n_selectors = bits.bits(15)? as usize;
    if n_selectors == 0 {
        return Err(malformed("block has no selectors"));
    }
    let selectors = read_selectors(bits, n_selectors, n_groups)?;
    let tables = read_tables(bits, n_groups, alpha)?;
    let tt = decode_symbols(bits, &tables, &selectors, &seq_to_unseq, block_max)?;

    if tt.is_empty() {
        return Ok(Vec::new());
    }
    if orig_ptr >= tt.len() {
        return Err(malformed("original pointer is outside the block"));
    }
    Ok(undo_initial_rle(&inverse_bwt(&tt, orig_ptr)))
}

/// Reads the bitmap of byte values that occur in the block, returning them in
/// ascending order.
fn read_symbol_map(bits: &mut Bits<'_>) -> Result<Vec<u8>> {
    let used_groups = bits.bits(16)?;
    let mut seq_to_unseq: Vec<u8> = Vec::with_capacity(256);
    for i in 0..16_u32 {
        if used_groups & (0x8000 >> i) != 0 {
            let used = bits.bits(16)?;
            for j in 0..16_u32 {
                if used & (0x8000 >> j) != 0 {
                    seq_to_unseq.push((i * 16 + j) as u8);
                }
            }
        }
    }
    if seq_to_unseq.is_empty() {
        return Err(malformed("block uses no byte values"));
    }
    Ok(seq_to_unseq)
}

/// Reads the selectors, which name the coding group for each run of 50 symbols.
/// They are MTF-coded in unary.
fn read_selectors(bits: &mut Bits<'_>, n_selectors: usize, n_groups: usize) -> Result<Vec<u8>> {
    let mut group_mtf: [u8; MAX_GROUPS] = [0, 1, 2, 3, 4, 5];
    let mut selectors = Vec::with_capacity(n_selectors);
    for _ in 0..n_selectors {
        let mut j = 0_usize;
        while bits.bit()? == 1 {
            j += 1;
            if j >= n_groups {
                return Err(malformed("selector index is out of range"));
            }
        }
        let value = group_mtf[j];
        group_mtf.copy_within(0..j, 1);
        group_mtf[0] = value;
        selectors.push(value);
    }
    Ok(selectors)
}

/// Reads the code lengths of each coding group, delta-coded per symbol.
fn read_tables(bits: &mut Bits<'_>, n_groups: usize, alpha: usize) -> Result<Vec<Huffman>> {
    let mut tables = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let mut len = bits.bits(5)? as i32;
        let mut lengths = Vec::with_capacity(alpha);
        for _ in 0..alpha {
            loop {
                if !(1..=MAX_CODE_LEN as i32).contains(&len) {
                    return Err(malformed("code length is out of range"));
                }
                if bits.bit()? == 0 {
                    break;
                }
                len += if bits.bit()? == 0 { 1 } else { -1 };
            }
            lengths.push(len as u8);
        }
        tables.push(Huffman::new(&lengths));
    }
    Ok(tables)
}

/// Decodes the Huffman symbols into the block's pre-BWT bytes (the BWT last column),
/// undoing the zero-run and move-to-front coding.
fn decode_symbols(
    bits: &mut Bits<'_>,
    tables: &[Huffman],
    selectors: &[u8],
    seq_to_unseq: &[u8],
    block_max: usize,
) -> Result<Vec<u8>> {
    // RUNA/RUNB encode zero-run lengths, 2..alpha-1 are MTF positions 1.., and
    // alpha-1 ends the block.
    let end_of_block = (seq_to_unseq.len() + 1) as u16;
    let mut mtf: Vec<u8> = (0..=255_u8).collect();
    let mut tt: Vec<u8> = Vec::with_capacity(block_max.min(1 << 20));
    let mut run_len = 0_usize;
    let mut run_weight = 1_usize;
    let mut group_index = 0_usize;
    let mut group_left = 0_usize;
    let mut table = &tables[0];

    loop {
        if group_left == 0 {
            let &sel =
                selectors.get(group_index).ok_or_else(|| malformed("ran out of selectors"))?;
            table = &tables[usize::from(sel)];
            group_index += 1;
            group_left = GROUP_SIZE;
        }
        group_left -= 1;
        let sym = table.decode(bits)?;

        if sym <= 1 {
            run_len += (usize::from(sym) + 1) * run_weight;
            run_weight <<= 1;
            if run_len > block_max || run_weight > block_max {
                return Err(malformed("zero run exceeds the block size"));
            }
            continue;
        }
        if run_len > 0 {
            let byte = seq_to_unseq[usize::from(mtf[0])];
            if tt.len() + run_len > block_max {
                return Err(malformed("block exceeds its size limit"));
            }
            tt.extend(std::iter::repeat_n(byte, run_len));
            run_len = 0;
            run_weight = 1;
        }
        if sym == end_of_block {
            break;
        }
        // Move the symbol's MTF position to the front.
        let index = usize::from(sym) - 1;
        let value = mtf[index];
        mtf.copy_within(0..index, 1);
        mtf[0] = value;
        if tt.len() >= block_max {
            return Err(malformed("block exceeds its size limit"));
        }
        tt.push(seq_to_unseq[usize::from(value)]);
    }
    Ok(tt)
}

/// Inverts the Burrows-Wheeler transform. `last` is the block's last column and
/// `orig_ptr` the row of the original text.
fn inverse_bwt(last: &[u8], orig_ptr: usize) -> Vec<u8> {
    let n = last.len();
    let mut starts = [0_u32; 256];
    for &b in last {
        starts[usize::from(b)] += 1;
    }
    let mut sum = 0_u32;
    for start in &mut starts {
        let count = *start;
        *start = sum;
        sum += count;
    }
    // Each link holds its successor's row in the high bits and its own byte in the
    // low 8 bits.
    let mut links = vec![0_u32; n];
    for (i, &b) in last.iter().enumerate() {
        let slot = &mut starts[usize::from(b)];
        links[*slot as usize] |= (i as u32) << 8;
        *slot += 1;
    }
    for (link, &b) in links.iter_mut().zip(last) {
        *link |= u32::from(b);
    }
    let mut pos = links[orig_ptr] >> 8;
    let mut text = Vec::with_capacity(n);
    for _ in 0..n {
        let entry = links[pos as usize];
        text.push((entry & 0xFF) as u8);
        pos = entry >> 8;
    }
    text
}

/// Reverses the initial run-length stage: after four equal bytes, the next byte is
/// the number of further copies (0 to 255).
fn undo_initial_rle(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len());
    let mut last = 0_u8;
    let mut run = 0_u32;
    let mut i = 0;
    while i < text.len() {
        let b = text[i];
        i += 1;
        if run == 4 {
            out.extend(std::iter::repeat_n(last, usize::from(b)));
            run = 0;
            continue;
        }
        if run > 0 && b == last {
            run += 1;
        } else {
            last = b;
            run = 1;
        }
        out.push(b);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_matches_the_bzip2_check_value() {
        // The standard CRC-32/BZIP2 check value for "123456789".
        assert_eq!(block_crc(b"123456789"), 0xFC89_1918);
    }

    #[test]
    fn rle_run_and_count() {
        // Four 'a' then a count of 3 gives seven 'a'.
        assert_eq!(undo_initial_rle(&[b'a', b'a', b'a', b'a', 3]), b"aaaaaaa");
        // A count of zero adds nothing, and the next run starts fresh.
        assert_eq!(undo_initial_rle(&[b'a', b'a', b'a', b'a', 0, b'a']), b"aaaaa");
    }

    #[test]
    fn rejects_bad_header() {
        assert!(decompress(b"XYZ9", 10).is_err());
        assert!(decompress(b"", 10).is_err());
    }

    // Fixtures were compressed by libbzip2 (Python's bz2 module) in
    // tests/fixtures/bzip2/generate.py. The plaintext below is rebuilt from the same
    // formulas.

    fn pattern(n: usize) -> Vec<u8> {
        (0..n)
            .map(|i| {
                let i = i as u64;
                b"abcdefghij"[(((i * i) / 17 + i / 64) % 10) as usize]
            })
            .collect()
    }

    fn runs() -> Vec<u8> {
        let mut out = Vec::new();
        for k in 0..20_usize {
            out.extend(std::iter::repeat_n(b'q', 300));
            out.extend(std::iter::repeat_n(b'b', k % 5));
            out.push(65 + k as u8);
        }
        out
    }

    #[test]
    fn decodes_empty_stream() {
        assert_eq!(
            decompress(include_bytes!("../../tests/fixtures/bzip2/empty.bz2"), 0).unwrap(),
            b""
        );
    }

    #[test]
    fn decodes_single_byte() {
        let input = include_bytes!("../../tests/fixtures/bzip2/one.bz2");
        assert_eq!(decompress(input, 1).unwrap(), b"x");
    }

    #[test]
    fn decodes_long_runs_through_the_initial_rle() {
        let input = include_bytes!("../../tests/fixtures/bzip2/runs.bz2");
        assert_eq!(decompress(input, 1 << 20).unwrap(), runs());
    }

    #[test]
    fn decodes_multi_group_pattern() {
        let input = include_bytes!("../../tests/fixtures/bzip2/pattern.bz2");
        assert_eq!(decompress(input, 1 << 20).unwrap(), pattern(3000));
    }

    #[test]
    fn decodes_several_blocks() {
        // Level 1 gives 100,000-byte blocks, so 250,000 bytes need three blocks and
        // the combined stream CRC.
        let input = include_bytes!("../../tests/fixtures/bzip2/blocks.bz2");
        assert_eq!(decompress(input, 1 << 20).unwrap(), pattern(250_000));
    }

    #[test]
    fn decodes_concatenated_streams() {
        let input = include_bytes!("../../tests/fixtures/bzip2/concat.bz2");
        let mut expected = pattern(1000);
        expected.extend_from_slice(b"second stream");
        assert_eq!(decompress(input, 1 << 20).unwrap(), expected);
    }

    #[test]
    fn corrupt_stream_is_an_error_not_a_panic() {
        let good = include_bytes!("../../tests/fixtures/bzip2/pattern.bz2");
        for at in [4, good.len() / 3, good.len() / 2, good.len() - 8] {
            let mut bad = good.to_vec();
            bad[at] ^= 0x10;
            assert!(decompress(&bad, 1 << 20).is_err(), "flip at {at} was accepted");
        }
    }

    #[test]
    fn truncated_stream_is_an_error() {
        let good = include_bytes!("../../tests/fixtures/bzip2/pattern.bz2");
        assert!(decompress(&good[..good.len() / 2], 1 << 20).is_err());
    }

    #[test]
    fn output_limit_is_enforced() {
        let input = include_bytes!("../../tests/fixtures/bzip2/pattern.bz2");
        assert!(decompress(input, 100).is_err());
    }
}
