//! Adaptive arithmetic decoding for TOK3 read-name streams (blob flag 1).
//!
//! Ported from htscodecs `arith_dynamic.c` (`arith_uncompress_to` and the order-0,
//! order-1 and run-length decoders), `c_simple_model.h` and `c_range_coder.h`. These
//! are the streams that samtools writes when a TOK3 blob sets its arithmetic flag,
//! for example with `--output-fmt-option use_arith=1` or the `archive` profile.
//!
//! Each coded stream starts with a flag byte: bits 0 and 1 are the order, bit 2 an
//! external (bzip2) codec, bit 3 a stripe of interleaved sub-streams, bit 4 means the
//! original size is not stored, bit 5 a stored copy, bit 6 run-length coding and bit 7
//! a bit-pack map. Unless bit 4 is set, a uint7 original size follows. Only the decode
//! side is implemented. The bzip2 form and order 2 are reported as
//! [`SeqError::Unsupported`].

use super::tok3::read_uint7;
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM TOK3 arithmetic stream";

const X_PACK: u8 = 0x80;
const X_RLE: u8 = 0x40;
const X_CAT: u8 = 0x20;
const X_NOSZ: u8 = 0x10;
const X_STRIPE: u8 = 0x08;
const X_EXT: u8 = 0x04;
const X_ORDER: u8 = 0x03;

/// Renormalisation threshold of the range coder (`TOP` in `c_range_coder.h`).
const TOP: u32 = 1 << 24;
/// Total frequency above which a model halves its counts (`MAX_FREQ`).
const MAX_FREQ: u32 = (1 << 16) - 17;
/// Increment applied to a symbol each time it is coded (`STEP`).
const STEP: u32 = 16;
/// Run-length symbols 0 to `MAX_RUN - 1`; the maximum value continues a run.
const MAX_RUN: usize = 4;
/// Symbols in the run-length model: the byte values, plus two run contexts.
const RUN_SYMBOLS: usize = 258;
const BYTE_SYMBOLS: usize = 256;
/// Largest original size accepted, as in the reference (`INT_MAX`).
const MAX_SIZE: usize = i32::MAX as usize;

/// Decodes one arithmetic-coded stream whose original size must not exceed `size`.
///
/// The caller compares the returned length with the size it expects. This follows
/// `arith_decode` in the reference, which passes the declared size as the buffer
/// capacity.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for the bzip2 and order-2 forms, and
/// [`SeqError::Malformed`] for truncated or inconsistent data.
pub(super) fn decode(input: &[u8], size: usize) -> Result<Vec<u8>> {
    uncompress(input, size)
}

fn malformed(reason: impl Into<String>) -> SeqError {
    SeqError::malformed(FORMAT, reason)
}

fn truncated() -> SeqError {
    malformed("stream is truncated")
}

/// Decodes a stream whose original size is at most `capacity` (`arith_uncompress_to`).
fn uncompress(input: &[u8], capacity: usize) -> Result<Vec<u8>> {
    let &flags = input.first().ok_or_else(truncated)?;
    if flags & X_STRIPE != 0 {
        return uncompress_striped(input, capacity);
    }
    let mut pos = 1;
    let size = if flags & X_NOSZ != 0 {
        capacity
    } else {
        let (value, used) = read_uint7(&input[pos..])?;
        pos += used;
        value
    };
    if size > capacity || size > MAX_SIZE {
        return Err(malformed(format!("stream declares {size} bytes, room for {capacity}")));
    }

    let mut packed_len = size;
    let mut pack = None;
    if flags & X_PACK != 0 {
        let (used, nsym, map) = unpack_meta(&input[pos..])?;
        pos += used;
        let (value, used) = read_uint7(&input[pos..])?;
        pos += used;
        if value > size {
            return Err(malformed("packed size exceeds the original size"));
        }
        packed_len = value;
        pack = Some((nsym, map));
    }

    let rest = &input[pos..];
    let coded = if rest.is_empty() {
        Vec::new()
    } else if flags & X_CAT != 0 {
        if packed_len > rest.len() {
            return Err(truncated());
        }
        rest[..packed_len].to_vec()
    } else if flags & X_EXT != 0 {
        return Err(SeqError::Unsupported(
            "TOK3 arithmetic streams with bzip2 (external) coding".into(),
        ));
    } else {
        match (flags & X_ORDER, flags & X_RLE != 0) {
            (0, false) => decode_order(rest, packed_len, false)?,
            (1, false) => decode_order(rest, packed_len, true)?,
            (0, true) => decode_run_length(rest, packed_len, false)?,
            (1, true) => decode_run_length(rest, packed_len, true)?,
            _ => return Err(SeqError::Unsupported("TOK3 arithmetic streams of order 2".into())),
        }
    };

    match pack {
        None => Ok(coded),
        Some((nsym, map)) => unpack(&coded, size, nsym, &map),
    }
}

/// Decodes a striped stream: N sub-streams whose bytes are dealt out round-robin.
fn uncompress_striped(input: &[u8], capacity: usize) -> Result<Vec<u8>> {
    let (size, used) = read_uint7(&input[1..])?;
    let mut pos = 1 + used;
    let &n = input.get(pos).ok_or_else(truncated)?;
    pos += 1;
    let n = usize::from(n);
    if n == 0 {
        return Err(malformed("stripe has no sub-streams"));
    }
    if size != capacity {
        return Err(malformed(format!("stripe declares {size} bytes, expected {capacity}")));
    }

    let mut lengths = Vec::with_capacity(n);
    for _ in 0..n {
        let (length, used) = read_uint7(input.get(pos..).ok_or_else(truncated)?)?;
        pos += used;
        if length == 0 {
            return Err(malformed("stripe sub-stream is empty"));
        }
        lengths.push(length);
    }
    let end = lengths
        .iter()
        .try_fold(pos, |acc, &length| acc.checked_add(length))
        .filter(|&end| end <= input.len())
        .ok_or_else(truncated)?;

    let mut streams = Vec::with_capacity(n);
    let mut start = pos;
    for (i, &length) in lengths.iter().enumerate() {
        let sub_size = size / n + usize::from(size % n > i);
        let sub = uncompress(&input[start..start + length], sub_size)?;
        if sub.len() != sub_size {
            return Err(malformed("stripe sub-stream has the wrong length"));
        }
        streams.push(sub);
        start += length;
    }
    debug_assert_eq!(start, end);

    let mut out = vec![0u8; size];
    let mut next = vec![0usize; n];
    for (j, byte) in out.iter_mut().enumerate() {
        let k = j % n;
        *byte = streams[k][next[k]];
        next[k] += 1;
    }
    Ok(out)
}

/// Reads the pack metadata (`hts_unpack_meta`). Returns the bytes used, the number of
/// packed symbols per byte class (0, 2, 4, 8, or 1 for no packing) and the symbol map.
fn unpack_meta(data: &[u8]) -> Result<(usize, usize, [u8; 256])> {
    let &first = data.first().ok_or_else(truncated)?;
    let mut map = [0u8; 256];
    let count = if first == 0 { 256 } else { usize::from(first) };
    let nsym = match count {
        1 => 0,
        2 => 8,
        3 | 4 => 4,
        5..=16 => 2,
        _ => return Ok((1, 1, map)),
    };
    if data.len() < 1 + count {
        return Err(malformed("pack map is truncated"));
    }
    map[..count].copy_from_slice(&data[1..=count]);
    Ok((1 + count, nsym, map))
}

/// Expands packed symbols back to `out_len` bytes (`hts_unpack`). Symbols are packed
/// least significant bits first.
fn unpack(data: &[u8], out_len: usize, nsym: usize, map: &[u8; 256]) -> Result<Vec<u8>> {
    let (per, bits) = match nsym {
        0 => return Ok(vec![map[0]; out_len]),
        1 => return Ok(data.to_vec()),
        8 => (8, 1),
        4 => (4, 2),
        2 => (2, 4),
        _ => return Err(malformed("invalid pack symbol count")),
    };
    if out_len.div_ceil(per) > data.len() {
        return Err(truncated());
    }
    let mask = (1u8 << bits) - 1;
    Ok((0..out_len)
        .map(|i| map[usize::from((data[i / per] >> ((i % per) * bits)) & mask)])
        .collect())
}

/// Reads the symbol count that opens an order-0 or order-1 coded stream.
fn symbol_limit(input: &[u8]) -> Result<usize> {
    let &first = input.first().ok_or_else(truncated)?;
    Ok(if first == 0 { 256 } else { usize::from(first) })
}

/// Order-0 (`order_one` false) or order-1 (`order_one` true) decoding of `len` bytes.
fn decode_order(input: &[u8], len: usize, order_one: bool) -> Result<Vec<u8>> {
    let max_sym = symbol_limit(input)?;
    let mut models = new_byte_models(order_one, max_sym);
    let mut rc = RangeDecoder::new(&input[1..]);
    let mut out = Vec::with_capacity(len);
    let mut last = 0usize;
    for _ in 0..len {
        let sym = models[if order_one { last } else { 0 }].decode(&mut rc)?;
        out.push(sym as u8);
        last = usize::from(sym);
    }
    rc.finish()?;
    Ok(out)
}

/// Run-length decoding of `len` bytes (`arith_uncompress_O0_RLE` and
/// `arith_uncompress_O1_RLE`). Each literal is followed by a run length, coded in a
/// context of the literal's value.
fn decode_run_length(input: &[u8], len: usize, order_one: bool) -> Result<Vec<u8>> {
    let max_sym = symbol_limit(input)?;
    let mut literal_models = new_byte_models(order_one, max_sym);
    let mut run_models: Vec<Model> =
        (0..RUN_SYMBOLS).map(|_| Model::new(RUN_SYMBOLS, MAX_RUN)).collect();
    let mut rc = RangeDecoder::new(&input[1..]);
    let mut out = vec![0u8; len];
    let mut prev = 0usize;
    let mut i = 0usize;
    while i < len {
        let lit = usize::from(literal_models[if order_one { prev } else { 0 }].decode(&mut rc)?);
        out[i] = lit as u8;

        let mut run = 0usize;
        let mut ctx = lit;
        loop {
            let r = usize::from(run_models[ctx].decode(&mut rc)?);
            if ctx == lit {
                ctx = BYTE_SYMBOLS;
            } else if ctx < RUN_SYMBOLS - 1 {
                ctx += 1;
            }
            run += r;
            if r != MAX_RUN - 1 || run >= len {
                break;
            }
        }

        let mut left = run;
        while left > 0 && i + 1 < len {
            left -= 1;
            i += 1;
            out[i] = lit as u8;
        }
        prev = lit;
        i += 1;
    }
    rc.finish()?;
    Ok(out)
}

/// One byte model per order-1 context, or a single model for order 0.
fn new_byte_models(order_one: bool, max_sym: usize) -> Vec<Model> {
    let count = if order_one { BYTE_SYMBOLS } else { 1 };
    (0..count).map(|_| Model::new(BYTE_SYMBOLS, max_sym)).collect()
}

/// Adaptive frequency model (`SIMPLE_MODEL` in `c_simple_model.h`).
struct Model {
    total: u32,
    nsym: usize,
    /// Index 0 is a sentinel with the largest frequency. Index `i + 1` holds the
    /// reference's `F[i]` as `(frequency, symbol)`. The model keeps the entries
    /// roughly sorted by frequency.
    freqs: Vec<(u32, u16)>,
}

impl Model {
    /// A model over `nsym` symbols where only the first `max_sym` are initially used.
    fn new(nsym: usize, max_sym: usize) -> Self {
        let mut freqs = vec![(0u32, 0u16); nsym + 2];
        freqs[0] = (MAX_FREQ, 0);
        for (i, entry) in freqs.iter_mut().skip(1).enumerate() {
            *entry = (u32::from(i < max_sym), i as u16);
        }
        Self { total: max_sym as u32, nsym, freqs }
    }

    /// Halves every frequency, stopping at the first zero entry.
    fn normalize(&mut self) {
        self.total = 0;
        for entry in &mut self.freqs[1..] {
            if entry.0 == 0 {
                break;
            }
            entry.0 -= entry.0 >> 1;
            self.total += entry.0;
        }
    }

    /// Decodes one symbol and adapts the model (`SIMPLE_MODEL_decodeSymbol`).
    fn decode(&mut self, rc: &mut RangeDecoder) -> Result<u16> {
        let freq = rc.get_freq(self.total);
        if freq > MAX_FREQ {
            return Err(malformed("arithmetic code is out of range"));
        }
        let mut idx = 1;
        let mut acc = 0u32;
        loop {
            if idx > self.nsym {
                return Err(malformed("arithmetic code is out of range"));
            }
            acc += self.freqs[idx].0;
            if acc > freq {
                break;
            }
            idx += 1;
        }
        let (symbol_freq, sym) = self.freqs[idx];
        rc.decode(acc - symbol_freq, symbol_freq);

        self.freqs[idx].0 += STEP;
        self.total += STEP;
        if self.total > MAX_FREQ {
            self.normalize();
        }
        if self.freqs[idx].0 > self.freqs[idx - 1].0 {
            self.freqs.swap(idx, idx - 1);
        }
        Ok(sym)
    }
}

/// Range decoder (`RangeCoder` in `c_range_coder.h`), reading bytes from a slice.
struct RangeDecoder<'a> {
    data: &'a [u8],
    pos: usize,
    code: u32,
    range: u32,
    failed: bool,
}

impl<'a> RangeDecoder<'a> {
    fn new(data: &'a [u8]) -> Self {
        let mut rc = Self { data, pos: 0, code: 0, range: u32::MAX, failed: false };
        if data.len() < 5 {
            rc.pos = data.len();
        } else {
            for _ in 0..5 {
                rc.code = (rc.code << 8) | u32::from(data[rc.pos]);
                rc.pos += 1;
            }
        }
        rc
    }

    /// Scales the range to `total` and returns the code's position within it.
    fn get_freq(&mut self, total: u32) -> u32 {
        if total != 0 && self.range >= total {
            self.range /= total;
            self.code / self.range
        } else {
            0
        }
    }

    /// Consumes the interval `[cum, cum + freq)` and renormalises.
    fn decode(&mut self, cum: u32, freq: u32) {
        self.code = self.code.wrapping_sub(cum.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(freq);
        while self.range < TOP {
            let Some(&byte) = self.data.get(self.pos) else {
                self.failed = true;
                return;
            };
            self.pos += 1;
            self.code = (self.code << 8) | u32::from(byte);
            self.range <<= 8;
        }
    }

    /// Fails if the decoder ran out of input (`RC_FinishDecode`).
    fn finish(&self) -> Result<()> {
        if self.failed { Err(truncated()) } else { Ok(()) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Stored (`X_CAT`) stream: flags, uint7 size, then the bytes.
    fn stored(data: &[u8]) -> Vec<u8> {
        let mut out = vec![X_CAT, data.len() as u8];
        out.extend_from_slice(data);
        out
    }

    #[test]
    fn packed_stored_stream_expands_one_bit_symbols() {
        // Map of two symbols, "A" and "C", packed one bit each, least significant first.
        // Bits 0,1,1,0,1,0,0,1 pack to 0x96.
        let mut data = vec![X_CAT | X_PACK, 8, 2, b'A', b'C', 1, 0x96];
        assert_eq!(decode(&data, 8).unwrap(), b"ACCACAAC");
        // A packed size larger than the original size is rejected.
        data[5] = 9;
        assert!(decode(&data, 8).is_err());
    }

    #[test]
    fn single_symbol_pack_fills_the_stream() {
        // Map of one symbol: the stream is empty and every byte is that symbol.
        let data = [X_PACK, 4, 1, b'x', 0];
        assert_eq!(decode(&data, 4).unwrap(), b"xxxx");
    }

    #[test]
    fn striped_stream_deals_sub_streams_round_robin() {
        // Five bytes in two stripes: "ace" and "bd", each stored raw.
        let sub0 = stored(b"ace");
        let sub1 = stored(b"bd");
        let mut data = vec![X_STRIPE, 5, 2, sub0.len() as u8, sub1.len() as u8];
        data.extend_from_slice(&sub0);
        data.extend_from_slice(&sub1);
        assert_eq!(decode(&data, 5).unwrap(), b"abcde");
        // The declared size must match what the caller expects.
        assert!(decode(&data, 6).is_err());
    }

    #[test]
    fn stored_stream_with_no_size_takes_the_capacity() {
        let data = [X_CAT | X_NOSZ, b'h', b'i'];
        assert_eq!(decode(&data, 2).unwrap(), b"hi");
    }

    #[test]
    fn bzip2_and_order_two_are_unsupported() {
        assert!(matches!(decode(&[X_EXT, 0, 0], 0), Err(SeqError::Unsupported(_))));
        assert!(matches!(decode(&[2, 0, 0], 0), Err(SeqError::Unsupported(_))));
    }

    #[test]
    fn truncated_coded_stream_is_an_error() {
        // Order-0 body of 2 bytes after the symbol count: asking for 200 bytes makes
        // the range decoder renormalise past the end of its input.
        assert!(decode(&[0, 200, 1, 0, 0], 200).is_err());
    }
}
