//! rANS Nx16 entropy codec (CRAM 3.1 codec method 5), following htscodecs.
//!
//! A stream starts with a flags byte. Lengths are big-endian base-128 varints: seven
//! bits per byte, most significant group first, with the high bit set on every byte
//! but the last. The bit values of the flags:
//!
//! | Bit  | Name      | Support |
//! |------|-----------|---------|
//! | 0x80 | `PACK`    | symbols packed 1, 2 or 4 to a byte, with a map table |
//! | 0x40 | `RLE`     | runs of the symbols in a map, with run lengths |
//! | 0x20 | `CAT`     | stored, uncompressed |
//! | 0x10 | `NO_SIZE` | size taken from the caller's expected length |
//! | 0x08 | `STRIPE`  | `N` interleaved sub-streams (each a complete stream) |
//! | 0x04 | `N32`     | 32 interleaved states instead of 4 |
//! | 0x01 | `ORDER`   | order-1 coder: the context is the previous symbol |
//!
//! The order-0 body is a run-length symbol list with a frequency for each symbol,
//! then the initial 32-bit states, then 16-bit renormalization words. Symbol `i` is
//! coded by state `i % lanes`. Frequencies are normalized by powers of two up to 4096.
//! The order-1 body starts with a byte holding the table shift (high nibble, 12 or
//! 10) and a flag for a compressed table (bit 0). Lane `z` codes the segment starting
//! at `z * (size / lanes)`, and the last lane also codes the remainder.
//! The PACK, RLE and CAT transforms apply in that order on decode: rANS (or CAT)
//! data, then RLE expansion, then unpacking.

use crate::cram::bytes::Cursor;
use crate::cram::rans::{Table, normalize, walk_symbols, write_symbols};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM rANS Nx16";
const FLAG_ORDER: u8 = 0x01;
const FLAG_N32: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

const TF_SHIFT: u32 = 12;
const TF_SHIFT_FAST: u32 = 10;
const TOTAL: u32 = 1 << TF_SHIFT;
/// Renormalization bound for 16-bit words.
const LOWER16: u32 = 1 << 15;
const N_STATES: usize = 4;
const N_STATES_32: usize = 32;
/// Largest output accepted, as a guard against corrupt sizes.
const MAX_OUTPUT: usize = 1 << 28;
/// Stripes and RLE and pack maps are capped, as a guard against corrupt counts.
const MAX_STRIPES: usize = 255;

/// Decodes an Nx16 payload. `expected_len` is the block's raw size. It is used
/// for `NO_SIZE` streams and must match the declared size otherwise.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for truncated input, size mismatches, or tables
/// that do not normalize.
pub fn decode(input: &[u8], expected_len: usize) -> Result<Vec<u8>> {
    let mut c = Cursor::new(input);
    decode_stream(&mut c, Some(expected_len))
}

/// Decodes one complete stream: a flags byte and its payload.
fn decode_stream(c: &mut Cursor<'_>, expected: Option<usize>) -> Result<Vec<u8>> {
    let flags = c.u8()?;
    if flags & FLAG_STRIPE != 0 {
        return decode_striped(c, expected);
    }
    let lanes = if flags & FLAG_N32 != 0 { N_STATES_32 } else { N_STATES };
    let order1 = flags & FLAG_ORDER != 0;
    let size = if flags & FLAG_NO_SIZE != 0 {
        expected.ok_or_else(|| SeqError::malformed(FORMAT, "NO_SIZE stream without a length"))?
    } else {
        let declared = read_varint(c)?;
        if expected.is_some_and(|e| e != declared) {
            return Err(SeqError::malformed(
                FORMAT,
                format!("stream declares {declared} bytes, block declares {expected:?}"),
            ));
        }
        declared
    };
    if size > MAX_OUTPUT {
        return Err(SeqError::malformed(FORMAT, format!("declared size {size} is too large")));
    }
    if size == 0 {
        return Ok(Vec::new());
    }

    // Transform metadata comes first: the pack map, then the RLE meta data.
    let pack = if flags & FLAG_PACK != 0 { Some(read_pack_map(c)?) } else { None };
    let packed_len = if pack.is_some() { read_varint(c)? } else { size };
    let rle = if flags & FLAG_RLE != 0 { Some(read_rle_meta(c)?) } else { None };

    // Stage one is rANS (or stored) data, of the length the transforms leave.
    let stage1_len = rle.as_ref().map_or(packed_len, |r| r.literal_len);
    let stage1 = if flags & FLAG_CAT != 0 {
        c.bytes(stage1_len)?.to_vec()
    } else {
        decode_body(c, stage1_len, order1, lanes)?
    };
    let stage2 = match &rle {
        Some(r) => {
            let target = if pack.is_some() { packed_len } else { size };
            unrle(&stage1, r, target)?
        }
        None => stage1,
    };
    let out = match &pack {
        Some(map) => unpack(&stage2, map, size)?,
        None => stage2,
    };
    if out.len() != size {
        return Err(SeqError::malformed(FORMAT, "decoded length does not match the header"));
    }
    Ok(out)
}

/// Decodes a striped stream: `N` sub-streams, interleaved byte by byte.
fn decode_striped(c: &mut Cursor<'_>, expected: Option<usize>) -> Result<Vec<u8>> {
    let total = read_varint(c)?;
    if expected.is_some_and(|e| e != total) {
        return Err(SeqError::malformed(FORMAT, "striped length does not match the block"));
    }
    if total > MAX_OUTPUT {
        return Err(SeqError::malformed(FORMAT, "striped length is too large"));
    }
    let n = usize::from(c.u8()?);
    if n == 0 || n > MAX_STRIPES {
        return Err(SeqError::malformed(FORMAT, "bad stripe count"));
    }
    let mut lengths = Vec::with_capacity(n);
    for k in 0..n {
        // Stripe k holds every n-th byte starting at k.
        lengths.push(total / n + usize::from(total % n > k));
    }
    let mut sizes = Vec::with_capacity(n);
    for _ in 0..n {
        sizes.push(read_varint(c)?);
    }
    let mut out = vec![0u8; total];
    for (k, (&clen, &ulen)) in sizes.iter().zip(&lengths).enumerate() {
        let body = c.bytes(clen)?;
        let stripe = decode_stream(&mut Cursor::new(body), Some(ulen))?;
        for (l, &b) in stripe.iter().enumerate() {
            out[l * n + k] = b;
        }
    }
    Ok(out)
}

/// The pack map: the symbols in use and how many fit in one byte.
struct PackMap {
    /// Symbols per byte: 0 (constant), 1 (no packing), 2, 4 or 8.
    per_byte: usize,
    symbols: Vec<u8>,
}

fn read_pack_map(c: &mut Cursor<'_>) -> Result<PackMap> {
    let n = match usize::from(c.u8()?) {
        0 => 256,
        n => n,
    };
    let per_byte = match n {
        1 => 0,
        2 => 8,
        3..=4 => 4,
        5..=16 => 2,
        _ => return Ok(PackMap { per_byte: 1, symbols: Vec::new() }),
    };
    Ok(PackMap { per_byte, symbols: c.bytes(n)?.to_vec() })
}

/// Unpacks `packed` into `out_len` symbols. Each byte holds `per_byte` symbols, least
/// significant bits first.
fn unpack(packed: &[u8], map: &PackMap, out_len: usize) -> Result<Vec<u8>> {
    match map.per_byte {
        1 => {
            if packed.len() != out_len {
                return Err(SeqError::malformed(FORMAT, "unpacked length does not match"));
            }
            Ok(packed.to_vec())
        }
        0 => {
            let sym = map.symbols.first().copied().unwrap_or(0);
            Ok(vec![sym; out_len])
        }
        per_byte => {
            let bits = 8 / per_byte as u32;
            let mask = (1u8 << bits) - 1;
            if packed.len() < out_len.div_ceil(per_byte) {
                return Err(SeqError::malformed(FORMAT, "packed data is too short"));
            }
            let mut out = Vec::with_capacity(out_len);
            'bytes: for &b in packed {
                for k in 0..per_byte {
                    if out.len() == out_len {
                        break 'bytes;
                    }
                    let idx = usize::from((b >> (k as u32 * bits)) & mask);
                    let sym =
                        map.symbols.get(idx).copied().ok_or_else(|| {
                            SeqError::malformed(FORMAT, "pack index out of range")
                        })?;
                    out.push(sym);
                }
            }
            Ok(out)
        }
    }
}

/// The run-length meta data: symbols that repeat, and their run lengths.
struct RleMeta {
    /// Length of the literal stream the rANS stage decodes to.
    literal_len: usize,
    /// Symbols that are run-length coded.
    symbols: Vec<u8>,
    /// Varint run lengths, one per coded symbol occurrence.
    runs: Vec<u8>,
}

fn read_rle_meta(c: &mut Cursor<'_>) -> Result<RleMeta> {
    let u_meta = read_varint(c)?;
    let literal_len = read_varint(c)?;
    let meta: Vec<u8> = if u_meta & 1 == 1 {
        // Stored meta data of u_meta / 2 bytes.
        c.bytes(u_meta / 2)?.to_vec()
    } else {
        // rANS-compressed meta data, which decodes to u_meta / 2 bytes.
        let c_meta = read_varint(c)?;
        let body = c.bytes(c_meta)?;
        decode_o0_body(&mut Cursor::new(body), u_meta / 2, N_STATES)?
    };
    let n = match meta.first().copied().map(usize::from) {
        None => return Err(SeqError::malformed(FORMAT, "empty RLE meta data")),
        Some(0) => 256,
        Some(n) => n,
    };
    if meta.len() < 1 + n {
        return Err(SeqError::malformed(FORMAT, "RLE symbol list is truncated"));
    }
    Ok(RleMeta { literal_len, symbols: meta[1..=n].to_vec(), runs: meta[1 + n..].to_vec() })
}

/// Expands the run-length coded literals. Each literal that is a run symbol is followed
/// by a varint `r`: the symbol then repeats `r + 1` times.
fn unrle(literals: &[u8], rle: &RleMeta, out_len: usize) -> Result<Vec<u8>> {
    let mut saved = [false; 256];
    for &s in &rle.symbols {
        saved[usize::from(s)] = true;
    }
    let mut runs = Cursor::new(&rle.runs);
    let mut out = Vec::with_capacity(out_len.min(1 << 20));
    for &b in literals {
        if saved[usize::from(b)] {
            let r = read_varint(&mut runs)?;
            if out.len() + r + 1 > out_len {
                return Err(SeqError::malformed(FORMAT, "run-length data overflows the output"));
            }
            out.extend(std::iter::repeat_n(b, r + 1));
        } else {
            if out.len() == out_len {
                return Err(SeqError::malformed(FORMAT, "run-length data overflows the output"));
            }
            out.push(b);
        }
    }
    if out.len() != out_len {
        return Err(SeqError::malformed(FORMAT, "run-length output is the wrong length"));
    }
    Ok(out)
}

/// Decodes a rANS body of `out_len` bytes with `lanes` interleaved states.
fn decode_body(c: &mut Cursor<'_>, out_len: usize, order1: bool, lanes: usize) -> Result<Vec<u8>> {
    if order1 { decode_o1_body(c, out_len, lanes) } else { decode_o0_body(c, out_len, lanes) }
}

/// A symbol lookup for frequencies that sum to `1 << shift`.
struct Lut {
    freq: [u32; 256],
    cum: [u32; 257],
    lookup: Vec<u8>,
}

impl Lut {
    fn new(freq: &[u32; 256], shift: u32) -> Self {
        let mut cum = [0u32; 257];
        for s in 0..256 {
            cum[s + 1] = cum[s] + freq[s];
        }
        let mut lookup = vec![0u8; 1 << shift];
        for s in 0..256 {
            for slot in cum[s]..cum[s + 1] {
                lookup[slot as usize] = s as u8;
            }
        }
        Self { freq: *freq, cum, lookup }
    }
}

/// Decodes one symbol from `state` with the table `lut` at `shift` bits, then
/// renormalizes the state from one 16-bit word if it fell below `LOWER16`.
fn decode_symbol(state: &mut u32, c: &mut Cursor<'_>, shift: u32, lut: &Lut) -> Result<u8> {
    let slot = *state & ((1 << shift) - 1);
    let sym = lut.lookup[slot as usize];
    let s = usize::from(sym);
    let mut x = lut.freq[s] * (*state >> shift) + slot - lut.cum[s];
    if x < LOWER16 {
        x = (x << 16) | read_word(c)?;
    }
    *state = x;
    Ok(sym)
}

/// Decodes an order-0 body of `out_len` bytes: table, initial states, then the data.
fn decode_o0_body(c: &mut Cursor<'_>, out_len: usize, lanes: usize) -> Result<Vec<u8>> {
    if out_len == 0 {
        return Ok(Vec::new());
    }
    let freq = read_normalized_freq(c)?;
    let lut = Lut::new(&freq, TF_SHIFT);
    let mut state = vec![0u32; lanes];
    for s in &mut state {
        *s = c.u32_le()?;
    }
    let mut out = Vec::with_capacity(out_len.min(1 << 20));
    for i in 0..out_len {
        out.push(decode_symbol(&mut state[i % lanes], c, TF_SHIFT, &lut)?);
    }
    Ok(out)
}

/// Decodes an order-1 body of `out_len` bytes. Lane `z` codes the segment starting at
/// `z * (out_len / lanes)`, and the last lane also codes the remainder. Each lane's
/// context is the symbol it coded last, starting from 0.
fn decode_o1_body(c: &mut Cursor<'_>, out_len: usize, lanes: usize) -> Result<Vec<u8>> {
    if out_len == 0 {
        return Ok(Vec::new());
    }
    let header = c.u8()?;
    let shift = u32::from(header >> 4);
    if shift != TF_SHIFT && shift != TF_SHIFT_FAST {
        return Err(SeqError::malformed(FORMAT, format!("order-1 shift {shift} is not 10 or 12")));
    }
    let luts = if header & 1 == 1 {
        // The frequency tables are an order-0 body of u_len bytes, stored compressed.
        let u_len = read_varint(c)?;
        let c_len = read_varint(c)?;
        let raw = decode_o0_body(&mut Cursor::new(c.bytes(c_len)?), u_len, N_STATES)?;
        read_o1_tables(&mut Cursor::new(&raw), shift)?
    } else {
        read_o1_tables(c, shift)?
    };
    let mut state = vec![0u32; lanes];
    for s in &mut state {
        *s = c.u32_le()?;
    }

    let seg = out_len / lanes;
    let mut out = vec![0u8; out_len];
    let mut ctx = vec![0u8; lanes];
    for t in 0..seg {
        for z in 0..lanes {
            let lut = luts[usize::from(ctx[z])].as_ref().ok_or_else(no_table)?;
            let sym = decode_symbol(&mut state[z], c, shift, lut)?;
            out[z * seg + t] = sym;
            ctx[z] = sym;
        }
    }
    let last = lanes - 1;
    for slot in &mut out[lanes * seg..] {
        let lut = luts[usize::from(ctx[last])].as_ref().ok_or_else(no_table)?;
        let sym = decode_symbol(&mut state[last], c, shift, lut)?;
        *slot = sym;
        ctx[last] = sym;
    }
    Ok(out)
}

fn no_table() -> SeqError {
    SeqError::malformed(FORMAT, "order-1 context has no frequency table")
}

/// Reads the order-1 tables: the symbols that can start a context, then for each of
/// them the frequencies of the symbols that follow it. Zero frequencies are run-length
/// coded: a zero frequency is followed by a count of further zeros. Each table is
/// scaled by powers of two to `1 << shift`. A context with no frequencies has no table.
fn read_o1_tables(c: &mut Cursor<'_>, shift: u32) -> Result<Vec<Option<Lut>>> {
    let mut present = [false; 256];
    walk_symbols(c, |_, s| {
        present[usize::from(s)] = true;
        Ok(())
    })?;
    let total = 1u64 << shift;
    let mut luts: Vec<Option<Lut>> = (0..256).map(|_| None).collect();
    for ctx in 0..256 {
        if !present[ctx] {
            continue;
        }
        let mut freq = [0u32; 256];
        let mut zeros = 0usize;
        for j in 0..256 {
            if !present[j] {
                continue;
            }
            if zeros > 0 {
                zeros -= 1;
            } else {
                let f = u32::try_from(read_varint(c)?)
                    .map_err(|_| SeqError::malformed(FORMAT, "frequency too large"))?;
                freq[j] = f;
                if f == 0 {
                    zeros = usize::from(c.u8()?);
                }
            }
        }
        let mut sum: u64 = freq.iter().map(|&f| u64::from(f)).sum();
        if sum == 0 {
            continue;
        }
        if sum > total {
            return Err(SeqError::malformed(FORMAT, "order-1 frequencies exceed the table size"));
        }
        while sum < total {
            sum <<= 1;
            for f in &mut freq {
                *f <<= 1;
            }
        }
        if sum != total {
            return Err(SeqError::malformed(FORMAT, "order-1 frequencies do not normalize"));
        }
        luts[ctx] = Some(Lut::new(&freq, shift));
    }
    Ok(luts)
}

/// Reads the alphabet and the varint frequencies, then scales the table by powers of two
/// until it sums to 4096. Any other total is malformed.
fn read_normalized_freq(c: &mut Cursor<'_>) -> Result<[u32; 256]> {
    // The alphabet comes first, then one frequency per symbol in ascending order.
    let mut syms = Vec::new();
    walk_symbols(c, |_, s| {
        syms.push(s);
        Ok(())
    })?;
    let mut freq = [0u32; 256];
    for &s in &syms {
        freq[usize::from(s)] = u32::try_from(read_varint(c)?)
            .map_err(|_| SeqError::malformed(FORMAT, "frequency too large"))?;
    }
    let mut sum: u64 = freq.iter().map(|&f| u64::from(f)).sum();
    if sum == 0 {
        return Err(SeqError::malformed(FORMAT, "empty frequency table"));
    }
    while sum < u64::from(TOTAL) {
        sum <<= 1;
        for f in &mut freq {
            *f <<= 1;
        }
    }
    if sum != u64::from(TOTAL) {
        return Err(SeqError::malformed(FORMAT, "frequencies do not normalize to 4096"));
    }
    Ok(freq)
}

/// Reads one 16-bit renormalization word, little-endian.
fn read_word(c: &mut Cursor<'_>) -> Result<u32> {
    let b = c.bytes(2)?;
    Ok(u32::from(u16::from_le_bytes([b[0], b[1]])))
}

/// Reads a big-endian base-128 varint.
fn read_varint(c: &mut Cursor<'_>) -> Result<usize> {
    let mut value = 0usize;
    for _ in 0..10 {
        let b = c.u8()?;
        value = value
            .checked_mul(128)
            .and_then(|v| v.checked_add(usize::from(b & 0x7F)))
            .ok_or_else(|| SeqError::malformed(FORMAT, "varint overflows"))?;
        if b & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(SeqError::malformed(FORMAT, "varint too long"))
}

/// Appends a big-endian base-128 varint.
fn put_varint(out: &mut Vec<u8>, v: usize) {
    let mut groups = Vec::new();
    let mut v = v;
    loop {
        groups.push((v & 0x7F) as u8);
        v >>= 7;
        if v == 0 {
            break;
        }
    }
    // Most significant group first; every byte but the last carries the high bit.
    let last = groups.len() - 1;
    for (i, g) in groups.iter().rev().enumerate() {
        out.push(if i < last { *g | 0x80 } else { *g });
    }
}

/// Encodes `data` as an order-0, four-state Nx16 stream with no transforms.
#[must_use]
pub fn encode(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8];
    put_varint(&mut out, data.len());
    if data.is_empty() {
        return out;
    }
    let mut counts = [0u64; 256];
    for &b in data {
        counts[usize::from(b)] += 1;
    }
    let freq = normalize(&counts);
    let table = Table::new(&freq);
    let syms: Vec<u8> = (0..256usize).filter(|&s| freq[s] > 0).map(|s| s as u8).collect();
    write_symbols(&syms, &mut out, |_, _| {});
    for &s in &syms {
        put_varint(&mut out, freq[usize::from(s)] as usize);
    }

    let mut state = [LOWER16; N_STATES];
    let mut emitted: Vec<u16> = Vec::new();
    for i in (0..data.len()).rev() {
        let s = usize::from(data[i]);
        let f = table.freq[s];
        let x = &mut state[i % N_STATES];
        let x_max = ((LOWER16 >> TF_SHIFT) << 16) * f;
        while *x >= x_max {
            emitted.push(*x as u16);
            *x >>= 16;
        }
        *x = ((*x / f) << TF_SHIFT) + (*x % f) + table.cum[s];
    }
    for x in &state {
        out.extend_from_slice(&x.to_le_bytes());
    }
    for w in emitted.iter().rev() {
        out.extend_from_slice(&w.to_le_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed | 1;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACGTN"[((s >> 33) % 5) as usize]
            })
            .collect()
    }

    #[test]
    fn varint_is_big_endian() {
        // 8964 is 70 << 7 | 4, so it is written c6 04, as in the samtools fixture.
        let mut out = Vec::new();
        put_varint(&mut out, 8964);
        assert_eq!(out, [0xC6, 0x04]);
        assert_eq!(read_varint(&mut Cursor::new(&[0xC6, 0x04])).unwrap(), 8964);
        out.clear();
        put_varint(&mut out, 0x7F);
        assert_eq!(out, [0x7F]);
    }

    #[test]
    fn roundtrip_order0() {
        let inputs = [
            Vec::new(),
            b"A".to_vec(),
            b"GAT".to_vec(),
            sample(4, 1),
            sample(33, 2),
            sample(1000, 3),
            (0..=255u8).cycle().take(5000).collect(),
            vec![b'C'; 777],
        ];
        for data in &inputs {
            let enc = encode(data);
            assert_eq!(&decode(&enc, data.len()).unwrap(), data, "len {}", data.len());
        }
    }

    #[test]
    fn hand_derived_order0_bytes() {
        // "AA": flags 0, size 2, then the table with A = 4096. 4096 is 32 << 7, so it
        // is written a0 00 (32 with the continuation bit, then 0). The terminator is 0.
        // Frequency 4096 leaves every state at LOWER16 = 0x8000.
        let enc = encode(b"AA");
        let mut expected = vec![0x00, 0x02, 65, 0x00, 0xA0, 0x00];
        for _ in 0..4 {
            expected.extend_from_slice(&[0x00, 0x80, 0x00, 0x00]);
        }
        assert_eq!(enc, expected);
        assert_eq!(decode(&enc, 2).unwrap(), b"AA");
    }

    #[test]
    fn cat_and_no_size_flags() {
        let mut cat = vec![FLAG_CAT, 6];
        cat.extend_from_slice(b"stored");
        assert_eq!(decode(&cat, 6).unwrap(), b"stored");
        let mut no_size = encode(b"ACGT");
        no_size[0] |= FLAG_NO_SIZE;
        no_size.remove(1);
        assert_eq!(decode(&no_size, 4).unwrap(), b"ACGT");
    }

    #[test]
    fn pack_with_cat_unpacks_bits_least_significant_first() {
        // flags PACK|CAT, size 8, map [A, C], packed length 1, then 0b1011_0010.
        let bytes = [0xA0, 8, 2, b'A', b'C', 1, 0b1011_0010];
        assert_eq!(decode(&bytes, 8).unwrap(), b"ACAACCAC");
    }

    #[test]
    fn rle_with_cat_expands_runs() {
        // flags RLE|CAT, size 4, stored meta of 3 bytes (u_meta 7), literal length 2.
        // The meta is one run symbol 'A' with run 2 (varint 2): "AAA" then "B".
        let bytes = [0x60, 4, 7, 2, 1, b'A', 2, b'A', b'B'];
        assert_eq!(decode(&bytes, 4).unwrap(), b"AAAB");
    }

    #[test]
    fn stripes_interleave_complete_streams() {
        // Two stripes of "AA" and "BB" give "ABAB".
        let a = encode(b"AA");
        let b = encode(b"BB");
        let mut bytes = vec![FLAG_STRIPE, 4, 2];
        bytes.push(a.len() as u8);
        bytes.push(b.len() as u8);
        bytes.extend_from_slice(&a);
        bytes.extend_from_slice(&b);
        assert_eq!(decode(&bytes, 4).unwrap(), b"ABAB");
    }

    /// The 53-byte text that the 32-state and order-1 vectors below decode to.
    const NX_TEXT: &[u8] = b"ACGTACGTTTGACCAGTAACCGTTAGGCATTTACGAACCGTAAAAAAAAGGGG";

    #[test]
    fn order1_stream_from_htscodecs() {
        // Produced by htscodecs rans_compress_to_4x16 with order 1 (flags 0x01) for the
        // 240 bytes "ACGTTGCAAC" repeated 24 times. The table shift is 10 (header a0).
        let stream = [
            0x01, 0x81, 0x70, 0xa0, 0x00, 0x41, 0x43, 0x47, 0x54, 0x00, 0x00, 0x00, 0x04, 0x00,
            0x02, 0x00, 0x00, 0x2a, 0x56, 0x00, 0x01, 0x00, 0x00, 0x56, 0x00, 0x00, 0x2a, 0x00,
            0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x20, 0x00, 0x02, 0x20, 0x20, 0x35, 0x76, 0xb4,
            0x00, 0x35, 0x76, 0xb4, 0x00, 0x35, 0x76, 0xb4, 0x00, 0x35, 0x76, 0xb4, 0x00, 0xd8,
            0x61, 0xd8, 0x61, 0xd8, 0x61, 0xd8, 0x61, 0x04, 0xc0, 0x04, 0xc0, 0x04, 0xc0, 0x04,
            0xc0, 0xc0, 0x42, 0xc0, 0x42, 0xc0, 0x42, 0xc0, 0x42,
        ];
        let text = b"ACGTTGCAAC".repeat(24);
        assert_eq!(decode(&stream, text.len()).unwrap(), text);
    }

    #[test]
    fn order0_32_states_from_htscodecs() {
        // Flags 0x04 (N32), size 53, then the body htscodecs rans_compress_O0_32x16
        // produced for NX_TEXT. The body is the table, 32 initial states, then words.
        let stream = [
            0x04, 0x35, 0x41, 0x43, 0x47, 0x54, 0x00, 0x18, 0x0c, 0x0f, 0x0d, 0x00, 0x82, 0x03,
            0x00, 0x00, 0x28, 0x0e, 0x00, 0x40, 0x3b, 0x09, 0x00, 0xc0, 0x8c, 0x06, 0x00, 0x00,
            0x82, 0x03, 0x00, 0x00, 0x28, 0x0e, 0x00, 0x40, 0x5a, 0x0b, 0x00, 0xc0, 0xad, 0x0a,
            0x00, 0xc0, 0x4d, 0x0c, 0x00, 0xc0, 0x8c, 0x06, 0x00, 0x80, 0xa9, 0x05, 0x00, 0x00,
            0x82, 0x03, 0x00, 0x00, 0x08, 0x07, 0x00, 0x00, 0x08, 0x07, 0x00, 0x00, 0x82, 0x03,
            0x00, 0x80, 0xa9, 0x05, 0x00, 0xc0, 0x8c, 0x06, 0x00, 0x80, 0xc1, 0x05, 0x00, 0x80,
            0xc1, 0x05, 0x00, 0x80, 0x87, 0x0b, 0x00, 0x80, 0x87, 0x0b, 0x00, 0x80, 0x29, 0x02,
            0x00, 0x00, 0x7e, 0x02, 0x00, 0x00, 0x7e, 0x02, 0x00, 0x00, 0x52, 0x01, 0x00, 0x80,
            0x29, 0x02, 0x00, 0x80, 0x29, 0x02, 0x00, 0x00, 0xa8, 0x02, 0x00, 0x00, 0x52, 0x01,
            0x00, 0x00, 0x7e, 0x02, 0x00, 0x00, 0x7e, 0x02, 0x00, 0x00, 0x7e, 0x02, 0x00,
        ];
        assert_eq!(decode(&stream, NX_TEXT.len()).unwrap(), NX_TEXT);
    }

    #[test]
    fn order1_32_states_from_htscodecs() {
        // Flags 0x05 (N32 and order 1), size 53, then the body htscodecs
        // rans_compress_O1_32x16 produced for NX_TEXT. The table shift is 10.
        let stream = [
            0x05, 0x35, 0xa0, 0x00, 0x41, 0x43, 0x47, 0x54, 0x00, 0x00, 0x00, 0x08, 0x07, 0x07,
            0x0a, 0x00, 0x00, 0x10, 0x0a, 0x05, 0x01, 0x00, 0x00, 0x03, 0x04, 0x09, 0x00, 0x00,
            0x00, 0x00, 0x02, 0x01, 0x05, 0x08, 0x00, 0x00, 0x08, 0x00, 0x00, 0x01, 0x07, 0x00,
            0x00, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x20, 0x4a, 0x02, 0x00, 0x40, 0x9b, 0x01,
            0x00, 0x00, 0x00, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x20, 0x4a, 0x02, 0x00, 0x40,
            0x9b, 0x01, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x20, 0x4a, 0x02,
            0x00, 0x00, 0x00, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x00,
            0x00, 0x02, 0x00, 0x20, 0x4a, 0x02, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x00, 0x00, 0x02,
            0x00, 0x00, 0x00, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x20,
            0x4a, 0x02, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x00, 0x00, 0x02,
            0x00, 0x20, 0x4a, 0x02, 0x00, 0x20, 0x4a, 0x02, 0x00, 0x40, 0x49, 0x02, 0x00, 0x00,
            0x00, 0x02, 0x00, 0x40, 0x9b, 0x01, 0x00, 0x40, 0x9b, 0x01, 0x00, 0xfa, 0xf3, 0xaf,
            0x2a, 0xc0, 0x01,
        ];
        assert_eq!(decode(&stream, NX_TEXT.len()).unwrap(), NX_TEXT);
    }

    #[test]
    fn order1_shift_other_than_10_or_12_is_malformed() {
        // Header 0x30 gives shift 3.
        let bytes = [0x01, 3, 0x30, 0, 0x41, 0];
        assert!(matches!(decode(&bytes, 3), Err(SeqError::Malformed { .. })));
    }

    /// The 420-byte text that the RLE and stripe vectors below decode to.
    fn runs_text() -> Vec<u8> {
        b"AAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCGGGGGGGGGGGGGGTTTTTTTTTTTTACGTAAAAAAAA".repeat(6)
    }

    #[test]
    fn rle_stream_from_htscodecs() {
        // Flags 0x40 (RLE, no CAT), from htscodecs rans_compress_to_4x16 with order 0x40.
        let stream = [
            0x40, 0x83, 0x24, 0x6c, 0x31, 0x33, 0x00, 0x04, 0x07, 0x0b, 0x0d, 0x17, 0x1f, 0x41,
            0x43, 0x47, 0x54, 0x00, 0x1f, 0x01, 0x08, 0x07, 0x07, 0x01, 0x05, 0x01, 0x01, 0x01,
            0x01, 0xda, 0xa7, 0x15, 0x00, 0x10, 0x6f, 0x6b, 0x00, 0x6d, 0x0f, 0x02, 0x24, 0xb6,
            0x5f, 0xe7, 0x56, 0xe9, 0xaf, 0xb1, 0x36, 0x80, 0x76, 0x00, 0x25, 0x00, 0xe3, 0x40,
            0x63, 0x41, 0x43, 0x47, 0x54, 0x00, 0x13, 0x0f, 0x0f, 0x0f, 0xa2, 0x71, 0x34, 0x00,
            0x46, 0x36, 0x16, 0x01, 0xb3, 0x08, 0x1d, 0x01, 0xb2, 0x2f, 0x23, 0x01, 0x00, 0x83,
            0x40, 0x86, 0x40, 0x1b, 0xc0, 0xec,
        ];
        assert_eq!(decode(&stream, 420).unwrap(), runs_text());
    }

    #[test]
    fn stripe_stream_from_htscodecs() {
        // Flags 0x08 (STRIPE, four stripes, each a complete order-0 stream), from
        // htscodecs rans_compress_to_4x16 with order 0x08.
        let stream = [
            0x08, 0x83, 0x24, 0x04, 0x2a, 0x30, 0x2c, 0x30, 0x10, 0x41, 0x43, 0x47, 0x54, 0x00,
            0x3f, 0x0e, 0x1d, 0x16, 0x23, 0xb1, 0x18, 0x47, 0xa0, 0x12, 0xfa, 0x51, 0xd3, 0x04,
            0x9e, 0x41, 0x85, 0x11, 0xf0, 0x1c, 0xc2, 0xe7, 0x8f, 0x40, 0x84, 0x71, 0x18, 0xa3,
            0xa0, 0xd2, 0x40, 0xe2, 0x80, 0x00, 0x20, 0xf9, 0x10, 0x41, 0x43, 0x47, 0x54, 0x00,
            0x3c, 0x12, 0x19, 0x19, 0x76, 0x10, 0x01, 0x00, 0xfa, 0x94, 0xc3, 0x54, 0x76, 0x22,
            0x01, 0x00, 0x66, 0x80, 0x00, 0x00, 0x82, 0x75, 0x8c, 0x79, 0x00, 0xbb, 0x00, 0x04,
            0xa8, 0x14, 0xe2, 0x97, 0x22, 0xc8, 0xa0, 0xb0, 0xc0, 0xb5, 0x20, 0xc8, 0xe0, 0x89,
            0x10, 0x41, 0x43, 0x47, 0x54, 0x00, 0x3f, 0x0e, 0x1d, 0x16, 0x02, 0x55, 0xee, 0x1b,
            0x37, 0x44, 0xcf, 0x24, 0x7f, 0x81, 0x00, 0x00, 0xe9, 0x56, 0x43, 0x58, 0xed, 0x48,
            0xd7, 0x8a, 0x8a, 0x04, 0x05, 0x90, 0x49, 0x19, 0xa0, 0x71, 0xc0, 0x06, 0x20, 0x15,
            0xe0, 0xa8, 0x10, 0x41, 0x43, 0x47, 0x54, 0x00, 0x3c, 0x12, 0x19, 0x19, 0xe8, 0xc4,
            0x00, 0x00, 0x5c, 0x82, 0x00, 0x00, 0x95, 0xa4, 0x01, 0x00, 0xa4, 0x84, 0x61, 0x61,
            0x5c, 0x12, 0x70, 0x14, 0xd3, 0xec, 0x7a, 0x44, 0x70, 0x14, 0x8e, 0xf9, 0x77, 0x2c,
            0x20, 0x05, 0x40, 0xb7, 0x40, 0x09, 0x20, 0x2c,
        ];
        assert_eq!(decode(&stream, 420).unwrap(), runs_text());
    }

    #[test]
    fn truncated_and_corrupted_streams_never_panic() {
        let data = sample(300, 9);
        let enc = encode(&data);
        assert_eq!(decode(&enc, data.len()).unwrap(), data);
        // Every proper prefix is an error: the decoder consumes every byte.
        for cut in 0..enc.len() {
            assert!(decode(&enc[..cut], data.len()).is_err(), "truncated to {cut}");
        }
        // Flipping any byte must not panic, whatever it returns.
        for pos in 0..enc.len() {
            let mut bad = enc.clone();
            bad[pos] ^= 0xFF;
            let _ = decode(&bad, data.len());
        }
    }

    #[test]
    fn random_garbage_never_panics() {
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        for _ in 0..2000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let len = (s % 48) as usize;
            let bytes: Vec<u8> = (0..len)
                .map(|_| {
                    s ^= s << 13;
                    s ^= s >> 7;
                    s ^= s << 17;
                    (s >> 24) as u8
                })
                .collect();
            let _ = decode(&bytes, 16);
            let _ = decode(&bytes, 0);
        }
    }

    #[test]
    fn size_mismatch_is_malformed() {
        let enc = encode(b"ACGT");
        assert!(matches!(decode(&enc, 5), Err(SeqError::Malformed { .. })));
    }

    #[test]
    fn bad_table_sum_is_malformed() {
        // The table for "AA" stores 4096 as a0 00. Changing 00 to 10 gives 4112.
        let mut enc = encode(b"AA");
        enc[5] = 0x10;
        assert!(matches!(decode(&enc, 2), Err(SeqError::Malformed { .. })));
    }
}
