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
//! | 0x04 | `N32`     | [`SeqError::Unsupported`] (32-state coder) |
//! | 0x01 | `ORDER`   | [`SeqError::Unsupported`] (order-1 coder) |
//!
//! The order-0 body is a run-length symbol list with a frequency for each symbol,
//! then four 32-bit initial states, then 16-bit renormalization words. Symbol `i` is
//! coded by state `i % 4`. Frequencies are normalized by powers of two up to 4096.
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
const TOTAL: u32 = 1 << TF_SHIFT;
/// Renormalization bound for 16-bit words.
const LOWER16: u32 = 1 << 15;
const N_STATES: usize = 4;
/// Largest output accepted, as a guard against corrupt sizes.
const MAX_OUTPUT: usize = 1 << 28;
/// Stripes and RLE and pack maps are capped, as a guard against corrupt counts.
const MAX_STRIPES: usize = 255;

/// Decodes an Nx16 payload. `expected_len` is the block's raw size. It is used
/// for `NO_SIZE` streams and must match the declared size otherwise.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for the 32-state and order-1 coders, and
/// [`SeqError::Malformed`] for truncated input, size mismatches, or tables that
/// do not normalize to 4096.
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
    if flags & FLAG_N32 != 0 {
        return Err(SeqError::Unsupported("rANS Nx16 with 32 states (N32)".into()));
    }
    if flags & FLAG_ORDER != 0 {
        return Err(SeqError::Unsupported("rANS Nx16 order 1".into()));
    }
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
        decode_o0_body(c, stage1_len)?
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
        decode_o0_body(&mut Cursor::new(body), u_meta / 2)?
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

/// Decodes an order-0 body of `out_len` bytes: table, initial states, then the data.
fn decode_o0_body(c: &mut Cursor<'_>, out_len: usize) -> Result<Vec<u8>> {
    if out_len == 0 {
        return Ok(Vec::new());
    }
    let freq = read_normalized_freq(c)?;
    let table = Table::new(&freq);
    let mut state = [0u32; N_STATES];
    for s in &mut state {
        *s = c.u32_le()?;
    }
    let mut out = Vec::with_capacity(out_len.min(1 << 20));
    for i in 0..out_len {
        let lane = i % N_STATES;
        let slot = state[lane] & (TOTAL - 1);
        let sym = table.symbol_at(slot);
        let s = usize::from(sym);
        let mut x = table.freq[s] * (state[lane] >> TF_SHIFT) + slot - table.cum[s];
        if x < LOWER16 {
            x = (x << 16) | read_word(c)?;
        }
        state[lane] = x;
        out.push(sym);
    }
    Ok(out)
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

    #[test]
    fn unsupported_coders_are_reported() {
        for flag in [FLAG_N32, FLAG_ORDER] {
            let err = decode(&[flag, 0], 0).unwrap_err();
            assert!(matches!(err, SeqError::Unsupported(_)), "flag {flag:#x}");
        }
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
