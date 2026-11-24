//! rANS Nx16 entropy codec (CRAM 3.1 codec method 5).
//!
//! Stream layout: a flags byte, then (unless `NO_SIZE`) the uncompressed size
//! as an unsigned LEB128 varint. With `CAT` the raw bytes follow directly.
//! Otherwise the body holds the frequency table(s), then the initial states as
//! little-endian `u32` (4, or 32 with `N32`), then the 16-bit renormalization
//! words. Symbol `i` is coded by state `i % N`.
//!
//! Flags (the bit values used by htscodecs):
//!
//! | Bit  | Name      | Support |
//! |------|-----------|---------|
//! | 0x01 | `ORDER`   | order 0 and order 1 (context = previous byte of the output) |
//! | 0x04 | `N32`     | 32 interleaved states |
//! | 0x08 | `STRIPE`  | [`SeqError::Unsupported`] |
//! | 0x10 | `NO_SIZE` | size taken from the caller's expected length |
//! | 0x20 | `CAT`     | stored, uncompressed |
//! | 0x40 | `RLE`     | [`SeqError::Unsupported`] |
//! | 0x80 | `PACK`    | [`SeqError::Unsupported`] |
//!
//! Frequency tables use the same run-length symbol list as [`crate::cram::rans`],
//! with each frequency as an unsigned LEB128 varint. Every non-empty table sums
//! to 4096. The encoder writes size and table in the same layout, so the codec
//! round-trips. The byte layout has not been checked against htslib output.

use crate::cram::bytes::Cursor;
use crate::cram::rans::{
    MAX_OUTPUT, TF_SHIFT, TOTAL, Table, normalize, read_table_with, walk_symbols, write_symbols,
};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM rANS Nx16";
const FLAG_ORDER: u8 = 0x01;
const FLAG_N32: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;
/// Renormalization bound for 16-bit words.
const LOWER16: u32 = 1 << 15;

/// Decodes an Nx16 payload. `expected_len` is the block's raw size. It is used
/// for `NO_SIZE` streams and must match the declared size otherwise.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for `STRIPE`, `RLE` and `PACK`, and
/// [`SeqError::Malformed`] for truncated input, size mismatches, or tables
/// that do not sum to 4096.
pub fn decode(input: &[u8], expected_len: usize) -> Result<Vec<u8>> {
    let mut c = Cursor::new(input);
    let flags = c.u8()?;
    for (bit, name) in [(FLAG_STRIPE, "STRIPE"), (FLAG_RLE, "RLE"), (FLAG_PACK, "PACK")] {
        if flags & bit != 0 {
            return Err(SeqError::Unsupported(format!("rANS Nx16 flag {name}")));
        }
    }
    let size = if flags & FLAG_NO_SIZE != 0 {
        expected_len
    } else {
        let declared = read_uleb(&mut c)?;
        if declared != expected_len {
            return Err(SeqError::malformed(
                FORMAT,
                format!("stream declares {declared} bytes, block declares {expected_len}"),
            ));
        }
        declared
    };
    if size > MAX_OUTPUT {
        return Err(SeqError::malformed(FORMAT, format!("declared size {size} is too large")));
    }
    if flags & FLAG_CAT != 0 {
        let raw = c.bytes(size)?;
        if !c.is_empty() {
            return Err(SeqError::malformed(FORMAT, "trailing bytes after stored data"));
        }
        return Ok(raw.to_vec());
    }
    let states = if flags & FLAG_N32 != 0 { 32 } else { 4 };
    if flags & FLAG_ORDER != 0 {
        decode_o1(&mut c, size, states)
    } else {
        decode_o0(&mut c, size, states)
    }
}

/// Encodes `data` as an Nx16 stream with the given order (0 or 1) and state count
/// (4, or 32 when `states32` is set).
#[must_use]
pub fn encode(data: &[u8], order1: bool, states32: bool) -> Vec<u8> {
    let n_states = if states32 { 32 } else { 4 };
    let mut out = Vec::new();
    let mut flags = 0u8;
    if order1 {
        flags |= FLAG_ORDER;
    }
    if states32 {
        flags |= FLAG_N32;
    }
    out.push(flags);
    put_uleb(&mut out, data.len());
    let mut body = Vec::new();
    let ctx_of = |k: usize| if order1 && k > 0 { usize::from(data[k - 1]) } else { 0 };
    let mut counts = vec![[0u64; 256]; if order1 { 256 } else { 1 }];
    for k in 0..data.len() {
        counts[ctx_of(k)][usize::from(data[k])] += 1;
    }
    let tables: Vec<Table> = counts.iter().map(|c| Table::new(&normalize(c))).collect();
    if order1 {
        let ctxs: Vec<u8> = (0..256usize)
            .filter(|&c| c == 0 || counts[c].iter().any(|&v| v > 0))
            .map(|c| c as u8)
            .collect();
        write_symbols(&ctxs, &mut body, |out, c| write_table(&tables[usize::from(c)], out));
    } else {
        write_table(&tables[0], &mut body);
    }
    let mut state = vec![LOWER16; n_states];
    let mut emitted: Vec<u16> = Vec::new();
    for i in (0..data.len()).rev() {
        let t = &tables[ctx_of(i)];
        let s = usize::from(data[i]);
        let f = t.freq[s];
        let x = &mut state[i % n_states];
        let x_max = ((LOWER16 >> TF_SHIFT) << 16) * f;
        while *x >= x_max {
            emitted.push(*x as u16);
            *x >>= 16;
        }
        *x = ((*x / f) << TF_SHIFT) + (*x % f) + t.cum[s];
    }
    for x in &state {
        body.extend_from_slice(&x.to_le_bytes());
    }
    for w in emitted.iter().rev() {
        body.extend_from_slice(&w.to_le_bytes());
    }
    out.extend_from_slice(&body);
    out
}

fn write_table(t: &Table, out: &mut Vec<u8>) {
    let syms: Vec<u8> = (0..256usize).filter(|&s| t.freq[s] > 0).map(|s| s as u8).collect();
    write_symbols(&syms, out, |out, s| put_uleb(out, t.freq[usize::from(s)] as usize));
}

fn read_uleb(c: &mut Cursor<'_>) -> Result<usize> {
    let mut value = 0usize;
    for shift in (0..64).step_by(7) {
        let b = c.u8()?;
        value |= usize::from(b & 0x7F) << shift;
        if b & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(SeqError::malformed(FORMAT, "varint too long"))
}

fn put_uleb(out: &mut Vec<u8>, mut v: usize) {
    loop {
        let b = (v & 0x7F) as u8;
        v >>= 7;
        if v == 0 {
            out.push(b);
            return;
        }
        out.push(b | 0x80);
    }
}

fn read_word(c: &mut Cursor<'_>) -> Result<u32> {
    let b = c.bytes(2)?;
    Ok(u32::from(u16::from_le_bytes([b[0], b[1]])))
}

/// One Nx16 decode step for slot `f` of `t`, with 16-bit renormalization.
fn step(cur: &mut Cursor<'_>, table: &Table, state: u32, slot: u32) -> Result<(u8, u32)> {
    let sym = table.symbol_at(slot);
    let s = usize::from(sym);
    let mut x = table.freq[s] * (state >> TF_SHIFT) + slot - table.cum[s];
    while x < LOWER16 {
        x = (x << 16) | read_word(cur)?;
    }
    Ok((sym, x))
}

fn read_states(c: &mut Cursor<'_>, n: usize) -> Result<Vec<u32>> {
    (0..n).map(|_| c.u32_le()).collect()
}

fn decode_o0(c: &mut Cursor<'_>, size: usize, n_states: usize) -> Result<Vec<u8>> {
    let table = read_table_with(c, size == 0, read_uleb_freq)?;
    let mut state = read_states(c, n_states)?;
    let mut out = Vec::with_capacity(size.min(1 << 20));
    for i in 0..size {
        let j = i % n_states;
        let (s, x) = step(c, &table, state[j], state[j] & (TOTAL - 1))?;
        out.push(s);
        state[j] = x;
    }
    Ok(out)
}

fn decode_o1(c: &mut Cursor<'_>, size: usize, n_states: usize) -> Result<Vec<u8>> {
    let mut tables: Vec<Option<Box<Table>>> = (0..256).map(|_| None).collect();
    walk_symbols(c, |c, ctx| {
        let t = read_table_with(c, size == 0, read_uleb_freq)?;
        tables[usize::from(ctx)] = Some(Box::new(t));
        Ok(())
    })?;
    let mut state = read_states(c, n_states)?;
    let mut out: Vec<u8> = Vec::with_capacity(size.min(1 << 20));
    for i in 0..size {
        let lane = i % n_states;
        let ctx = if i == 0 { 0 } else { usize::from(out[i - 1]) };
        let table = tables[ctx]
            .as_deref()
            .ok_or_else(|| SeqError::malformed(FORMAT, "missing order-1 context table"))?;
        let (sym, x) = step(c, table, state[lane], state[lane] & (TOTAL - 1))?;
        out.push(sym);
        state[lane] = x;
    }
    Ok(out)
}

fn read_uleb_freq(c: &mut Cursor<'_>) -> Result<u32> {
    u32::try_from(read_uleb(c)?).map_err(|_| SeqError::malformed(FORMAT, "frequency too large"))
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
    fn roundtrip_all_order_and_state_combinations() {
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
            for order1 in [false, true] {
                for states32 in [false, true] {
                    let enc = encode(data, order1, states32);
                    let dec = decode(&enc, data.len()).unwrap();
                    assert_eq!(&dec, data, "len {} order1 {order1} n32 {states32}", data.len());
                }
            }
        }
    }

    #[test]
    fn hand_derived_order0_bytes() {
        // "AA" with 4 states: flags 0, size 2, then the table with A = 4096.
        // 4096 as a LEB128 varint is [0x80, 0x20]; then the terminator 0.
        // Frequency 4096 leaves every state unchanged at LOWER16 = 0x8000.
        let enc = encode(b"AA", false, false);
        let mut expected = vec![0x00, 0x02, 65, 0x80, 0x20, 0x00];
        for _ in 0..4 {
            expected.extend_from_slice(&[0x00, 0x80, 0x00, 0x00]);
        }
        assert_eq!(enc, expected);
        assert_eq!(decode(&enc, 2).unwrap(), b"AA");
    }

    #[test]
    fn cat_and_no_size_flags() {
        let raw = b"stored".to_vec();
        let mut cat = vec![FLAG_CAT, 6];
        cat.extend_from_slice(&raw);
        assert_eq!(decode(&cat, 6).unwrap(), raw);
        let mut no_size = encode(b"ACGT", false, false);
        no_size[0] |= FLAG_NO_SIZE;
        no_size.remove(1);
        assert_eq!(decode(&no_size, 4).unwrap(), b"ACGT");
    }

    #[test]
    fn unsupported_flags_are_reported() {
        for flag in [FLAG_STRIPE, FLAG_RLE, FLAG_PACK] {
            let err = decode(&[flag, 0], 0).unwrap_err();
            assert!(matches!(err, SeqError::Unsupported(_)), "flag {flag:#x}");
        }
    }

    #[test]
    fn truncated_and_corrupted_streams_never_panic() {
        use crate::cram::rans;
        let data = sample(300, 9);
        let streams: Vec<(&str, Vec<u8>)> = vec![
            ("rans o0", rans::encode_o0(&data)),
            ("rans o1", rans::encode_o1(&data)),
            ("nx16 o0", encode(&data, false, false)),
            ("nx16 o1", encode(&data, true, false)),
            ("nx16 o1 n32", encode(&data, true, true)),
        ];
        for (name, enc) in &streams {
            let decode_any = |bytes: &[u8]| {
                if name.starts_with("rans") {
                    rans::decode(bytes)
                } else {
                    decode(bytes, data.len())
                }
            };
            assert_eq!(decode_any(enc).unwrap(), data, "{name} baseline");
            // Every proper prefix is an error: the decoder consumes every byte.
            for cut in 0..enc.len() {
                assert!(decode_any(&enc[..cut]).is_err(), "{name} truncated to {cut}");
            }
            // Flipping any byte must not panic, whatever it returns.
            for pos in 0..enc.len() {
                let mut bad = enc.clone();
                bad[pos] ^= 0xFF;
                let _ = decode_any(&bad);
            }
        }
    }

    #[test]
    fn random_garbage_never_panics() {
        use crate::cram::rans;
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
            let _ = rans::decode(&bytes);
            let _ = decode(&bytes, 16);
            let _ = decode(&bytes, 0);
        }
    }

    #[test]
    fn size_mismatch_is_malformed() {
        let enc = encode(b"ACGT", false, false);
        assert!(matches!(decode(&enc, 5), Err(SeqError::Malformed { .. })));
    }

    #[test]
    fn bad_table_sum_is_malformed() {
        // The table for "AA" stores 4096 as varint bytes 0x80 0x20 at [3..5].
        // Changing 0x20 to 0x10 gives 2048, which is not 4096.
        let mut enc = encode(b"AA", false, false);
        enc[4] = 0x10;
        assert!(matches!(decode(&enc, 2), Err(SeqError::Malformed { .. })));
    }
}
