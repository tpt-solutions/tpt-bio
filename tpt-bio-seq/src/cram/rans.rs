//! rANS 4x8 entropy codec (CRAM codec method 4).
//!
//! Container layout: one order byte (0 or 1), the compressed size and the
//! uncompressed size as little-endian `u32`, then the body. The body holds a
//! frequency table, four little-endian `u32` initial states and the renormalized
//! byte stream. Only order 0 is implemented. Order 1 returns
//! [`SeqError::Unsupported`].
//!
//! The frequency-table reader follows the htscodecs run-length convention
//! (a symbol followed by `symbol + 1` is followed by a run count). The encoder
//! in this module writes the same convention, so the codec round-trips. It has
//! not been checked against samtools-written CRAM files.

use crate::cram::bytes::Cursor;
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM rANS";
const TF_SHIFT: u32 = 12;
const TOTAL: u32 = 1 << TF_SHIFT;
const LOWER: u32 = 1 << 23;

/// Decodes a rANS 4x8 block payload.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for order 1, and [`SeqError::Malformed`]
/// for inconsistent sizes or a frequency table that does not sum to 4096.
pub fn decode(input: &[u8]) -> Result<Vec<u8>> {
    let mut c = Cursor::new(input);
    let order = c.u8()?;
    let comp_len = c.u32_le()? as usize;
    let raw_len = c.u32_le()? as usize;
    if comp_len > input.len().saturating_sub(9) {
        return Err(SeqError::malformed(FORMAT, "compressed size exceeds input"));
    }
    match order {
        0 => decode_o0(&input[9..9 + comp_len], raw_len),
        1 => Err(SeqError::Unsupported("rANS 4x8 order-1 codec".into())),
        _ => Err(SeqError::malformed(FORMAT, format!("unknown order {order}"))),
    }
}

fn read_freq(c: &mut Cursor<'_>) -> Result<u32> {
    let b = u32::from(c.u8()?);
    if b >= 128 { Ok(((b & 0x7F) << 8) | u32::from(c.u8()?)) } else { Ok(b) }
}

/// Reads the order-0 frequency table. Returns frequencies indexed by symbol.
fn read_table(c: &mut Cursor<'_>, allow_empty: bool) -> Result<[u32; 256]> {
    let mut freq = [0u32; 256];
    let mut sym = c.u8()?;
    let mut rle = 0u32;
    loop {
        freq[usize::from(sym)] = read_freq(c)?;
        if rle > 0 {
            rle -= 1;
            sym = sym.wrapping_add(1);
        } else {
            let next = c.u8()?;
            if u16::from(next) == u16::from(sym) + 1 {
                rle = u32::from(c.u8()?);
            }
            sym = next;
            if sym == 0 {
                break;
            }
        }
    }
    let sum: u32 = freq.iter().sum();
    if sum != TOTAL && !(sum == 0 && allow_empty) {
        return Err(SeqError::malformed(FORMAT, format!("frequencies sum to {sum}, not 4096")));
    }
    Ok(freq)
}

fn decode_o0(body: &[u8], raw_len: usize) -> Result<Vec<u8>> {
    let mut cur = Cursor::new(body);
    let freq = read_table(&mut cur, raw_len == 0)?;
    let mut cum = [0u32; 257];
    for s in 0..256 {
        cum[s + 1] = cum[s] + freq[s];
    }
    // Slot-to-symbol lookup over the 4096 cumulative slots.
    let mut lookup = vec![0u8; TOTAL as usize];
    for s in 0..256 {
        for slot in cum[s]..cum[s + 1] {
            lookup[slot as usize] = s as u8;
        }
    }
    let mut state = [0u32; 4];
    for st in &mut state {
        *st = cur.u32_le()?;
    }
    let mut out = Vec::with_capacity(raw_len);
    for i in 0..raw_len {
        let j = i & 3;
        let x = state[j];
        let f = x & (TOTAL - 1);
        let s = lookup[f as usize];
        out.push(s);
        let s = usize::from(s);
        let mut x = freq[s] * (x >> TF_SHIFT) + f - cum[s];
        while x < LOWER {
            x = (x << 8) | u32::from(cur.u8()?);
        }
        state[j] = x;
    }
    Ok(out)
}

/// Normalizes raw counts so that every nonzero count is at least 1 and the
/// total is exactly 4096.
fn normalize(counts: &[u64; 256]) -> [u32; 256] {
    let total: u64 = counts.iter().sum();
    let mut freq = [0u32; 256];
    if total == 0 {
        return freq;
    }
    let mut sum = 0i64;
    for s in 0..256 {
        if counts[s] > 0 {
            let f = (counts[s] * u64::from(TOTAL) / total).max(1) as u32;
            freq[s] = f;
            sum += i64::from(f);
        }
    }
    // Put the rounding difference on the most frequent symbol.
    let top = (0..256).max_by_key(|&s| counts[s]).unwrap_or(0);
    freq[top] = (i64::from(freq[top]) + i64::from(TOTAL) - sum) as u32;
    freq
}

fn write_table(freq: &[u32; 256], out: &mut Vec<u8>) {
    let syms: Vec<usize> = (0..256).filter(|&s| freq[s] > 0).collect();
    if syms.is_empty() {
        // Empty input: one zero-frequency entry for symbol 0, then the terminator.
        out.extend_from_slice(&[0, 0, 0]);
        return;
    }
    for (idx, &s) in syms.iter().enumerate() {
        if idx == 0 {
            out.push(s as u8);
        }
        let f = freq[s];
        if f < 128 {
            out.push(f as u8);
        } else {
            out.push(0x80 | (f >> 8) as u8);
            out.push(f as u8);
        }
        match syms.get(idx + 1) {
            Some(&next) if next == s + 1 => {
                out.push(next as u8);
                out.push(0); // run count of zero: the next symbol is read normally
            }
            Some(&next) => out.push(next as u8),
            None => out.push(0), // terminator
        }
    }
}

/// Encodes `data` with rANS 4x8 order 0. Used by tests and by the writer.
#[must_use]
pub fn encode_o0(data: &[u8]) -> Vec<u8> {
    let mut counts = [0u64; 256];
    for &b in data {
        counts[usize::from(b)] += 1;
    }
    let freq = normalize(&counts);
    let mut cum = [0u32; 257];
    for s in 0..256 {
        cum[s + 1] = cum[s] + freq[s];
    }
    let mut body = Vec::new();
    write_table(&freq, &mut body);
    let mut state = [LOWER; 4];
    let mut stream = Vec::new();
    for i in (0..data.len()).rev() {
        let j = i & 3;
        let s = usize::from(data[i]);
        let f = freq[s];
        let x_max = ((LOWER >> TF_SHIFT) << 8) * f;
        let mut x = state[j];
        while x >= x_max {
            stream.push(x as u8);
            x >>= 8;
        }
        state[j] = ((x / f) << TF_SHIFT) + (x % f) + cum[s];
    }
    // Flush the states so that after reversal they appear as LE u32 in order 0..4.
    for j in (0..4).rev() {
        stream.extend_from_slice(&state[j].to_be_bytes());
    }
    stream.reverse();
    body.extend_from_slice(&stream);

    let mut out = Vec::with_capacity(body.len() + 9);
    out.push(0);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out.extend_from_slice(&body);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order0_roundtrip_various() {
        let cases: Vec<Vec<u8>> = vec![
            b"ACGTACGTTTGA".to_vec(),
            vec![b'A'; 1000],
            (0..=255u8).cycle().take(5000).collect(),
            b"x".to_vec(),
            Vec::new(),
        ];
        for data in cases {
            let enc = encode_o0(&data);
            let dec = decode(&enc).unwrap();
            let first_bad = dec.iter().zip(&data).position(|(a, b)| a != b);
            assert_eq!(first_bad, None, "len {} first mismatch", data.len());
            assert_eq!(dec.len(), data.len());
        }
    }

    #[test]
    fn order1_is_unsupported() {
        let mut enc = encode_o0(b"ACGT");
        enc[0] = 1;
        assert!(matches!(decode(&enc), Err(SeqError::Unsupported(_))));
    }

    #[test]
    fn bad_size_is_malformed() {
        let mut enc = encode_o0(b"ACGT");
        enc[1..5].copy_from_slice(&1_000_000u32.to_le_bytes());
        assert!(matches!(decode(&enc), Err(SeqError::Malformed { .. })));
    }
}
