//! rANS 4x8 entropy codec (CRAM codec method 4), orders 0 and 1.
//!
//! Container layout: one order byte (0 or 1), the compressed size and the
//! uncompressed size as little-endian `u32`, then the body. The body holds the
//! frequency table(s), four little-endian `u32` initial states and the
//! renormalized byte stream.
//!
//! Order 0 uses one 256-entry frequency table. Order 1 uses one table per
//! context, where the context of a byte is the previous byte of the same
//! interleaved stream (0 at the start of each stream). The four interleaved
//! streams split the data into quarters; stream 3 also takes the remainder.
//!
//! Frequency tables use the htscodecs run-length convention: a symbol followed
//! by `symbol + 1` is followed by a run count of further consecutive symbols,
//! and a zero byte ends the list. Frequencies below 128 take one byte, others
//! two (`0x80 | high`, `low`). Every non-empty table sums to 4096.
//!
//! The encoders in this module write the same conventions, so they round-trip.
//! Round-trip and hand-derived byte tests cover the layout. The byte layout has
//! not been checked against htscodecs or samtools output.

use crate::cram::bytes::Cursor;
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM rANS";
pub(crate) const TF_SHIFT: u32 = 12;
pub(crate) const TOTAL: u32 = 1 << TF_SHIFT;
pub(crate) const LOWER: u32 = 1 << 23;
/// Largest output a single block may declare. A one-symbol table consumes no
/// input bytes, so the output size cannot be bounded by the input size.
pub(crate) const MAX_OUTPUT: usize = 1 << 28;

/// Decodes a rANS 4x8 block payload (order 0 or 1).
///
/// # Errors
/// Returns [`SeqError::Malformed`] for unknown orders, inconsistent sizes,
/// truncated input, or a frequency table that does not sum to 4096.
pub fn decode(input: &[u8]) -> Result<Vec<u8>> {
    let mut c = Cursor::new(input);
    let order = c.u8()?;
    let comp_len = c.u32_le()? as usize;
    let raw_len = c.u32_le()? as usize;
    if comp_len > input.len().saturating_sub(9) {
        return Err(SeqError::malformed(FORMAT, "compressed size exceeds input"));
    }
    if raw_len > MAX_OUTPUT {
        return Err(SeqError::malformed(FORMAT, format!("declared size {raw_len} is too large")));
    }
    let body = &input[9..9 + comp_len];
    match order {
        0 => decode_o0(body, raw_len),
        1 => decode_o1(body, raw_len),
        _ => Err(SeqError::malformed(FORMAT, format!("unknown order {order}"))),
    }
}

/// A frequency table with its cumulative offsets and a slot-to-symbol lookup.
#[derive(Debug, Clone)]
pub(crate) struct Table {
    pub(crate) freq: [u32; 256],
    pub(crate) cum: [u32; 257],
    lookup: Vec<u8>,
}

impl Table {
    /// Builds the lookup for a table whose frequencies sum to 4096 or are all zero.
    pub(crate) fn new(freq: &[u32; 256]) -> Self {
        let freq = *freq;
        let mut cum = [0u32; 257];
        for s in 0..256 {
            cum[s + 1] = cum[s] + freq[s];
        }
        let mut lookup = vec![0u8; TOTAL as usize];
        for s in 0..256 {
            for slot in cum[s]..cum[s + 1] {
                lookup[slot as usize] = s as u8;
            }
        }
        Self { freq, cum, lookup }
    }

    pub(crate) fn symbol_at(&self, slot: u32) -> u8 {
        self.lookup[slot as usize]
    }
}

/// Walks a run-length coded symbol list (the layout shared by frequency tables
/// and the order-1 context list). `each` is called once per listed symbol, in
/// increasing order, and may read that symbol's payload from the cursor.
pub(crate) fn walk_symbols(
    c: &mut Cursor<'_>,
    mut each: impl FnMut(&mut Cursor<'_>, u8) -> Result<()>,
) -> Result<()> {
    let mut sym = c.u8()?;
    let mut run = 0u32;
    loop {
        each(c, sym)?;
        if run > 0 {
            if sym == 255 {
                return Err(SeqError::malformed(FORMAT, "symbol run passes 255"));
            }
            run -= 1;
            sym += 1;
        } else {
            let next = c.u8()?;
            if u16::from(next) == u16::from(sym) + 1 {
                run = u32::from(c.u8()?);
            }
            sym = next;
            if sym == 0 {
                return Ok(());
            }
        }
    }
}

/// Writes a run-length coded symbol list, calling `body` for each symbol
/// after its header byte. The list ends with a zero byte. An empty list is
/// written as a single zero-valued entry (its body is written for symbol 0)
/// followed by the terminator, which is how the reader expects an empty table.
pub(crate) fn write_symbols(
    syms: &[u8],
    out: &mut Vec<u8>,
    mut body: impl FnMut(&mut Vec<u8>, u8),
) {
    if syms.is_empty() {
        out.push(0);
        body(out, 0);
        out.push(0);
        return;
    }
    for (idx, &s) in syms.iter().enumerate() {
        if idx == 0 {
            out.push(s);
        }
        body(out, s);
        match syms.get(idx + 1) {
            Some(&next) if next == s + 1 => {
                out.push(next);
                out.push(0); // run count of zero: the next symbol follows as a header
            }
            Some(&next) => out.push(next),
            None => out.push(0),
        }
    }
}

fn write_freq(f: u32, out: &mut Vec<u8>) {
    if f < 128 {
        out.push(f as u8);
    } else {
        out.push(0x80 | (f >> 8) as u8);
        out.push(f as u8);
    }
}

fn read_freq(c: &mut Cursor<'_>) -> Result<u32> {
    let b = u32::from(c.u8()?);
    if b >= 128 { Ok(((b & 0x7F) << 8) | u32::from(c.u8()?)) } else { Ok(b) }
}

/// Reads one 4x8 table. Empty tables are accepted only when `allow_empty` is
/// set; any other total must be exactly 4096.
pub(crate) fn read_table(c: &mut Cursor<'_>, allow_empty: bool) -> Result<Table> {
    read_table_with(c, allow_empty, read_freq)
}

/// Reads one table whose frequencies are decoded by `read_one`.
pub(crate) fn read_table_with(
    c: &mut Cursor<'_>,
    allow_empty: bool,
    mut read_one: impl FnMut(&mut Cursor<'_>) -> Result<u32>,
) -> Result<Table> {
    let mut freq = [0u32; 256];
    walk_symbols(c, |c, s| {
        freq[usize::from(s)] = read_one(c)?;
        Ok(())
    })?;
    let sum: u32 = freq.iter().sum();
    if sum != TOTAL && !(sum == 0 && allow_empty) {
        return Err(SeqError::malformed(FORMAT, format!("frequencies sum to {sum}, not 4096")));
    }
    Ok(Table::new(&freq))
}

fn decode_o0(body: &[u8], raw_len: usize) -> Result<Vec<u8>> {
    let mut cur = Cursor::new(body);
    let table = read_table(&mut cur, raw_len == 0)?;
    let mut state = [0u32; 4];
    for st in &mut state {
        *st = cur.u32_le()?;
    }
    let mut out = Vec::with_capacity(raw_len.min(1 << 20));
    for i in 0..raw_len {
        let j = i & 3;
        let x = state[j];
        let f = x & (TOTAL - 1);
        let s = table.symbol_at(f);
        out.push(s);
        state[j] = step(&mut cur, &table, s, x, f)?;
    }
    Ok(out)
}

/// One rANS decode step for symbol `s`, followed by byte renormalization.
fn step(cur: &mut Cursor<'_>, table: &Table, s: u8, x: u32, f: u32) -> Result<u32> {
    let s = usize::from(s);
    let mut x = table.freq[s] * (x >> TF_SHIFT) + f - table.cum[s];
    while x < LOWER {
        x = (x << 8) | u32::from(cur.u8()?);
    }
    Ok(x)
}

fn decode_o1(body: &[u8], raw_len: usize) -> Result<Vec<u8>> {
    let mut cur = Cursor::new(body);
    let mut tables: Vec<Option<Box<Table>>> = (0..256).map(|_| None).collect();
    walk_symbols(&mut cur, |c, ctx| {
        let table = read_table(c, raw_len == 0)?;
        tables[usize::from(ctx)] = Some(Box::new(table));
        Ok(())
    })?;
    let mut state = [0u32; 4];
    for st in &mut state {
        *st = cur.u32_le()?;
    }
    let isz4 = raw_len / 4;
    let starts = [0, isz4, 2 * isz4, 3 * isz4];
    let mut out = vec![0u8; raw_len];
    let mut ctx = [0u8; 4];
    for t in 0..isz4 {
        for (j, &start) in starts.iter().enumerate() {
            let idx = start + t;
            decode_o1_at(&mut cur, &tables, &mut ctx, &mut state, &mut out, j, idx)?;
        }
    }
    for idx in starts[3] + isz4..raw_len {
        decode_o1_at(&mut cur, &tables, &mut ctx, &mut state, &mut out, 3, idx)?;
    }
    Ok(out)
}

/// Decodes one order-1 symbol for stream `j` at output position `idx`.
fn decode_o1_at(
    cur: &mut Cursor<'_>,
    tables: &[Option<Box<Table>>],
    ctx: &mut [u8; 4],
    state: &mut [u32; 4],
    out: &mut [u8],
    j: usize,
    idx: usize,
) -> Result<()> {
    let table = tables[usize::from(ctx[j])]
        .as_deref()
        .ok_or_else(|| SeqError::malformed(FORMAT, "missing order-1 context table"))?;
    let x = state[j];
    let f = x & (TOTAL - 1);
    let s = table.symbol_at(f);
    out[idx] = s;
    state[j] = step(cur, table, s, x, f)?;
    ctx[j] = s;
    Ok(())
}

/// Normalizes raw counts so that every nonzero count is at least 1 and the
/// total is exactly 4096.
pub(crate) fn normalize(counts: &[u64; 256]) -> [u32; 256] {
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
    // Shave units from the largest entries if the minimum-of-one rule overshot.
    while sum > i64::from(TOTAL) {
        let s = (0..256).max_by_key(|&s| freq[s]).unwrap_or(0);
        freq[s] -= 1;
        sum -= 1;
    }
    // Put the rounding difference on the most frequent symbol.
    let top = (0..256).max_by_key(|&s| counts[s]).unwrap_or(0);
    freq[top] = (i64::from(freq[top]) + i64::from(TOTAL) - sum) as u32;
    freq
}

fn write_table(freq: &[u32; 256], out: &mut Vec<u8>) {
    let syms: Vec<u8> = (0..256usize).filter(|&s| freq[s] > 0).map(|s| s as u8).collect();
    write_symbols(&syms, out, |out, s| write_freq(freq[usize::from(s)], out));
}

/// Frames a body as a rANS 4x8 payload with the given order.
fn frame(order: u8, body: &[u8], raw_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(body.len() + 9);
    out.push(order);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(&(raw_len as u32).to_le_bytes());
    out.extend_from_slice(body);
    out
}

/// Flushes the four states (as the decoder reads them) and reverses the
/// emitted byte stream so that the decoder reads it forward.
fn finish(state: [u32; 4], emitted: &mut Vec<u8>, body: &mut Vec<u8>) {
    for j in (0..4).rev() {
        emitted.extend_from_slice(&state[j].to_be_bytes());
    }
    emitted.reverse();
    body.extend_from_slice(emitted);
}

/// Encodes one symbol into a state, emitting renormalization bytes.
fn put(emitted: &mut Vec<u8>, table: &Table, s: u8, x: &mut u32) {
    let s = usize::from(s);
    let f = table.freq[s];
    let x_max = ((LOWER >> TF_SHIFT) << 8) * f;
    while *x >= x_max {
        emitted.push(*x as u8);
        *x >>= 8;
    }
    *x = ((*x / f) << TF_SHIFT) + (*x % f) + table.cum[s];
}

/// Encodes `data` with rANS 4x8 order 0.
#[must_use]
pub fn encode_o0(data: &[u8]) -> Vec<u8> {
    let mut counts = [0u64; 256];
    for &b in data {
        counts[usize::from(b)] += 1;
    }
    let table = Table::new(&normalize(&counts));
    let mut body = Vec::new();
    write_table(&table.freq, &mut body);
    let mut state = [LOWER; 4];
    let mut emitted = Vec::new();
    for i in (0..data.len()).rev() {
        put(&mut emitted, &table, data[i], &mut state[i & 3]);
    }
    finish(state, &mut emitted, &mut body);
    frame(0, &body, data.len())
}

/// Encodes `data` with rANS 4x8 order 1.
#[must_use]
pub fn encode_o1(data: &[u8]) -> Vec<u8> {
    let n = data.len();
    let isz4 = n / 4;
    let starts = [0, isz4, 2 * isz4, 3 * isz4];
    // The context of byte `k` is the previous byte, or 0 at a stream start.
    let ctx_of = |k: usize| if starts.contains(&k) { 0 } else { usize::from(data[k - 1]) };
    let mut counts = vec![[0u64; 256]; 256];
    for k in 0..n {
        counts[ctx_of(k)][usize::from(data[k])] += 1;
    }
    // Context 0 is always listed: every stream starts with it.
    let ctxs: Vec<u8> = (0..256usize)
        .filter(|&c| c == 0 || counts[c].iter().any(|&v| v > 0))
        .map(|c| c as u8)
        .collect();
    let tables: Vec<Table> = counts.iter().map(|c| Table::new(&normalize(c))).collect();
    let mut body = Vec::new();
    write_symbols(&ctxs, &mut body, |out, c| write_table(&tables[usize::from(c)].freq, out));

    let mut state = [LOWER; 4];
    let mut emitted = Vec::new();
    // Mirror the decoder schedule in reverse: the stream-3 tail first, then
    // quarters from the last position to the first, states 3 down to 0.
    for idx in (starts[3] + isz4..n).rev() {
        put(&mut emitted, &tables[ctx_of(idx)], data[idx], &mut state[3]);
    }
    for t in (0..isz4).rev() {
        for j in (0..4).rev() {
            let idx = starts[j] + t;
            put(&mut emitted, &tables[ctx_of(idx)], data[idx], &mut state[j]);
        }
    }
    finish(state, &mut emitted, &mut body);
    frame(1, &body, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_inputs() -> Vec<Vec<u8>> {
        let mut seed = 0x2545_F491_4F6C_DD1Du64;
        let mut random = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    (seed >> 24) as u8
                })
                .collect()
        };
        vec![
            Vec::new(),
            b"x".to_vec(),
            b"ACG".to_vec(),
            random(1),
            random(3),
            random(4),
            random(7),
            random(1000),
            random(65_537),
            b"ACGTACGTTTGA".to_vec(),
            vec![b'A'; 1000],
            (0..=255u8).cycle().take(5000).collect(),
        ]
    }

    #[test]
    fn order0_roundtrip_various() {
        for data in sample_inputs() {
            let enc = encode_o0(&data);
            assert_eq!(decode(&enc).unwrap(), data, "order 0, len {}", data.len());
        }
    }

    #[test]
    fn order1_roundtrip_various() {
        for data in sample_inputs() {
            let enc = encode_o1(&data);
            assert_eq!(decode(&enc).unwrap(), data, "order 1, len {}", data.len());
        }
    }

    #[test]
    fn order0_hand_derived_bytes() {
        // "AB": two symbols at 2048 each. Table: 65 (0x88 0x00), run to 66,
        // (0x88 0x00), terminator. Data is coded from the end: 'B' takes state 1
        // to 2^24 + 2048, 'A' takes state 0 to 2^24. States 2 and 3 stay at LOWER.
        let enc = encode_o0(b"AB");
        let expected: Vec<u8> = [
            &[0u8][..],
            &[24, 0, 0, 0],
            &[2, 0, 0, 0],
            &[65, 0x88, 0x00, 66, 0x00, 0x88, 0x00, 0x00],
            &[0x00, 0x00, 0x00, 0x01],
            &[0x00, 0x08, 0x00, 0x01],
            &[0x00, 0x00, 0x80, 0x00],
            &[0x00, 0x00, 0x80, 0x00],
        ]
        .concat();
        assert_eq!(enc, expected);
        assert_eq!(decode(&enc).unwrap(), b"AB");
    }

    #[test]
    fn order1_hand_derived_bytes() {
        // "AAAA": one symbol per stream, each with frequency 4096, which leaves
        // every state at LOWER. The context list is [0]; the context-0 table
        // is A = 4096 (0x90 0x00) followed by a terminator.
        let enc = encode_o1(b"AAAA");
        let mut expected = vec![1u8, 22, 0, 0, 0, 4, 0, 0, 0];
        expected.extend_from_slice(&[0, 65, 0x90, 0x00, 0, 0]);
        for _ in 0..4 {
            expected.extend_from_slice(&[0x00, 0x00, 0x80, 0x00]);
        }
        assert_eq!(enc, expected);
        assert_eq!(decode(&enc).unwrap(), b"AAAA");
    }

    #[test]
    fn bad_size_is_malformed() {
        let mut enc = encode_o0(b"ACGT");
        enc[1..5].copy_from_slice(&1_000_000u32.to_le_bytes());
        assert!(matches!(decode(&enc), Err(SeqError::Malformed { .. })));
    }

    #[test]
    fn unknown_order_is_malformed() {
        let mut enc = encode_o1(b"ACGT");
        enc[0] = 2;
        assert!(matches!(decode(&enc), Err(SeqError::Malformed { .. })));
    }

    #[test]
    fn bad_table_sum_is_malformed() {
        // Order 0 table for "AB": raise 65's frequency from 2048 to 2049, sum 4097.
        let mut enc = encode_o0(b"AB");
        enc[9 + 2] = 0x01;
        assert!(matches!(decode(&enc), Err(SeqError::Malformed { .. })));
    }
}
