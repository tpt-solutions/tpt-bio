//! TOK3 read-name tokeniser, decoding only (CRAM block method 8).
//!
//! Ported from htscodecs `tokenise_name3.c` (`tok3_decode_names` and `decode_name`),
//! which samtools and htslib use for the CRAM 3.1 read-name block. Names are split
//! into tokens (alphanumeric runs, digit runs, single characters). Token position `n`
//! has a type stream and one value stream per type. Each name is rebuilt from its
//! earlier name `pnum`: token 0 is either a duplicate of `pnum` or a name that
//! refers back to it by distance, and each later token position takes a type
//! (literal, match against `pnum`, or a numeric delta).
//!
//! Blob layout: raw length (u32 LE), name count (u32 LE), an arithmetic-coder flag
//! (u8), then descriptors. A descriptor byte has bit 7 for a new token position, bit
//! 6 for a copy of an earlier stream (two bytes name the source), and the low nibble
//! for the stream type. A coded stream is a uint7 length followed by the stream's
//! body, which is a rANS Nx16 stream when the flag is 0 and an adaptive arithmetic
//! stream (see [`super::arith`]) when it is 1. uint7 is big-endian base-128.

use crate::cram::{arith, rans_nx16};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM TOK3 read names";
const HEADER_LEN: usize = 9;
const MAX_TOKENS: usize = 128;
const MAX_STREAMS: usize = MAX_TOKENS << 4;
/// Largest name count accepted, as in the reference (a guard against corrupt sizes).
const MAX_NAMES: usize = 10_000_000;

// Token types, numbered as in the reference `enum name_type`. Types 0 (type
// stream of token 0), 5 (`N_DUP`) and 6 (`N_DIFF`) are not token kinds here.
const N_ALPHA: u8 = 1;
const N_CHAR: u8 = 2;
const N_DIGITS0: u8 = 3;
const N_DZLEN: u8 = 4;
const N_DUP: u8 = 5;
const N_DIGITS: u8 = 7;
const N_DDELTA: u8 = 8;
const N_DDELTA0: u8 = 9;
const N_MATCH: u8 = 10;
const N_NOP: u8 = 11;
const N_END: u8 = 12;

/// Decodes a TOK3 blob into the NUL-terminated read names it holds.
///
/// `raw_size` is the block's declared raw size, which must equal the blob's own
/// declared length.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for the bzip2 and order-2 forms of the
/// arithmetic streams, and [`SeqError::Malformed`] for truncated or inconsistent data.
pub fn decode(blob: &[u8], raw_size: usize) -> Result<Vec<u8>> {
    if blob.len() < HEADER_LEN {
        return Err(malformed("header truncated"));
    }
    let ulen = u32_le(&blob[0..4]);
    let nreads = u32_le(&blob[4..8]);
    let arithmetic = match blob[8] {
        0 => false,
        1 => true,
        flag => return Err(malformed(format!("arithmetic flag {flag} is not 0 or 1"))),
    };
    if ulen != raw_size {
        return Err(malformed(format!("blob declares {ulen} bytes, block declares {raw_size}")));
    }
    if nreads == 0 || nreads > MAX_NAMES {
        return Err(malformed(format!("name count {nreads} is out of range")));
    }

    let (streams, max_tok) = read_streams(blob, nreads, arithmetic)?;
    let mut dec = Decoder { streams, max_tok, nreads };
    let mut names: Vec<Name> = Vec::with_capacity(nreads.min(1 << 16));
    let mut out = Vec::with_capacity(ulen.min(1 << 24));
    while let Some(name) = dec.next_name(&names)? {
        out.extend_from_slice(&name.bytes);
        names.push(name);
        if out.len() > ulen {
            return Err(malformed("names exceed the declared length"));
        }
    }
    if out.len() != ulen {
        return Err(malformed(format!("decoded {} bytes, blob declares {ulen}", out.len())));
    }
    Ok(out)
}

fn malformed(reason: impl Into<String>) -> SeqError {
    SeqError::malformed(FORMAT, reason)
}

fn u32_le(b: &[u8]) -> usize {
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize
}

/// Reads a big-endian base-128 unsigned integer of at most five bytes.
/// Returns the value and the number of bytes it occupied.
pub(super) fn read_uint7(b: &[u8]) -> Result<(usize, usize)> {
    let mut value: u64 = 0;
    for (i, &byte) in b.iter().take(5).enumerate() {
        value = (value << 7) | u64::from(byte & 0x7f);
        if byte & 0x80 == 0 {
            let value = usize::try_from(value).map_err(|_| malformed("length overflows"))?;
            return Ok((value, i + 1));
        }
    }
    Err(malformed("length varint is truncated or too long"))
}

/// Unpacks the descriptors into the token streams, indexed by `(token << 4) | type`.
/// Returns the streams and the number of token positions defined.
fn read_streams(
    blob: &[u8],
    nreads: usize,
    arithmetic: bool,
) -> Result<(Vec<Option<Stream>>, usize)> {
    let mut streams: Vec<Option<Stream>> = (0..MAX_STREAMS).map(|_| None).collect();
    let mut max_tok = 1usize;
    let mut tnum: Option<usize> = None;
    let mut o = HEADER_LEN;
    while o < blob.len() {
        let ttype = blob[o];
        o += 1;
        let low = usize::from(ttype & 15);
        let copy = ttype & 64 != 0;
        let mut source = 0usize;
        if copy {
            let pair = blob.get(o..o + 2).ok_or_else(|| malformed("descriptor truncated"))?;
            source = (usize::from(pair[0]) << 4) + usize::from(pair[1]);
            o += 2;
        }
        if ttype & 128 != 0 {
            let next = tnum.map_or(0, |t| t + 1);
            if next >= MAX_TOKENS {
                return Err(malformed("too many token positions"));
            }
            tnum = Some(next);
            max_tok = next + 1;
            for s in &mut streams[next << 4..(next << 4) + 16] {
                *s = None;
            }
            // A new position whose stream type is not 0 starts with a constant type
            // stream: every name is a match except the first, which has this type.
            if low != 0 {
                let mut data = vec![N_MATCH; nreads];
                data[0] = low as u8;
                streams[next << 4] = Some(Stream { data, pos: 0 });
            }
        }
        let t = tnum.ok_or_else(|| malformed("descriptor before the first token position"))?;
        let id = (t << 4) | low;
        if copy {
            if source >= id {
                return Err(malformed("copy descriptor names a later stream"));
            }
            let data = streams[source]
                .as_ref()
                .ok_or_else(|| malformed("copy descriptor names a missing stream"))?
                .data
                .clone();
            streams[id] = Some(Stream { data, pos: 0 });
            continue;
        }
        let rest = &blob[o..];
        let (clen, nb) = read_uint7(rest)?;
        let (ulen, _) = read_uint7(rest.get(nb + 1..).unwrap_or(&[]))?;
        let body = rest
            .get(nb..nb.saturating_add(clen))
            .ok_or_else(|| malformed("coded stream truncated"))?;
        let data =
            if arithmetic { arith::decode(body, ulen)? } else { rans_nx16::decode(body, ulen)? };
        if data.len() != ulen {
            return Err(malformed("coded stream length mismatch"));
        }
        streams[id] = Some(Stream { data, pos: 0 });
        o += nb + clen;
    }
    Ok((streams, max_tok))
}

/// One token stream: its bytes and the read position.
struct Stream {
    data: Vec<u8>,
    pos: usize,
}

impl Stream {
    fn next_byte(&mut self) -> Option<u8> {
        let b = *self.data.get(self.pos)?;
        self.pos += 1;
        Some(b)
    }

    fn next_u32(&mut self) -> Option<u32> {
        let b = self.data.get(self.pos..self.pos + 4)?;
        let v = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        self.pos += 4;
        Some(v)
    }
}

/// A decoded token. `int` is the value (a digit run's number, a character code, a
/// delta base) and `aux` is the field width (digit runs) or the offset in the name
/// (alphanumeric runs), as in the reference `last_context_tok`.
#[derive(Clone, Copy, Default)]
struct Tok {
    kind: u8,
    int: u32,
    aux: usize,
}

/// A decoded name with its NUL terminator, and its tokens indexed by position.
#[derive(Clone, Default)]
struct Name {
    bytes: Vec<u8>,
    toks: Vec<Tok>,
}

impl Name {
    /// The last token position of the name; zero for the empty name.
    fn last_ntok(&self) -> usize {
        self.toks.len().saturating_sub(1)
    }
}

struct Decoder {
    streams: Vec<Option<Stream>>,
    max_tok: usize,
    nreads: usize,
}

/// Stream index for token position `ntok` and type `kind`.
fn sid(ntok: usize, kind: u8) -> usize {
    (ntok << 4) | usize::from(kind)
}

impl Decoder {
    fn stream(&mut self, id: usize) -> Result<&mut Stream> {
        self.streams[id].as_mut().ok_or_else(|| malformed("token stream missing"))
    }

    /// The next token type at position `ntok`, or `None` when that stream is exhausted.
    fn type_at(&mut self, ntok: usize) -> Option<u8> {
        self.streams[ntok << 4].as_mut()?.next_byte()
    }

    fn u32_at(&mut self, id: usize) -> Result<u32> {
        self.stream(id)?.next_u32().ok_or_else(|| malformed("integer stream exhausted"))
    }

    fn byte_at(&mut self, id: usize) -> Result<u8> {
        self.stream(id)?.next_byte().ok_or_else(|| malformed("byte stream exhausted"))
    }

    /// A NUL-terminated string, without its terminator.
    fn alpha_at(&mut self, id: usize) -> Result<Vec<u8>> {
        let s = self.stream(id)?;
        let mut out = Vec::new();
        loop {
            let c = s.next_byte().ok_or_else(|| malformed("string stream exhausted"))?;
            if c == 0 {
                return Ok(out);
            }
            out.push(c);
        }
    }

    /// Decodes the name after the names in `names`, or `None` at the end of the names.
    fn next_name(&mut self, names: &[Name]) -> Result<Option<Name>> {
        let cnum = names.len();
        if cnum > self.nreads {
            return Err(malformed("more names than the blob declares"));
        }
        // Token 0 is the duplicate-or-distance type; an exhausted stream ends the names.
        let t0 = match self.type_at(0) {
            Some(t) => usize::from(t),
            None => return Ok(None),
        };
        if t0 >= self.max_tok * 16 {
            return Ok(None);
        }
        let dist = self.u32_at(t0)?;
        if u64::from(dist) > cnum as u64 {
            return Err(malformed("name distance points before the first name"));
        }
        let pnum = cnum - dist as usize;
        let empty = Name::default();
        let prev = if pnum < cnum { &names[pnum] } else { &empty };

        if t0 == usize::from(N_DUP) {
            if pnum == cnum {
                return Err(malformed("duplicate name refers to itself"));
            }
            return Ok(Some(prev.clone()));
        }

        let mut bytes = Vec::new();
        let mut toks = vec![Tok::default()];
        for ntok in 1..self.max_tok.min(MAX_TOKENS) {
            // An exhausted or unknown type ends the name, as the reference's elided END.
            let kind = self.type_at(ntok).unwrap_or(N_END);
            let tok = match kind {
                N_CHAR => {
                    let c = self.byte_at(sid(ntok, N_CHAR))?;
                    bytes.push(c);
                    Tok { kind: N_CHAR, int: u32::from(c), aux: 0 }
                }
                N_ALPHA => {
                    let s = self.alpha_at(sid(ntok, N_ALPHA))?;
                    let aux = bytes.len();
                    let int = s.len() as u32;
                    bytes.extend_from_slice(&s);
                    Tok { kind: N_ALPHA, int, aux }
                }
                N_DIGITS0 => {
                    let width = self.byte_at(sid(ntok, N_DZLEN))?;
                    let v = self.u32_at(sid(ntok, N_DIGITS0))?;
                    push_fixed(&mut bytes, v, usize::from(width));
                    Tok { kind: N_DIGITS0, int: v, aux: usize::from(width) }
                }
                N_DDELTA0 => {
                    check_match_position(ntok, prev)?;
                    let delta = self.byte_at(sid(ntok, N_DDELTA0))?;
                    let base = prev.toks[ntok];
                    let v = u32::from(delta).wrapping_add(base.int);
                    push_fixed(&mut bytes, v, base.aux);
                    Tok { kind: N_DIGITS0, int: v, aux: base.aux }
                }
                N_DIGITS => {
                    let v = self.u32_at(sid(ntok, N_DIGITS))?;
                    push_var(&mut bytes, v);
                    Tok { kind: N_DIGITS, int: v, aux: 0 }
                }
                N_DDELTA => {
                    check_match_position(ntok, prev)?;
                    let delta = self.byte_at(sid(ntok, N_DDELTA))?;
                    let v = u32::from(delta).wrapping_add(prev.toks[ntok].int);
                    push_var(&mut bytes, v);
                    Tok { kind: N_DIGITS, int: v, aux: 0 }
                }
                N_NOP => Tok { kind: N_NOP, int: 0, aux: 0 },
                N_MATCH => copy_match(&mut bytes, prev, ntok)?,
                _ => {
                    bytes.push(0);
                    toks.push(Tok { kind: N_END, int: 0, aux: 0 });
                    return Ok(Some(Name { bytes, toks }));
                }
            };
            toks.push(tok);
        }
        Err(malformed("name has no end token"))
    }
}

/// A match token: copies the token at `ntok` of the previous name `prev`.
fn copy_match(bytes: &mut Vec<u8>, prev: &Name, ntok: usize) -> Result<Tok> {
    check_match_position(ntok, prev)?;
    let base = prev.toks[ntok];
    Ok(match base.kind {
        N_CHAR => {
            bytes.push(base.int as u8);
            Tok { kind: N_CHAR, int: base.int, aux: 0 }
        }
        N_ALPHA => {
            let len = base.int as usize;
            let src = prev
                .bytes
                .get(base.aux..base.aux.saturating_add(len))
                .ok_or_else(|| malformed("match runs past the previous name"))?;
            let aux = bytes.len();
            bytes.extend_from_slice(src);
            Tok { kind: N_ALPHA, int: base.int, aux }
        }
        N_DIGITS => {
            push_var(bytes, base.int);
            Tok { kind: N_DIGITS, int: base.int, aux: 0 }
        }
        N_DIGITS0 => {
            push_fixed(bytes, base.int, base.aux);
            Tok { kind: N_DIGITS0, int: base.int, aux: base.aux }
        }
        _ => return Err(malformed("match against an unmatchable token")),
    })
}

fn check_match_position(ntok: usize, prev: &Name) -> Result<()> {
    if ntok >= prev.last_ntok() {
        return Err(malformed("token position is past the previous name"));
    }
    Ok(())
}

/// Appends `v` zero-padded to `width` digits.
fn push_fixed(out: &mut Vec<u8>, v: u32, width: usize) {
    out.extend_from_slice(format!("{v:0width$}").as_bytes());
}

/// Appends `v` in decimal. The reference prints nothing for zero, and so does this.
fn push_var(out: &mut Vec<u8>, v: u32) {
    if v != 0 {
        out.extend_from_slice(v.to_string().as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// uint7 encoding of `v`, as the reference writes it.
    fn put_uint7(out: &mut Vec<u8>, v: usize) {
        let mut groups = vec![(v & 0x7f) as u8];
        let mut rest = v >> 7;
        while rest > 0 {
            groups.push(((rest & 0x7f) as u8) | 0x80);
            rest >>= 7;
        }
        groups.reverse();
        out.extend_from_slice(&groups);
    }

    /// A coded stream: uint7 length, then the rANS Nx16 encoding of `data`.
    fn coded(data: &[u8]) -> Vec<u8> {
        let enc = rans_nx16::encode(data);
        let mut out = Vec::new();
        put_uint7(&mut out, enc.len());
        out.extend_from_slice(&enc);
        out
    }

    /// Assembles a blob from (descriptor, body) pairs, with the header for `names`.
    fn blob(ulen: usize, nreads: u32, parts: &[(u8, Vec<u8>)]) -> Vec<u8> {
        let mut out = (ulen as u32).to_le_bytes().to_vec();
        out.extend_from_slice(&nreads.to_le_bytes());
        out.push(0);
        for (desc, body) in parts {
            out.push(*desc);
            out.extend_from_slice(body);
        }
        out
    }

    /// Descriptors for names "a01" and "a02": token 0 is DIFF (6) and the second name
    /// refers back by distance 1. Token 1 is ALPHA "a" then a match. Token 2 is
    /// DIGITS0 width 2 (value 1) then a DDELTA0 of +1, which gives "02".
    fn zero_padded_parts() -> Vec<(u8, Vec<u8>)> {
        zero_padded_parts_with(coded)
    }

    /// The parts of [`zero_padded_parts`], with each stream encoded by `code`.
    fn zero_padded_parts_with(code: fn(&[u8]) -> Vec<u8>) -> Vec<(u8, Vec<u8>)> {
        vec![
            (0x80, code(&[6, 6])),                   // token 0 types
            (0x06, code(&[0, 0, 0, 0, 1, 0, 0, 0])), // DIFF distances (LE)
            (0x80, code(&[1, 10])),                  // token 1 types: ALPHA, MATCH
            (0x01, code(b"a\0")),                    // ALPHA "a"
            (0x80, code(&[3, 9])),                   // token 2 types: DIGITS0, DDELTA0
            (0x03, code(&[1, 0, 0, 0])),             // DIGITS0 value 1
            (0x04, code(&[2])),                      // DIGITS0 width 2
            (0x09, code(&[1])),                      // DDELTA0 +1
            (0x80, code(&[12, 12])),                 // token 3 types: END, END
        ]
    }

    /// An arithmetic stream stored raw (flag bit 5, `X_CAT`), preceded by its uint7
    /// length. The reference writes this form for tiny streams.
    fn stored_arith(data: &[u8]) -> Vec<u8> {
        let mut body = vec![0x20];
        put_uint7(&mut body, data.len());
        body.extend_from_slice(data);
        let mut out = Vec::new();
        put_uint7(&mut out, body.len());
        out.extend_from_slice(&body);
        out
    }

    /// Blob with the arithmetic flag set.
    fn arith_blob(ulen: usize, nreads: u32, parts: &[(u8, Vec<u8>)]) -> Vec<u8> {
        let mut out = blob(ulen, nreads, parts);
        out[8] = 1;
        out
    }

    #[test]
    fn decodes_zero_padded_numeric_tokens() {
        let want = b"a01\0a02\0";
        let parts = zero_padded_parts();
        let data = blob(want.len(), 2, &parts);
        assert_eq!(decode(&data, want.len()).unwrap(), want);
    }

    #[test]
    fn decodes_duplicate_name() {
        // Name 0 is "x". Name 1 is a duplicate of name 0 (type DUP, distance 1).
        let want = b"x\0x\0";
        let parts = vec![
            (0x80, coded(&[6, 5])),       // token 0 types: DIFF, DUP
            (0x06, coded(&[0, 0, 0, 0])), // DIFF distance 0 for name 0
            (0x05, coded(&[1, 0, 0, 0])), // DUP distance 1 for name 1
            (0x80, coded(&[1])),          // token 1: ALPHA for name 0 only
            (0x01, coded(b"x\0")),
            (0x80, coded(&[12])), // token 2: END for name 0 only
        ];
        let data = blob(want.len(), 2, &parts);
        assert_eq!(decode(&data, want.len()).unwrap(), want);
    }

    #[test]
    fn decodes_constant_type_stream() {
        // Names "r1" and "r2" with a constant type stream at token 1 (ALPHA "r" for
        // name 0, MATCH for name 1), and DIGITS 1 then DDELTA +1 at token 2.
        let want = b"r1\0r2\0";
        let parts = vec![
            (0x80, coded(&[6, 6])),
            (0x06, coded(&[0, 0, 0, 0, 1, 0, 0, 0])),
            (0x81, coded(b"r\0")),        // constant type ALPHA: only name 0
            (0x80, coded(&[7, 8])),       // token 2: DIGITS, DDELTA
            (0x07, coded(&[1, 0, 0, 0])), // DIGITS 1
            (0x08, coded(&[1])),          // DDELTA +1
            (0x80, coded(&[12, 12])),
        ];
        let data = blob(want.len(), 2, &parts);
        assert_eq!(decode(&data, want.len()).unwrap(), want);
    }

    #[test]
    fn rejects_length_mismatch_and_truncation() {
        let want = b"a01\0a02\0";
        let parts = zero_padded_parts();
        let data = blob(want.len(), 2, &parts);
        assert!(decode(&data, want.len() + 1).is_err());
        assert!(decode(&data[..data.len() - 3], want.len()).is_err());
        assert!(decode(&data[..5], want.len()).is_err());
    }

    #[test]
    fn decodes_stored_arithmetic_streams() {
        let want = b"a01\0a02\0";
        let parts = zero_padded_parts_with(stored_arith);
        let data = arith_blob(want.len(), 2, &parts);
        assert_eq!(decode(&data, want.len()).unwrap(), want);
    }

    #[test]
    fn rejects_unknown_arithmetic_flag() {
        let mut data = blob(0, 1, &[]);
        data[8] = 2;
        assert!(matches!(decode(&data, 0), Err(SeqError::Malformed { .. })));
    }
}
