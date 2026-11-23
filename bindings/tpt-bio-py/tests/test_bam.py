"""Indexed BAM tests. The fixture BAM is written by the small encoder below, so
the tests need no binary fixtures in the repository."""

import struct
import zlib

import pytest

tpt_bio = pytest.importorskip("tpt_bio")

SEQ_CODES = "=ACMGRSVTWYHKDBN"


def _bgzf_block(payload: bytes) -> bytes:
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = co.compress(payload) + co.flush()
    total = 18 + len(cdata) + 8
    header = (
        b"\x1f\x8b\x08\x04"
        + b"\x00\x00\x00\x00"
        + b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<HH", 2, total - 1)
    )
    trailer = struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))
    return header + cdata + trailer


def _reg2bin(beg: int, end: int) -> int:
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def _encode_record(name: bytes, ref_id: int, pos: int, seq: bytes, qual: bytes, cigar_len: int) -> bytes:
    name_z = name + b"\x00"
    packed = bytearray((len(seq) + 1) // 2)
    for i, base in enumerate(seq):
        code = SEQ_CODES.index(chr(base).upper())
        if i % 2 == 0:
            packed[i // 2] |= code << 4
        else:
            packed[i // 2] |= code
    cigar = struct.pack("<I", (cigar_len << 4) | 0)  # cigar_len M
    body = struct.pack(
        "<iiBBHHHiiii",
        ref_id,
        pos,
        len(name_z),
        60,
        _reg2bin(pos, pos + cigar_len),
        1,
        0,
        len(seq),
        -1,
        -1,
        0,
    )
    body += name_z + cigar + bytes(packed) + bytes(qual)
    return struct.pack("<i", len(body)) + body


def write_bam(path, refs, records):
    """refs: list of (name, length). records: list of (name, ref_id, pos, seq, qual)."""
    text = b"@HD\tVN:1.6\tSO:coordinate\n"
    payload = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(refs))
    for name, length in refs:
        nz = name.encode() + b"\x00"
        payload += struct.pack("<i", len(nz)) + nz + struct.pack("<i", length)
    for name, ref_id, pos, seq, qual in records:
        payload += _encode_record(name, ref_id, pos, seq, qual, len(seq))
    # Split into blocks that stay well under the 64 KiB BGZF limit.
    data = b""
    for i in range(0, len(payload), 0xFF00):
        data += _bgzf_block(payload[i : i + 0xFF00])
    data += _bgzf_block(b"")  # standard EOF marker
    path.write_bytes(data)


@pytest.fixture
def bam_pair(tmp_path):
    bam = tmp_path / "t.bam"
    bai = tmp_path / "t.bam.bai"
    # Reads of 10 bases each. Query chr1:100-200 (1-based) overlaps only "hit".
    records = [
        (b"early", 0, 0, b"ACGTACGTAC", b"\x1e" * 10),
        (b"hit", 0, 120, b"GGGGCCCCAA", b"\x28" * 10),
        (b"late", 0, 300, b"TTTTTTTTTT", b"\x28" * 10),
        (b"chr2read", 1, 5, b"ACGTACGTAC", b"\x28" * 10),
    ]
    write_bam(bam, [("chr1", 2000), ("chr2", 500)], records)
    tpt_bio.build_bai(str(bam), str(bai))
    return str(bam), str(bai)


def test_build_bai_writes_index(bam_pair):
    _, bai = bam_pair
    with open(bai, "rb") as f:
        assert f.read(4) == b"BAI\x01"


def test_references(bam_pair):
    bam, bai = bam_pair
    reader = tpt_bio.IndexedBam(bam, bai)
    assert reader.references == [("chr1", 2000), ("chr2", 500)]


def test_fetch_region_returns_overlapping_records(bam_pair):
    bam, bai = bam_pair
    reader = tpt_bio.IndexedBam(bam, bai)
    hits = reader.fetch("chr1:100-200")
    assert [r.name for r in hits] == [b"hit"]
    rec = hits[0]
    assert rec.pos == 120
    assert rec.ref_id == 0
    assert rec.cigar == "10M"
    assert rec.seq == b"GGGGCCCCAA"
    assert rec.qual == bytes([0x28] * 10)


def test_fetch_whole_contig(bam_pair):
    bam, bai = bam_pair
    reader = tpt_bio.IndexedBam(bam, bai)
    names = [r.name for r in reader.fetch("chr2")]
    assert names == [b"chr2read"]


def test_fetch_unknown_reference_raises(bam_pair):
    bam, bai = bam_pair
    reader = tpt_bio.IndexedBam(bam, bai)
    with pytest.raises(ValueError):
        reader.fetch("chrX:1-10")


def test_fetch_bad_region_raises(bam_pair):
    bam, bai = bam_pair
    reader = tpt_bio.IndexedBam(bam, bai)
    with pytest.raises(ValueError):
        reader.fetch("chr1:200-100")
