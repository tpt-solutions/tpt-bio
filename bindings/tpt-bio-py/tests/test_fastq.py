"""FASTQ reader tests for the tpt_bio extension."""

import gzip

import pytest

tpt_bio = pytest.importorskip("tpt_bio")

FASTQ = b"@r1 extra\nACGT\n+\nIIII\n@r2\nGGN\n+\n!!#\n"


def test_fastq_iterates_records(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_bytes(FASTQ)
    recs = list(tpt_bio.FastqReader(str(path)))
    assert [r.name for r in recs] == [b"r1 extra", b"r2"]
    assert [r.seq for r in recs] == [b"ACGT", b"GGN"]
    assert [r.qual for r in recs] == [b"IIII", b"!!#"]
    assert all(isinstance(r.seq, bytes) for r in recs)


def test_fastq_reads_gzip(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(FASTQ))
    recs = list(tpt_bio.FastqReader(str(path)))
    assert len(recs) == 2
    assert recs[0].seq == b"ACGT"


def test_fastq_next_returns_none_style_stop(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_bytes(FASTQ)
    reader = tpt_bio.FastqReader(str(path))
    assert next(reader).name == b"r1 extra"
    assert next(reader).name == b"r2"
    with pytest.raises(StopIteration):
        next(reader)


def test_fastq_read_batch(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_bytes(FASTQ)
    reader = tpt_bio.FastqReader(str(path))
    first = reader.read_batch(1)
    assert [r.name for r in first] == [b"r1 extra"]
    rest = reader.read_batch(10)
    assert [r.name for r in rest] == [b"r2"]
    assert reader.read_batch(5) == []


def test_fastq_missing_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        tpt_bio.FastqReader(str(tmp_path / "missing.fastq"))


def test_fastq_truncated_record_raises_value_error(tmp_path):
    path = tmp_path / "bad.fastq"
    path.write_bytes(b"@r1\nACGT\n")
    with pytest.raises(ValueError):
        list(tpt_bio.FastqReader(str(path)))
