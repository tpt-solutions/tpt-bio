"""Smoke tests for the tpt_bio extension. Build the module first, for example
with `maturin develop` or by copying the cdylib to a file named tpt_bio.pyd / .so.
"""

import pytest

tpt_bio = pytest.importorskip("tpt_bio")


def test_reverse_complement_preserves_case():
    assert tpt_bio.reverse_complement(b"ACGTN") == b"NACGT"
    assert tpt_bio.reverse_complement(b"acgRYm") == b"kRYcgt"


def test_reverse_complement_rejects_invalid():
    with pytest.raises(ValueError):
        tpt_bio.reverse_complement(b"AC?")


def test_gc_content():
    assert tpt_bio.gc_content(b"GGCCAT") == pytest.approx(4 / 6)
    assert tpt_bio.gc_content(b"NNNN") is None


def test_canonical_kmers():
    assert tpt_bio.canonical_kmers(b"AAAC", 3) == [(0, 0), (1, 1)]


def test_canonical_kmers_rejects_bad_k():
    with pytest.raises(ValueError):
        tpt_bio.canonical_kmers(b"ACGT", 0)
