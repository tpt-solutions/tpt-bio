"""Type stubs for the tpt_bio extension module."""

__version__: str

def reverse_complement(seq: bytes) -> bytes:
    """Reverse complement of an IUPAC DNA sequence. Case is preserved.

    Raises ValueError on the first byte that is not an IUPAC code.
    """
    ...

def gc_content(seq: bytes) -> float | None:
    """Fraction of unambiguous bases that are G or C, or None if there are none."""
    ...

def canonical_kmers(seq: bytes, k: int) -> list[tuple[int, int]]:
    """Canonical k-mers as (offset, packed 2-bit value) pairs. Windows with
    non-ACGT bytes are skipped. Raises ValueError if k is not in 1..=32.
    """
    ...
