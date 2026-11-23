"""Type stubs for the tpt_bio extension module."""

from collections.abc import Iterator
from typing import Literal, overload

import numpy as np
from numpy.typing import NDArray

__version__: str

# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------------

class FastqRecord:
    """One FASTQ record. `name` excludes the leading '@'; `qual` is Phred+33 text."""

    @property
    def name(self) -> bytes: ...
    @property
    def seq(self) -> bytes: ...
    @property
    def qual(self) -> bytes: ...

class FastqReader:
    """Streaming FASTQ reader. Gzip-compressed files are detected automatically.

    Iterating yields FastqRecord objects. Raises OSError if the file cannot be
    opened and ValueError on malformed input.
    """

    def __init__(self, path: str) -> None: ...
    def __iter__(self) -> FastqReader: ...
    def __next__(self) -> FastqRecord: ...
    def read_batch(self, n: int) -> list[FastqRecord]:
        """Reads up to n records. Returns fewer, or an empty list, at end of file."""
        ...

# ---------------------------------------------------------------------------
# Indexed BAM
# ---------------------------------------------------------------------------

class BamRecord:
    """One alignment record. `pos` is 0-based (-1 when unset). `seq` is ASCII
    bases. `qual` holds raw Phred scores (not +33); empty means unavailable.
    """

    @property
    def name(self) -> bytes: ...
    @property
    def flags(self) -> int: ...
    @property
    def ref_id(self) -> int: ...
    @property
    def pos(self) -> int: ...
    @property
    def mapq(self) -> int: ...
    @property
    def cigar(self) -> str: ...
    @property
    def next_ref_id(self) -> int: ...
    @property
    def next_pos(self) -> int: ...
    @property
    def tlen(self) -> int: ...
    @property
    def seq(self) -> bytes: ...
    @property
    def qual(self) -> bytes: ...

class IndexedBam:
    """A coordinate-sorted BAM file with its BAI index."""

    def __init__(self, bam_path: str, bai_path: str) -> None: ...
    @property
    def references(self) -> list[tuple[str, int]]:
        """Reference names and lengths, in header order."""
        ...
    def fetch(self, region: str) -> list[BamRecord]:
        """Records overlapping region: 'contig', 'contig:start' or
        'contig:start-end' (1-based, inclusive).

        Raises ValueError for an unknown reference or a malformed region.
        """
        ...

def build_bai(bam_path: str, bai_path: str) -> None:
    """Builds a BAI index for a coordinate-sorted BAM and writes it to bai_path."""
    ...

# ---------------------------------------------------------------------------
# Vision (2D images as NumPy arrays, shape (rows, cols))
# ---------------------------------------------------------------------------

Numeric2D = NDArray[np.uint8] | NDArray[np.uint16] | NDArray[np.uint32] | NDArray[np.int8] | NDArray[np.int16] | NDArray[np.int32] | NDArray[np.float32] | NDArray[np.float64]

def otsu_threshold(image: Numeric2D, bins: int = 256) -> float | None:
    """Otsu threshold of a 2D image, or None if the image is empty or constant."""
    ...

def apply_threshold(image: Numeric2D, threshold_value: float) -> NDArray[np.bool_]:
    """Boolean foreground mask: True where image > threshold_value."""
    ...

def label_components(
    mask: NDArray[np.bool_], connectivity: Literal[4, 8] = 8
) -> tuple[NDArray[np.uint32], int]:
    """Labels connected foreground components of a 2D boolean mask.

    Returns (labels, count). Background is 0; components are 1..=count.
    connectivity is 4 (edge-sharing) or 8 (edge- and corner-sharing).
    """
    ...

def watershed_2d(
    image: Numeric2D,
    markers: NDArray[np.uint32],
    mask: NDArray[np.bool_] | None = None,
) -> NDArray[np.uint32]:
    """Marker-controlled watershed. markers holds seed labels (0 = unlabelled).
    mask, if given, restricts flooding to True pixels.
    """
    ...

def load_image(path: str, channel: int = 0) -> tuple[NDArray[np.generic], int]:
    """Loads one channel of a PNG or TIFF file as a 2D array, keeping its sample type
    (uint8, uint16 or float32). Returns (array, channel_count).

    Raises OSError if the file cannot be opened, ValueError if it cannot be decoded,
    is a volume, or channel is out of range.
    """
    ...

def object_features(
    labels: NDArray[np.uint32],
    intensity: Numeric2D | None = None,
) -> dict[str, list[float]]:
    """Per-object features of a 2D label image (background 0, objects 1..=n).

    Returns a dict mapping column name to one value per object, in label order.
    Shape columns are label, area, perimeter, eccentricity and solidity. When
    intensity is given, mean, std, min, max, q25, median, q75 and integrated are
    added. Intensity must have the same shape as labels.
    """
    ...

# ---------------------------------------------------------------------------
# Pipelines
# ---------------------------------------------------------------------------

PipelineOutput = (
    NDArray[np.float64]
    | NDArray[np.uint64]
    | list[FastqRecord]
    | list[BamRecord]
    | NDArray[np.bool_]
    | tuple[NDArray[np.uint32], int]
    | dict[str, list[float]]
    | str
)

class Pipeline:
    """A pipeline loaded from a TOML definition. See the pipeline crate docs for the
    node kinds. Outputs are returned in definition order.
    """

    @staticmethod
    def from_toml(path: str) -> Pipeline:
        """Reads and validates a TOML pipeline file.

        Raises OSError if the file cannot be read and ValueError if it is invalid.
        """
        ...

    def run(self) -> list[PipelineOutput]:
        """Runs every node and returns the outputs as Python objects:

        - table: float64 array, shape (rows, cols)
        - labels: uint64 array, one label per row
        - fastq: list of FastqRecord
        - records: list of BamRecord
        - file: str path written by a sink
        - image and mask: 2D arrays (float64 and bool)
        - objects: (uint32 labels, count) tuple
        - frame: dict of column name to list of floats

        Raises RuntimeError if a node fails.
        """
        ...
