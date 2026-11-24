# tpt-bio (Python)

Python bindings for the tpt-bio Rust crates: sequence utilities, FASTQ and
indexed BAM readers, and 2D image segmentation helpers. Heavy operations
release the GIL.

## Install

The package is not on PyPI yet. Build a wheel from this repository with
[maturin](https://www.maturin.rs/) and install it:

```sh
pip install maturin
cd bindings/tpt-bio-py
maturin build --release
pip install target/wheels/tpt_bio-*.whl
```

Requires Python 3.9 or later. Built as an abi3 wheel with PyO3. See
[docs/tutorials/python.md](../../docs/tutorials/python.md) for a worked example.

## Example

```python
import tpt_bio
import numpy as np

tpt_bio.reverse_complement(b"ACGT")          # b"ACGT"

for rec in tpt_bio.FastqReader("reads.fastq.gz"):
    print(rec.name, len(rec.seq))

img = np.array([[0, 9, 9], [0, 9, 0]], dtype=np.uint8)
t = tpt_bio.otsu_threshold(img)
labels, count = tpt_bio.label_components(tpt_bio.apply_threshold(img, t), 8)
```

## Licence

MIT OR Apache-2.0.
