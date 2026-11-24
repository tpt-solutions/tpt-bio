# Python quickstart

The `tpt_bio` package exposes the Rust crates to Python through PyO3. It covers
FASTQ reading, indexed BAM queries, 2D image thresholding, labelling and per-object
features, and running TOML pipelines. Heavy image and sequence work happens in Rust, and
arrays are passed as NumPy arrays.

The package is not on PyPI yet. This tutorial builds a wheel from the repository and
installs it into a virtual environment. Nothing is uploaded.

## 1. Build a wheel

You need Python 3.9 or later, Rust 1.95 (from `rust-toolchain.toml`) and maturin:

```sh
pip install maturin
cd bindings/tpt-bio-py
maturin build --release --out ../../target/wheels
cd ../..
```

The wheel is an `abi3` build, so one wheel works for every Python from 3.9 up on
its platform. On Windows the file is named like `tpt_bio-0.1.0-cp39-abi3-win_amd64.whl`.

## 2. Install it into a virtual environment

```sh
python -m venv .venv
source .venv/Scripts/activate     # Windows (Git Bash); use .venv/bin/activate on Linux and macOS
pip install --force-reinstall numpy ./target/wheels/tpt_bio-*.whl
python -c "import tpt_bio; print(tpt_bio.reverse_complement(b'ACGT'))"
```

`--force-reinstall` matters when you rebuild the wheel. Without it, pip can keep a
previously installed copy with the same version number.

## 3. Generate an input image

The Python example reads a PNG. Create one with the vision example from the vision
tutorial, run from the repository root:

```sh
cargo run -q -p tpt-bio-vision --example vision_quickstart -- vision-demo
```

## 4. Run the quickstart script

The script below is `bindings/tpt-bio-py/examples/quickstart.py`. It loads the image,
thresholds and labels it, prints per-object features, and runs a TOML pipeline. It is
reproduced here exactly as it is in the repository.

```python
"""tpt-bio quickstart: load an image, measure its objects and run a pipeline.

Usage: python quickstart.py <scene.png> [work-dir]

Create the input with the vision example from the repository root:
    cargo run -q -p tpt-bio-vision --example vision_quickstart -- vision-demo
    python quickstart.py vision-demo/scene.png python-demo
"""

import sys
from pathlib import Path

import tpt_bio

png = Path(sys.argv[1])
work = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
work.mkdir(parents=True, exist_ok=True)

# 1. Load channel 0. The array keeps the file's sample type (uint8 for this PNG).
image, channels = tpt_bio.load_image(str(png), 0)
print(f"loaded {image.shape} {image.dtype}, {channels} channel(s)")

# 2. Threshold, label the foreground and measure every object.
threshold = tpt_bio.otsu_threshold(image)
if threshold is None:
    raise SystemExit("image is empty or constant")
mask = tpt_bio.apply_threshold(image, threshold)
labels, count = tpt_bio.label_components(mask, 8)
features = tpt_bio.object_features(labels, image)
print(f"threshold {threshold:.1f}, {count} objects")
print("columns:", sorted(features))
print("areas:", features["area"])
print("mean intensity:", [round(v, 1) for v in features["mean"]])

# 3. Run a TOML pipeline: a small statistics flow over an in-memory table.
stats = work / "stats.toml"
stats.write_text(
    """
[[nodes]]
name = "table"
type = "table_source"
params = { values = [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]] }

[[nodes]]
name = "scaled"
type = "normalize"
params = { method = "zscore" }
inputs = ["table"]

[[nodes]]
name = "reduced"
type = "pca"
params = { components = 1 }
inputs = ["scaled"]

[[nodes]]
name = "groups"
type = "kmeans"
params = { k = 2, seed = 1 }
inputs = ["reduced"]
"""
)
pipeline = tpt_bio.Pipeline.from_toml(str(stats))
outputs = pipeline.run()  # one entry per node, in definition order
for i, output in enumerate(outputs):
    print(f"output {i}: {type(output).__name__} {getattr(output, 'shape', '')}")
print("cluster labels:", outputs[-1].tolist())
```

Run it with the scene image and a work directory for the pipeline file:

```sh
python bindings/tpt-bio-py/examples/quickstart.py vision-demo/scene.png python-demo
```

Expected output (abbreviated):

```text
loaded (64, 96) uint8, 1 channel(s)
threshold 18.1, 4 objects
areas: [113.0, 49.0, 149.0, 29.0]
output 3: ndarray (4,)
cluster labels: [0, 0, 1, 1]
```

## What each call returns

| Call | Returns |
| --- | --- |
| `load_image(path, channel)` | `(array, channel_count)`; the array keeps uint8, uint16 or float32 |
| `otsu_threshold(image)` | a float, or `None` for an empty or constant image |
| `apply_threshold(image, t)` | boolean NumPy mask, `True` where the image is above `t` |
| `label_components(mask, 8)` | `(labels, count)`; background is 0, objects are `1..=count` |
| `object_features(labels, image)` | a dict of column name to list of floats; intensity columns need `image` |
| `Pipeline.from_toml(path).run()` | a list with one output per node, in definition order |

The stub file `bindings/tpt-bio-py/tpt_bio.pyi` has the full signatures and docstrings.

## Sequence data

`FastqReader` iterates over FASTQ records (gzip detected). `IndexedBam(bam, bai)` opens
an indexed BAM, and `fetch("chr1:1000-2000")` returns the overlapping records.
`build_bai(bam, bai)` writes an index. These use the same code as the `tptbio seq`
commands. The BAM must be coordinate-sorted, as in the sequence tutorial.

## Limits

- Only the functions listed in `tpt_bio.pyi` are exposed. The 3D renderer, VCF, BCF
  and CRAM are not available from Python.
- The vision functions take 2D arrays only.
- Only the sequence utilities release the GIL. Image and pipeline calls hold it.
- The tests in `bindings/tpt-bio-py/tests` run with pytest against a locally built module.
