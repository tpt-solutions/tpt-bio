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
