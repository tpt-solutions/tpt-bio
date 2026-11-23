"""Tests for Pipeline.from_toml(path).run(). Inputs are generated under pytest's tmp_path."""

import numpy as np
import pytest

tpt_bio = pytest.importorskip("tpt_bio")

from test_images import write_png  # noqa: E402  (same directory; pytest puts it on sys.path)


def _squares(h=32, w=32):
    img = np.zeros((h, w), dtype=np.uint8)
    img[4:10, 4:10] = 200  # 6x6 interior object
    img[20:26, 18:24] = 200  # second interior object
    img[12:14, 28:30] = 200  # third object, 2x2
    return img


def _write_definition(tmp_path, body: str):
    path = tmp_path / "pipeline.toml"
    path.write_text(body, encoding="utf-8")
    return path


def test_vision_pipeline_returns_python_objects(tmp_path):
    png = tmp_path / "plate.png"
    write_png(png, _squares())
    toml = _write_definition(
        tmp_path,
        f"""
[[nodes]]
name = "raw"
type = "load_image"
params = {{ path = '{png}' }}

[[nodes]]
name = "mask"
type = "otsu"
inputs = ["raw"]

[[nodes]]
name = "cells"
type = "label"
inputs = ["mask"]

[[nodes]]
name = "features"
type = "measure"
inputs = ["cells", "raw"]

[[nodes]]
name = "count"
type = "count_objects"
inputs = ["cells"]
""",
    )
    outputs = tpt_bio.Pipeline.from_toml(str(toml)).run()
    assert len(outputs) == 5

    raw, mask, cells, features, count = outputs
    assert isinstance(raw, np.ndarray) and raw.shape == (32, 32)
    assert mask.dtype == np.bool_ and mask.shape == (32, 32)
    labels, n = cells
    assert labels.dtype == np.uint32 and n == 3
    assert isinstance(features, dict)
    assert features["area"] == [36.0, 4.0, 36.0]  # raster order of first pixel
    assert "integrated" in features
    assert count.shape == (1, 1) and count[0, 0] == 3.0


def test_filter_and_label_outputs_track_object_count(tmp_path):
    png = tmp_path / "plate.png"
    write_png(png, _squares())
    toml = _write_definition(
        tmp_path,
        f"""
[[nodes]]
name = "raw"
type = "load_image"
params = {{ path = '{png}' }}

[[nodes]]
name = "mask"
type = "otsu"
inputs = ["raw"]

[[nodes]]
name = "cells"
type = "label"
inputs = ["mask"]

[[nodes]]
name = "big"
type = "filter_objects"
params = {{ min_area = 10 }}
inputs = ["cells"]
""",
    )
    _, _, _, kept = tpt_bio.Pipeline.from_toml(str(toml)).run()
    labels, count = kept
    assert count == 2
    assert sorted(np.unique(labels).tolist()) == [0, 1, 2]


def test_fastq_pipeline_returns_records_and_count(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nIIII\n", encoding="utf-8")
    toml = _write_definition(
        tmp_path,
        f"""
[[nodes]]
name = "reads"
type = "read_fastq"
params = {{ path = '{fastq}' }}

[[nodes]]
name = "n"
type = "count_reads"
inputs = ["reads"]
""",
    )
    reads, count = tpt_bio.Pipeline.from_toml(str(toml)).run()
    assert [r.name for r in reads] == [b"r1", b"r2"]
    assert reads[0].seq == b"ACGT"
    assert count[0, 0] == 2.0


def test_frame_csv_export_returns_written_path(tmp_path):
    png = tmp_path / "plate.png"
    write_png(png, _squares())
    out_csv = tmp_path / "features.csv"
    toml = _write_definition(
        tmp_path,
        f"""
[[nodes]]
name = "raw"
type = "load_image"
params = {{ path = '{png}' }}

[[nodes]]
name = "mask"
type = "otsu"
inputs = ["raw"]

[[nodes]]
name = "cells"
type = "label"
inputs = ["mask"]

[[nodes]]
name = "features"
type = "measure"
inputs = ["cells", "raw"]

[[nodes]]
name = "export"
type = "write_frame_csv"
params = {{ path = '{out_csv}' }}
inputs = ["features"]
""",
    )
    *_, written = tpt_bio.Pipeline.from_toml(str(toml)).run()
    assert written == str(out_csv)
    assert out_csv.read_text(encoding="utf-8").startswith("label,area,")


def test_invalid_definition_is_value_error(tmp_path):
    toml = _write_definition(tmp_path, '[[nodes]]\nname = "a"\ntype = "no_such_kind"\n')
    with pytest.raises(ValueError, match="no_such_kind"):
        tpt_bio.Pipeline.from_toml(str(toml))


def test_missing_definition_file_is_os_error(tmp_path):
    with pytest.raises(OSError):
        tpt_bio.Pipeline.from_toml(str(tmp_path / "absent.toml"))


def test_failing_node_is_runtime_error(tmp_path):
    toml = _write_definition(
        tmp_path,
        """
[[nodes]]
name = "data"
type = "table_source"
params = { values = [[0.0, 0.0], [1.0, 2.0]] }

[[nodes]]
name = "reduced"
type = "pca"
params = { components = 5 }
inputs = ["data"]
""",
    )
    with pytest.raises(RuntimeError, match="reduced"):
        tpt_bio.Pipeline.from_toml(str(toml)).run()


def test_table_pipeline_returns_float_array(tmp_path):
    toml = _write_definition(
        tmp_path,
        """
[[nodes]]
name = "data"
type = "table_source"
params = { values = [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]] }

[[nodes]]
name = "scaled"
type = "normalize"
params = { method = "zscore" }
inputs = ["data"]

[[nodes]]
name = "clusters"
type = "kmeans"
params = { k = 2, seed = 1 }
inputs = ["scaled"]
""",
    )
    _, table, labels = tpt_bio.Pipeline.from_toml(str(toml)).run()
    assert table.dtype == np.float64 and table.shape == (4, 2)
    assert labels.dtype == np.uint64
    assert labels[0] == labels[1] != labels[2]
