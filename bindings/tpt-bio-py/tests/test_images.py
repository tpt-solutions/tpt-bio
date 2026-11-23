"""Tests for load_image and object_features. PNG fixtures are written by the small
encoder below, so the tests need no image library and no binary files."""

import struct
import zlib

import numpy as np
import pytest

tpt_bio = pytest.importorskip("tpt_bio")


def _chunk(kind: bytes, body: bytes) -> bytes:
    crc = zlib.crc32(kind + body) & 0xFFFFFFFF
    return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", crc)


def write_png(path, pixels: np.ndarray, bit_depth: int = 8) -> None:
    """Writes an (h, w) grey or (h, w, 3) RGB array as an uncompressed-filter PNG."""
    if pixels.ndim == 2:
        colour_type, channels = 0, 1
    else:
        colour_type, channels = 2, 3
    h, w = pixels.shape[:2]
    rows = pixels.reshape(h, w * channels)
    if bit_depth == 8:
        raw = rows.astype(np.uint8)
    else:
        raw = rows.astype(">u2")
    scan = b"".join(b"\x00" + raw[y].tobytes() for y in range(h))
    header = struct.pack(">IIBBBBB", w, h, bit_depth, colour_type, 0, 0, 0)
    data = b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", header)
    data += _chunk(b"IDAT", zlib.compress(scan)) + _chunk(b"IEND", b"")
    path.write_bytes(data)


def test_load_grey_png_keeps_uint8(tmp_path):
    pixels = np.arange(12, dtype=np.uint8).reshape(3, 4) * 10
    path = tmp_path / "grey.png"
    write_png(path, pixels)
    arr, channels = tpt_bio.load_image(str(path))
    assert channels == 1
    assert arr.dtype == np.uint8
    assert arr.shape == (3, 4)
    np.testing.assert_array_equal(arr, pixels)


def test_load_16_bit_png_gives_uint16(tmp_path):
    pixels = np.array([[0, 1000], [40000, 65535]], dtype=np.uint16)
    path = tmp_path / "deep.png"
    write_png(path, pixels, bit_depth=16)
    arr, channels = tpt_bio.load_image(str(path))
    assert channels == 1
    assert arr.dtype == np.uint16
    np.testing.assert_array_equal(arr, pixels)


def test_rgb_png_reports_three_channels_and_selects_one(tmp_path):
    rgb = np.zeros((2, 2, 3), dtype=np.uint8)
    rgb[..., 0] = 1
    rgb[..., 1] = 2
    rgb[..., 2] = 3
    path = tmp_path / "rgb.png"
    write_png(path, rgb)
    blue, channels = tpt_bio.load_image(str(path), 2)
    assert channels == 3
    np.testing.assert_array_equal(blue, np.full((2, 2), 3, dtype=np.uint8))


def test_channel_out_of_range_is_value_error(tmp_path):
    path = tmp_path / "grey.png"
    write_png(path, np.zeros((2, 2), dtype=np.uint8))
    with pytest.raises(ValueError):
        tpt_bio.load_image(str(path), 1)


def test_unsupported_extension_is_value_error(tmp_path):
    path = tmp_path / "image.jpg"
    path.write_bytes(b"not an image")
    with pytest.raises(ValueError):
        tpt_bio.load_image(str(path))


def test_missing_file_is_os_error(tmp_path):
    with pytest.raises(OSError):
        tpt_bio.load_image(str(tmp_path / "absent.png"))


def test_object_features_shape_columns():
    labels = np.zeros((5, 6), dtype=np.uint32)
    labels[0:2, 0:2] = 1  # area 4, a square
    labels[2:5, 3:6] = 7  # area 9, label 7 renumbered to 2
    feats = tpt_bio.object_features(labels)
    assert set(feats) == {"label", "area", "perimeter", "eccentricity", "solidity"}
    assert feats["area"] == [4.0, 9.0]
    assert len(feats["label"]) == 2


def test_object_features_adds_intensity_columns():
    labels = np.zeros((4, 4), dtype=np.uint32)
    labels[0:2, 0:2] = 1
    intensity = np.zeros((4, 4), dtype=np.float64)
    intensity[0:2, 0:2] = 10.0
    feats = tpt_bio.object_features(labels, intensity)
    assert feats["mean"] == [10.0]
    assert feats["integrated"] == [40.0]
    assert feats["median"] == [10.0]


def test_object_features_rejects_shape_mismatch():
    labels = np.zeros((2, 2), dtype=np.uint32)
    labels[0, 0] = 1
    with pytest.raises(ValueError):
        tpt_bio.object_features(labels, np.zeros((3, 3), dtype=np.float64))
