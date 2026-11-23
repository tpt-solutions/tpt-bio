"""Vision binding tests: thresholding, labelling and watershed on NumPy arrays."""

import numpy as np
import pytest

tpt_bio = pytest.importorskip("tpt_bio")


def test_otsu_separates_bimodal_image():
    img = np.array([[0, 0, 10, 10], [0, 0, 10, 10]], dtype=np.uint8)
    t = tpt_bio.otsu_threshold(img)
    assert t is not None
    assert 0 <= t < 10


def test_otsu_constant_image_is_none():
    assert tpt_bio.otsu_threshold(np.zeros((3, 3), dtype=np.float64)) is None


def test_otsu_rejects_bad_dtype():
    with pytest.raises(ValueError):
        tpt_bio.otsu_threshold(np.zeros((3, 3), dtype=np.bool_))


def test_apply_threshold_returns_bool_mask():
    img = np.array([[1.0, 5.0], [9.0, 2.0]])
    mask = tpt_bio.apply_threshold(img, 4.0)
    assert mask.dtype == np.bool_
    assert mask.tolist() == [[False, True], [True, False]]


def test_label_components_counts_blobs():
    mask = np.array(
        [
            [True, True, False, False],
            [False, False, False, True],
            [False, False, False, True],
        ]
    )
    labels, count = tpt_bio.label_components(mask, 4)
    assert count == 2
    assert labels.dtype == np.uint32
    assert labels.shape == (3, 4)
    assert labels[0, 0] == labels[0, 1] != 0
    assert labels[0, 2] == 0
    assert labels[1, 3] == labels[2, 3] != labels[0, 0]


def test_label_connectivity_matters_for_diagonals():
    mask = np.array([[True, False], [False, True]])
    _, four = tpt_bio.label_components(mask, 4)
    _, eight = tpt_bio.label_components(mask, 8)
    assert four == 2
    assert eight == 1


def test_label_components_rejects_bad_connectivity():
    with pytest.raises(ValueError):
        tpt_bio.label_components(np.ones((2, 2), dtype=np.bool_), 6)


def test_label_components_default_is_eight_connected():
    mask = np.array([[True, False], [False, True]])
    _, count = tpt_bio.label_components(mask)
    assert count == 1


def test_watershed_splits_between_markers():
    # A valley at column 2 separates two basins; each marker floods its side.
    img = np.array([[0.0, 0.0, 5.0, 1.0, 1.0]])
    markers = np.array([[1, 0, 0, 0, 2]], dtype=np.uint32)
    out = tpt_bio.watershed_2d(img, markers)
    assert out.dtype == np.uint32
    assert out[0, 0] == 1
    assert out[0, 4] == 2


def test_watershed_mask_restricts_flooding():
    img = np.zeros((1, 4))
    markers = np.array([[1, 0, 0, 0]], dtype=np.uint32)
    mask = np.array([[True, True, False, False]])
    out = tpt_bio.watershed_2d(img, markers, mask)
    assert out.tolist() == [[1, 1, 0, 0]]


def test_non_contiguous_input_is_accepted():
    base = np.arange(16, dtype=np.float64).reshape(4, 4)
    view = base[:, ::2]  # strided view
    assert tpt_bio.otsu_threshold(view) is not None
