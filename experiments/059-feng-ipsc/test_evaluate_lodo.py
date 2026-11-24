"""Unit tests for low-memory NPY access and Feng blending."""

from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "evaluate_feng_lodo", HERE / "evaluate_lodo.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_npy_pread_matches_numpy() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "matrix.npy"
        expected = np.arange(60, dtype=np.float32).reshape(6, 10)
        np.save(path, expected)
        layout = MODULE.NpyLayout.inspect(path)
        descriptor = os.open(path, os.O_RDONLY)
        try:
            np.testing.assert_array_equal(
                layout.read_block(descriptor, 2, 5), expected[2:5]
            )
            np.testing.assert_array_equal(layout.read_row(descriptor, 4), expected[4])
        finally:
            os.close(descriptor)


def test_blend_preserves_unavailable_coordinates() -> None:
    baseline = np.arange(20, dtype=np.float64).reshape(4, 5) + 1
    feng = np.ones((4, 5), dtype=np.float64)
    available = np.zeros((4, 5), dtype=np.bool_)
    available[:2, :3] = True
    active = np.asarray([True, False, False, False])
    candidate, scaling = MODULE.blend_feng(baseline, feng, available, active, 0.1)
    assert scaling["feng_scale"] > 0
    np.testing.assert_array_equal(candidate[1:], baseline[1:])
    np.testing.assert_array_equal(candidate[0, 3:], baseline[0, 3:])
    assert not np.array_equal(candidate[0, :3], baseline[0, :3])
