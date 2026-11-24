#!/usr/bin/env python3
"""Synthetic tests for the TeloHAEC direct-source pre-gate."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("source_gate", HERE / "source_gate.py")
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "source_gate.py")
GATE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GATE
SPEC.loader.exec_module(GATE)


def test_gallery_is_fixed_nonpanel_and_complete() -> None:
    gallery = GATE.frozen_gallery()
    assert len(gallery) == GATE.GALLERY_SIZE
    assert len(set(gallery)) == len(gallery)
    assert not (set(gallery) & GATE.PANEL_SET)


def test_response_axis_excludes_gallery_and_panel() -> None:
    gallery = GATE.frozen_gallery()
    axis = GATE.common_gene_axis(gallery)
    assert not (set(axis) & set(gallery))
    assert not (set(axis) & GATE.PANEL_SET)


def test_equal_source_blend_is_exact_mean() -> None:
    left = np.asarray([[1.0, -1.0], [3.0, 5.0]], dtype=np.float32)
    right = np.asarray([[3.0, 1.0], [1.0, -1.0]], dtype=np.float32)
    observed = GATE.equal_source_blend(left, right)
    np.testing.assert_array_equal(observed, (left + right) / 2)


def test_pds_identity_is_one() -> None:
    real = np.eye(4)
    np.testing.assert_array_equal(GATE.pds_per_target(real, real), 1)


def test_permutation_seeds_are_deterministic() -> None:
    first = np.random.default_rng(GATE.stable_seed(GATE.EXPERIMENT, "null", 3))
    second = np.random.default_rng(GATE.stable_seed(GATE.EXPERIMENT, "null", 3))
    np.testing.assert_array_equal(first.permutation(100), second.permutation(100))


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
