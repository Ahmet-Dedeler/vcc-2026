#!/usr/bin/env python3
"""Synthetic and frozen-data tests for experiment 108."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp108_run", HERE / "run.py")
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "run.py")
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_gallery_and_axis_exclude_all_target_coordinates() -> None:
    gallery = RUN.frozen_gallery()
    axis = RUN.response_axis(gallery)
    assert len(gallery) == RUN.GALLERY_SIZE
    assert not (set(gallery) & RUN.PANEL_SET)
    assert not (set(axis) & set(gallery))
    assert not (set(axis) & RUN.PANEL_SET)


def test_regulon_is_signed_normalized_and_nonempty() -> None:
    gallery = RUN.frozen_gallery()
    axis = RUN.response_axis(gallery)
    regulon, tf_names, support, audit = RUN.build_regulon(RUN.NETWORK, axis)
    assert len(tf_names) > 50
    assert audit["signed_rows"] > 1_000
    assert support.sum() > 1_000
    norms = np.linalg.norm(regulon, axis=1)
    np.testing.assert_allclose(norms, 1.0, rtol=1e-6, atol=1e-6)
    assert np.any(regulon > 0)
    assert np.any(regulon < 0)


def test_projection_is_in_regulon_row_space() -> None:
    regulon = np.asarray([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    effects = np.asarray([[2.0, -3.0, 9.0]])
    projected = RUN.project_regulon(effects, regulon, ridge=0.0)
    np.testing.assert_allclose(projected, [[2.0, -3.0, 0.0]])


def test_reweight_preserves_sign_zero_support_and_norm() -> None:
    base = np.asarray([[2.0, -1.0, 4.0, 0.0]])
    reconstruction = np.asarray([[1.0, 3.0, 0.0, 9.0]])
    support = np.asarray([True, True, False, True])
    observed = RUN.sign_preserving_reweight(
        base, reconstruction, support, beta=0.1
    )
    np.testing.assert_array_equal(np.sign(observed), np.sign(base))
    np.testing.assert_allclose(np.linalg.norm(observed), np.linalg.norm(base))
    assert observed[0, 2] != 0
    assert observed[0, 3] == 0


def test_direction_blend_restores_base_norm() -> None:
    base = np.asarray([[1.0, 2.0, -1.0]])
    reconstruction = np.asarray([[2.0, -1.0, 4.0]])
    observed = RUN.direction_blend(base, reconstruction, weight=0.1)
    np.testing.assert_allclose(np.linalg.norm(observed), np.linalg.norm(base))


def test_pds_identity_is_one() -> None:
    real = np.eye(4)
    np.testing.assert_array_equal(RUN.pds_per_target(real, real), 1)


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
