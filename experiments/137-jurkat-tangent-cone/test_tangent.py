from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp137_tangent", HERE / "tangent.py")
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def angle_degrees(left: np.ndarray, right: np.ndarray) -> float:
    cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def test_zero_cap_is_bit_exact() -> None:
    base = np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    source = np.asarray([[2.0, -1.0], [-4.0, 3.0]], dtype=np.float64)
    candidate, effective, angle = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True, True]),
        np.asarray([True, True]),
        0.0,
    )
    assert np.array_equal(candidate.view(np.uint8), base.view(np.uint8))
    assert not effective.any()
    assert not angle.any()


def test_rotation_hits_cap_and_preserves_norm() -> None:
    base = np.asarray([[1.0, 0.0]], dtype=np.float64)
    source = np.asarray([[0.0, 2.0]], dtype=np.float64)
    candidate, effective, angle = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True]),
        np.asarray([True, True]),
        3.0,
    )
    assert effective.tolist() == [True]
    assert np.allclose(np.linalg.norm(candidate, axis=1), np.linalg.norm(base, axis=1))
    assert np.isclose(angle[0], 3.0)
    assert np.isclose(angle_degrees(base[0], candidate[0]), 3.0)


def test_short_source_arc_is_not_overshot() -> None:
    theta = np.radians(1.0)
    base = np.asarray([[1.0, 0.0]], dtype=np.float64)
    source = np.asarray([[np.cos(theta), np.sin(theta)]], dtype=np.float64)
    candidate, effective, angle = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True]),
        np.asarray([True, True]),
        6.0,
    )
    assert effective.tolist() == [True]
    assert np.isclose(angle[0], 1.0)
    assert np.allclose(candidate, source)


def test_fallback_rows_and_coordinates_are_bit_exact() -> None:
    base = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64)
    source = np.asarray([[0.0, 3.0, 8.0], [1.0, 0.0, 7.0]], dtype=np.float64)
    candidate, effective, _ = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True, False]),
        np.asarray([True, True, False]),
        3.0,
    )
    assert effective.tolist() == [True, False]
    assert np.array_equal(candidate[1].view(np.uint8), base[1].view(np.uint8))
    assert np.array_equal(
        np.ascontiguousarray(candidate[:, 2]).view(np.uint8),
        np.ascontiguousarray(base[:, 2]).view(np.uint8),
    )


def test_collinear_antipodal_and_zero_rows_fall_back() -> None:
    base = np.asarray(
        [[1.0, 0.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64
    )
    source = np.asarray(
        [[2.0, 0.0], [-2.0, 0.0], [0.0, 1.0]], dtype=np.float64
    )
    candidate, effective, angle = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True, True, True]),
        np.asarray([True, True]),
        3.0,
    )
    assert not effective.any()
    assert not angle.any()
    assert np.array_equal(candidate.view(np.uint8), base.view(np.uint8))


def test_near_parallel_angle_uses_robust_atan2() -> None:
    epsilon = 1e-10
    base = np.asarray([[1.0, 0.0]], dtype=np.float64)
    source = np.asarray([[1.0, epsilon]], dtype=np.float64)
    candidate, effective, angle = M.spherical_tangent_candidate(
        base,
        source,
        np.asarray([True]),
        np.asarray([True, True]),
        6.0,
    )
    assert effective.tolist() == [True]
    assert angle[0] > 0.0
    assert np.allclose(candidate / np.linalg.norm(candidate), source / np.linalg.norm(source))


def test_pds_midrank_and_zero_prediction() -> None:
    truth = np.eye(3, dtype=np.float64)
    assert np.array_equal(M.pds_per_target(truth, truth), np.ones(3))
    zero = np.zeros_like(truth)
    assert np.array_equal(M.pds_per_target(zero, truth), np.full(3, 0.5))


def test_nested_selection_uses_smallest_cap_on_tie() -> None:
    scores = {
        cap: tuple(np.full(63, 0.5 + (0.01 if cap == 3.0 else 0.0)) for _ in range(5))
        for cap in M.CAP_DEGREES
    }
    selected, held = M.select_cap(scores)
    assert selected == (3.0,) * 5
    assert np.allclose(held, 0.51)
    tied = {
        cap: tuple(np.full(63, 0.5) for _ in range(5)) for cap in M.CAP_DEGREES
    }
    selected, held = M.select_cap(tied)
    assert selected == (0.0,) * 5
    assert np.allclose(held, 0.5)


def test_final_cap_rejects_nonfinite_or_wrong_fold_shape() -> None:
    valid = {
        cap: tuple(np.full(63, 0.5) for _ in range(5)) for cap in M.CAP_DEGREES
    }
    assert M.final_cap(valid) == 0.0
    broken = dict(valid)
    broken[1.5] = (*broken[1.5][:-1], np.full(62, 0.5))
    try:
        M.final_cap(broken)
    except M.GeometryError:
        pass
    else:
        raise AssertionError("wrong fold shape was accepted")
    nonfinite = dict(valid)
    bad = np.full(63, 0.5)
    bad[0] = np.nan
    nonfinite[3.0] = (bad, *nonfinite[3.0][1:])
    try:
        M.final_cap(nonfinite)
    except M.GeometryError:
        pass
    else:
        raise AssertionError("nonfinite fold score was accepted")


if __name__ == "__main__":
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} tests passed")
