from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("dual_moment", HERE / "dual_moment.py")
assert SPEC and SPEC.loader
DM = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DM
SPEC.loader.exec_module(DM)


def fixture() -> tuple[sp.csr_matrix, sp.csr_matrix]:
    rng = np.random.default_rng(147)
    incumbent = rng.poisson(
        np.array([10.0, 4.0, 1.0, 8.0, 2.0, 5.0])[None, :],
        size=(8, 6),
    ).astype(np.int64)
    direction = incumbent.copy()
    direction[:, 0] += np.arange(8)
    direction[:, 1] = np.maximum(0, direction[:, 1] - np.arange(8) // 2)
    incumbent += 1
    direction += 1
    return sp.csr_matrix(incumbent), sp.csr_matrix(direction)


def test_exact_margins_pseudobulk_and_count_contract() -> None:
    incumbent, direction = fixture()
    emitted, diagnostics = DM.emit_dual_moment(
        incumbent,
        direction,
        blend_lambda=0.5,
        depth_contrast=0.5,
        quasi_spread=0.25,
        seed=14700,
    )
    np.testing.assert_array_equal(
        np.asarray(emitted.sum(axis=0)).ravel(),
        np.asarray(incumbent.sum(axis=0)).ravel(),
    )
    assert np.all(np.asarray(emitted.sum(axis=1)).ravel() > 0)
    np.testing.assert_array_equal(
        DM.scorer_pseudobulk(emitted),
        DM.scorer_pseudobulk(incumbent),
    )
    assert np.issubdtype(emitted.dtype, np.integer)
    assert emitted.data.min(initial=0) >= 0
    assert diagnostics.max_row_total <= 1_000_000


def test_bounded_projection_is_exact() -> None:
    raw = np.array([-5.0, 2.2, 9.1, 100.0])
    caps = np.array([0, 5, 10, 20], dtype=np.int64)
    projected, integer = DM._project_box_sum(
        raw,
        caps,
        17,
        np.random.default_rng(1),
    )
    assert np.isclose(projected.sum(), 17.0, atol=1e-8)
    assert integer.sum() == 17
    assert np.all(integer >= 0)
    assert np.all(integer <= caps)


def test_deterministic_and_directionally_effective() -> None:
    incumbent, direction = fixture()
    first, first_diag = DM.emit_dual_moment(
        incumbent,
        direction,
        blend_lambda=1.0,
        seed=101,
    )
    second, second_diag = DM.emit_dual_moment(
        incumbent,
        direction,
        blend_lambda=1.0,
        seed=101,
    )
    assert (first != second).nnz == 0
    assert first_diag == second_diag

    base = DM.arithmetic_mean_cpm(incumbent)
    target = DM.arithmetic_mean_cpm(direction)
    achieved = DM.arithmetic_mean_cpm(first)
    assert np.abs(achieved - target).sum() < np.abs(base - target).sum()


def test_row_margin_shape_and_depth_contrast() -> None:
    rows, low = DM.two_depth_row_margins(80_003, 400, 0.5, seed=4)
    assert rows.sum() == 80_003
    assert low.sum() == 200
    assert rows[low].max() - rows[low].min() <= 1
    assert rows[~low].max() - rows[~low].min() <= 1
    assert rows[low].mean() < rows[~low].mean()
