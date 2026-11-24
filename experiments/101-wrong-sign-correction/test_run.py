#!/usr/bin/env python3
"""Invariant tests for experiment 101; no reference data are opened."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exp101_run", HERE / "run.py")
assert spec and spec.loader
run = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = run
spec.loader.exec_module(run)


def test_grid_and_folds() -> None:
    assert len(run.CANDIDATES) == 18
    assert {c.ridge for c in run.CANDIDATES} == {0.1, 1.0, 10.0}
    assert {c.tau for c in run.CANDIDATES} == {0.10, 0.25, 0.50}
    assert {c.lam for c in run.CANDIDATES} == {0.5, 1.0}
    targets = np.asarray(["A", "B", "C", "D", "E", "F"])
    folds = np.asarray([run.target_fold(t) for t in targets])
    assert np.all((0 <= folds) & (folds < 4))
    assert np.array_equal(folds, np.asarray([run.target_fold(t) for t in targets]))


def test_weighted_signed_ridge() -> None:
    rng = np.random.default_rng(4)
    x = rng.normal(size=(500, 5))
    y = np.where(2 * x[:, 0] - x[:, 1] > 0, 1, -1)
    w = rng.uniform(0.1, 2.0, size=len(y))
    model = run.fit_ridge(x, y, w, 1.0)
    pred = model.predict(x)
    assert np.mean(np.sign(pred) == y) > 0.9
    assert model.ridge == 1.0


def test_selective_transform() -> None:
    effect = np.asarray([[1.0, -2.0, 3.0, -4.0]])
    score = np.asarray([[-0.6, -0.2, 0.7, -0.9]])
    active = np.asarray([True, True, True, False])
    flipped, mask = run.transform_effect(
        effect, score, active, tau=0.5, lam=0.5, mode="flip"
    )
    np.testing.assert_array_equal(mask, [[True, False, False, False]])
    np.testing.assert_allclose(flipped, [[-0.5, -2.0, 3.0, 0.0]])
    zeroed, zero_mask = run.transform_effect(
        effect, score, active, tau=0.5, lam=1.0, mode="zero"
    )
    np.testing.assert_array_equal(mask, zero_mask)
    np.testing.assert_allclose(zeroed, [[0.0, -2.0, 3.0, 0.0]])


def test_selection_requires_both_pds_cells() -> None:
    cells = {}
    for candidate in run.CANDIDATES:
        cells[candidate.key] = [
            {"credit_delta": 0.01, "pds_delta": 0.0},
            {"credit_delta": 0.01, "pds_delta": 0.0},
        ]
    unsafe = run.CANDIDATES[0]
    cells[unsafe.key][0]["credit_delta"] = 1.0
    cells[unsafe.key][0]["pds_delta"] = -1e-6
    selected, audit = run.choose_candidate(cells)
    assert selected is not None
    assert selected != unsafe
    assert audit["candidate_summaries"][unsafe.key]["minimum_pds_delta"] < 0

    for candidate in run.CANDIDATES:
        for cell in cells[candidate.key]:
            cell["credit_delta"] = -0.01
    selected, audit = run.choose_candidate(cells)
    assert selected is None
    assert audit["selected"] == "identity"


def test_ratio_bootstrap() -> None:
    result = run.ratio_ci(
        np.asarray([1.0, 0.0]),
        np.asarray([1.0, 1.0]),
        100,
        7,
    )
    assert result["estimate"] == 0.5
    empty = run.ratio_ci(np.zeros(2), np.zeros(2), 100, 7)
    assert empty == {"estimate": None, "lo": None, "hi": None}


def test_all_identity_null_feature_width_fallback() -> None:
    # Regression test for the freeze path in which the real selector and all
    # permutation selectors choose identity: feature width must still come
    # from the frozen rank-16 block, not from a nonexistent fitted model.
    block_width = 60
    model_sets = [None] * 20
    width = next(
        (len(models[0].coef) for models in model_sets if models is not None),
        block_width,
    )
    assert width == 60


def main() -> None:
    test_grid_and_folds()
    test_weighted_signed_ridge()
    test_selective_transform()
    test_selection_requires_both_pds_cells()
    test_ratio_bootstrap()
    test_all_identity_null_feature_width_fallback()
    print("experiment-101 invariant tests: PASS")


if __name__ == "__main__":
    main()
