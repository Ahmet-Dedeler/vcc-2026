from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp057", HERE / "run.py")
RUN = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_pds_uses_complete_gallery() -> None:
    truth = np.eye(4, dtype=np.float64)
    prediction = truth.copy()
    values = RUN.pds_per_query(prediction, truth)
    np.testing.assert_allclose(values, 1.0)
    prediction[[0, 1]] = prediction[[1, 0]]
    values = RUN.pds_per_query(prediction, truth)
    assert values[0] < 1.0
    assert values[1] < 1.0


def test_query_folds_partition_queries_but_not_gallery() -> None:
    targets = [f"P{i}" for i in range(10)]
    rows = RUN.query_fold_summary(np.arange(10, dtype=float), targets, 5)
    assert sum(row["n_queries"] for row in rows) == 10
    assert {row["gallery_size"] for row in rows} == {10}


def test_fixed_blend_preserves_missing_auxiliary_genes() -> None:
    baseline = np.asarray([[1.0, 2.0, 3.0]])
    auxiliary = np.asarray([[4.0, 5.0, 6.0]])
    available = np.asarray([True, False, True])
    candidate = baseline.copy()
    candidate[:, available] = (1.0 - RUN.KOLF_FRACTION) * baseline[
        :, available
    ] + RUN.KOLF_FRACTION * auxiliary[:, available]
    assert candidate[0, 1] == baseline[0, 1]
    np.testing.assert_allclose(candidate[0, [0, 2]], [1.3, 3.3])


def test_whitening_positive_scalar_cannot_change_pds() -> None:
    rng = np.random.default_rng(57)
    prediction = rng.normal(size=(12, 8))
    truth = rng.normal(size=(12, 8))
    expected = RUN.pds_per_query(prediction, truth)
    observed = RUN.pds_per_query(prediction * 9.7, truth)
    np.testing.assert_allclose(observed, expected)


def test_axis_hash_is_order_sensitive() -> None:
    assert RUN.sha256_lines(["A", "B"]) != RUN.sha256_lines(["B", "A"])
