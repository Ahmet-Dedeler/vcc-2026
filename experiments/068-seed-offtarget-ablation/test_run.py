#!/usr/bin/env python3
"""Small deterministic checks for experiment 068's bounded primitives."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp068_run", HERE / "run.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import run.py")
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_stored_npz_transform() -> None:
    path = HERE / ".synthetic_operator.npz"
    rng = np.random.default_rng(68)
    raw = rng.normal(size=(9, 9))
    covariance = raw.T @ raw / len(raw)
    evals, evecs = np.linalg.eigh(covariance)
    evals = np.maximum(evals, evals.max() * 1e-12)
    np.savez(path, evals=evals, evecs=np.asarray(evecs, dtype=np.float32))
    try:
        stored = RUN.stored_npy_member(path, "evecs.npy")
        operator = RUN.Operator(
            name="synthetic",
            path=path,
            evals=evals,
            evecs=stored,
            wcols=np.arange(9, dtype=np.int64),
            preserve_scale=1.0,
            metadata={},
        )
        values = rng.normal(size=(5, 9))
        observed = RUN.transform_rows(values, operator, block_rows=3)
        weights = evals ** (-RUN.WHITEN_GAMMA / 2.0)
        expected = ((values @ evecs.astype(np.float32)) * weights) @ evecs.astype(
            np.float32
        ).T
        np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)
    finally:
        path.unlink(missing_ok=True)


def test_rank_query() -> None:
    held = RUN.PreparedMain(
        name="held",
        effects=np.asarray(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.5, 0.5, 0.0],
            ],
            dtype=np.float64,
        ),
        available=np.ones(3, dtype=bool),
        mean_on_axis=np.zeros(3, dtype=np.float64),
        raw_median_norm=1.0,
        mapped_gene_count=3,
        duplicate_symbols=(),
    )
    result = RUN.rank_query(
        np.asarray([1.0, 0.0, 0.0]),
        held,
        np.arange(3, dtype=np.int64),
        0,
        np.arange(3, dtype=np.int64),
        np.arange(3, dtype=np.int64),
    )
    assert result["rank_1based"] == 1.0
    assert result["pds"] == 1.0


def test_fallback_ablation() -> None:
    effect = np.asarray([[1.0, 2.0], [3.0, 4.0]])
    available = np.ones(2, dtype=bool)
    companion = RUN.PreparedMain(
        name="companion",
        effects=effect.copy(),
        available=available.copy(),
        mean_on_axis=np.asarray([2.0, 1.0]),
        raw_median_norm=1.0,
        mapped_gene_count=2,
        duplicate_symbols=(),
    )
    result = RUN.build_baseline_rows(
        [0],
        effect,
        available,
        np.asarray([1.0, 2.0]),
        companion,
        removed=frozenset({"k562_slot", "companion"}),
    )
    expected = RUN.EFFECT_SCALE * np.asarray([1.5, 1.5]) / np.linalg.norm(
        [1.5, 1.5]
    )
    np.testing.assert_allclose(result[0], expected)


if __name__ == "__main__":
    test_stored_npz_transform()
    test_rank_query()
    test_fallback_ablation()
    print("3 synthetic checks passed")
