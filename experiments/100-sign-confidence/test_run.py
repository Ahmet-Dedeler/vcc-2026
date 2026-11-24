#!/usr/bin/env python3
"""Small invariant tests for experiment 100; no real evaluation is opened."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exp100_run", HERE / "run.py")
assert spec and spec.loader
run = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = run
spec.loader.exec_module(run)


def main() -> None:
    targets = np.asarray(["A", "B", "C", "D", "E", "F", "G", "H"])
    folds = np.asarray([run.target_fold(t) for t in targets])
    assert np.all((0 <= folds) & (folds < run.N_FOLDS))
    assert all(run.target_fold(t) == run.target_fold(t) for t in targets)

    rng = np.random.default_rng(0)
    effect = rng.normal(size=(3, 20))
    probability = rng.uniform(0.1, 0.9, size=(3, 20))
    active = np.ones(20, dtype=bool)
    active[[2, 7]] = False
    effect[:, ~active] = 0
    transformed = run.transform_effect(effect, probability, active)
    assert np.all(np.sign(transformed[:, active]) == np.sign(effect[:, active]))
    assert np.all(transformed[:, ~active] == 0)
    assert np.allclose(
        np.linalg.norm(transformed[:, active], axis=1),
        np.linalg.norm(effect[:, active], axis=1),
        atol=1e-10,
    )

    pred = np.eye(5)
    real = np.eye(5)
    pds = run.pds_per_target(pred, real, np.ones(5, dtype=bool))
    assert np.allclose(pds, 1.0)

    correct = np.asarray([1, 0, 1, 0, 1], dtype=float)
    confidence = np.asarray([5, 4, 3, 2, 1], dtype=float)
    weight = np.ones(5)
    assert run.weighted_top_precision(correct, confidence, weight, 0.4) == 0.5

    print("experiment-100 invariant tests: PASS")


if __name__ == "__main__":
    main()
