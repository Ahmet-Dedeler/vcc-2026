from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("context_adaptation", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_pds_identity_is_perfect() -> None:
    rng = np.random.default_rng(1)
    truth = MODULE.row_normalize(rng.normal(size=(30, 50)))
    scores = truth @ truth.T
    assert MODULE.pds_from_scores(scores) == 1.0


def test_diagonal_map_recovers_held_out_direction() -> None:
    rng = np.random.default_rng(2)
    x = MODULE.row_normalize(rng.normal(size=(180, 80)))
    scale = np.linspace(0.35, 1.8, x.shape[1])
    y = MODULE.row_normalize(x * scale)
    candidate = MODULE.Candidate(
        "diagonal",
        ridge=0.3,
        eta=1.0,
        clip="positive2",
    )
    model = MODULE.fit_models(x[:120], y[:120], [candidate], seed=3)[candidate.key]
    baseline_cosine = float(np.mean(np.einsum("ij,ij->i", x[120:], y[120:])))
    adapted_vectors = MODULE.row_normalize(
        MODULE.apply_local(x[120:], model)
    )
    adapted_cosine = float(
        np.mean(np.einsum("ij,ij->i", adapted_vectors, y[120:]))
    )
    assert adapted_cosine > baseline_cosine


def test_candidate_grid_is_unique_and_low_rank() -> None:
    candidates = MODULE.candidate_grid()
    keys = [candidate.key for candidate in candidates]
    assert len(keys) == len(set(keys))
    assert max(candidate.rank for candidate in candidates) <= 32
    assert max(candidate.residual_rank for candidate in candidates) <= 16


def test_fast_scores_match_explicit_adapted_vectors() -> None:
    rng = np.random.default_rng(4)
    x = MODULE.row_normalize(rng.normal(size=(90, 60)))
    y = MODULE.row_normalize(rng.normal(size=(90, 60)))
    candidate = MODULE.Candidate(
        "pca_ridge",
        ridge=0.3,
        rank=8,
        eta=0.5,
    )
    basis, _ = np.linalg.qr(rng.normal(size=(60, 8)))
    mapping = 0.1 * rng.normal(size=(8, 8))
    model = MODULE.Fitted(
        candidate,
        basis=basis,
        mapping=mapping,
    )
    pair = MODULE.Pair(
        x[60:],
        y[60:],
        x[60:],
        y[60:],
        np.ones(30),
    )
    # Apple's BLAS occasionally leaves floating-status flags set under pytest
    # even though the resulting small matrices are finite; validate the values
    # directly and suppress those spurious status warnings here.
    with np.errstate(all="ignore"):
        fast = MODULE.model_scores(pair, model)
        explicit = MODULE.apply_local(x[60:], model) @ y[60:].T
    assert np.isfinite(fast).all()
    assert np.isfinite(explicit).all()
    np.testing.assert_allclose(fast, explicit, rtol=1e-10, atol=1e-10)
