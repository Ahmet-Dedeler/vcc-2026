from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("source_only_subspace", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_pds_identity_is_perfect() -> None:
    rng = np.random.default_rng(1)
    truth = MODULE.row_normalize(rng.normal(size=(30, 50)))
    assert MODULE.pds_from_scores(truth @ truth.T) == 1.0


def test_fast_spectral_scores_match_explicit_transform() -> None:
    rng = np.random.default_rng(2)
    q = MODULE.row_normalize(rng.normal(size=(20, 40)))
    t = MODULE.row_normalize(rng.normal(size=(20, 40)))
    basis, _ = np.linalg.qr(rng.normal(size=(40, 8)))
    random_basis, _ = np.linalg.qr(rng.normal(size=(40, 8)))
    subspace = MODULE.Subspace(
        basis.astype(np.float32),
        np.linspace(5.0, 1.0, 8),
        random_basis.astype(np.float32),
    )
    candidate = MODULE.Candidate(
        "source", "spectral", rank=8, gamma=0.65, ridge=0.1
    )
    pair = MODULE.Pair(q, t, q, t)
    _, _, matrices = MODULE.score_candidates(
        pair, subspace, [candidate], matrix_keys={candidate.key}
    )
    weights = MODULE.candidate_weights(candidate, subspace.source_evals)
    with np.errstate(all="ignore"):
        explicit_q = q + ((q @ basis) * (weights - 1.0)) @ basis.T
        explicit_scores = explicit_q @ t.T
    assert np.isfinite(matrices[candidate.key]).all()
    assert np.isfinite(explicit_scores).all()
    np.testing.assert_allclose(
        matrices[candidate.key],
        explicit_scores,
        rtol=1e-6,
        atol=1e-7,
    )


def test_covariance_subspace_is_row_permutation_invariant() -> None:
    rng = np.random.default_rng(3)
    x = rng.normal(size=(180, 80)).astype(np.float32)
    first = MODULE.fit_subspace(x, max_rank=24, seed=9)
    second = MODULE.fit_subspace(
        x[rng.permutation(len(x))], max_rank=24, seed=9
    )
    with np.errstate(all="ignore"):
        op_first = (
            first.source_basis
            @ np.diag(first.source_evals)
            @ first.source_basis.T
        )
        op_second = (
            second.source_basis
            @ np.diag(second.source_evals)
            @ second.source_basis.T
        )
    assert np.isfinite(op_first).all()
    assert np.isfinite(op_second).all()
    np.testing.assert_allclose(op_first, op_second, rtol=2e-4, atol=2e-5)


def test_selection_never_selects_random_null() -> None:
    source, null = MODULE.build_candidates()
    source_summary = {
        candidate.key: {
            "delta": 0.01,
            "positive_folds": 10,
            "donor_delta": {"a": 0.01, "b": 0.01},
        }
        for candidate in source
    }
    selected = MODULE.select_candidate(
        source_summary,
        source,
        min_gain=0.001,
        min_positive_folds=6,
        tie_tolerance=5e-4,
        require_each_donor=True,
    )
    assert selected in {candidate.key for candidate in source}
    assert selected not in {candidate.key for candidate in null}


def test_train_and_outer_design_is_disjoint() -> None:
    _, train_perts, outer_perts, _ = MODULE.validate_design_files()
    assert len(train_perts) == 600
    assert len(outer_perts) == 259
    assert not set(train_perts) & set(outer_perts)
