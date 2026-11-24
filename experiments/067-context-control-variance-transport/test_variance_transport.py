from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp067", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_stored_npz_member_direct_roundtrip() -> None:
    rng = np.random.default_rng(1)
    expected = rng.standard_normal((11, 7)).astype(np.float32)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "operator.npz"
        np.savez(path, evecs=expected, evals=np.arange(7, dtype=np.float64))
        layout = RUN.stored_npy_layout(path, "evecs.npy")
        assert layout.shape == expected.shape
        actual = RUN.load_stored_npy_member(path, "evecs.npy")
        np.testing.assert_array_equal(actual, expected)


def test_blockwise_operator_matches_dense() -> None:
    rng = np.random.default_rng(2)
    basis = rng.standard_normal((9, 9))
    covariance = basis.T @ basis / len(basis)
    evals, evecs = np.linalg.eigh(covariance)
    evals = np.maximum(evals, evals.max() * 1e-12)
    x = rng.standard_normal((5, 9)).astype(np.float32)
    expected = (
        (x @ evecs) * np.power(evals, -RUN.WHITEN_GAMMA / 2.0)
    ) @ evecs.T
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "operator.npz"
        np.savez(path, evecs=evecs.astype(np.float32))
        layout = RUN.stored_npy_layout(path, "evecs.npy")
        actual = RUN.apply_operator_rows(
            x,
            evals,
            layout,
            block_rows=3,
            max_rss_gib=0,
            label="synthetic",
        )
    np.testing.assert_allclose(actual, expected, rtol=2e-5, atol=2e-5)


def test_direct_npy_row_blocks_roundtrip() -> None:
    expected = np.arange(35, dtype=np.float32).reshape(5, 7)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "matrix.npy"
        np.save(path, expected)
        blocks = list(RUN.iter_npy_row_blocks(path, chunk_rows=2))
    assert [offset for offset, _ in blocks] == [0, 2, 4]
    np.testing.assert_array_equal(
        np.vstack([block for _, block in blocks]),
        expected,
    )


def test_pre_post_identity_and_cosine_score_algebra() -> None:
    rng = np.random.default_rng(3)
    q = rng.standard_normal((8, 12)).astype(np.float32)
    t = RUN.row_normalize(rng.standard_normal((8, 12)).astype(np.float32))
    diagonal = np.linspace(0.75, 1.25, 12)
    qw = q * diagonal
    tw = t * diagonal
    baseline = np.einsum("ik,jk->ij", qw, t)
    reflected = np.einsum("ik,jk->ij", q, tw)
    np.testing.assert_allclose(baseline, reflected, rtol=1e-9, atol=1e-9)
    identity = np.ones(q.shape[1], dtype=np.float32)
    np.testing.assert_allclose(
        np.einsum("ik,jk->ij", q * identity, tw), baseline
    )
    np.testing.assert_allclose(
        np.einsum("ik,jk->ij", qw * identity, t), baseline
    )


def test_ratio_is_clipped_and_missing_values_are_identity() -> None:
    columns = np.arange(5)
    controls = {
        "target_a_sd": np.array([4.0, 0.1, np.nan, 2.0, 1.0]),
        "source_a_sd": np.array([1.0, 2.0, 3.0, np.nan, 1.0]),
    }
    ratio, metadata = RUN.variance_ratio(
        controls,
        "target",
        ("source",),
        "sd",
        "a",
        columns,
    )
    np.testing.assert_allclose(ratio, [2.0, 0.5, 1.0, 1.0, 1.0])
    assert metadata["valid_genes"] == 3
    assert metadata["identity_genes"] == 2
    assert metadata["clipped_low"] == 1
    assert metadata["clipped_high"] == 1


def test_source_only_selection_can_choose_identity() -> None:
    candidates = RUN.candidate_grid()
    inner = []
    for _ in range(2):
        rows = {}
        for candidate in candidates:
            delta = 0.0 if candidate.eta == 0 else -0.01
            rows[candidate.key] = {"delta": delta}
        inner.append(rows)
    selected = RUN.choose_candidate(inner, candidates, tie_tolerance=1e-6)
    assert selected["candidate"] == "identity:eta=0"


def test_control_scale_estimators_are_finite() -> None:
    rng = np.random.default_rng(5)
    sample = rng.poisson(2.0, size=(128, 17)).astype(np.float32)
    stats = RUN.robust_scales(sample, gene_block=4)
    assert set(stats) == {"sd", "mad"}
    assert all(values.shape == (17,) for values in stats.values())
    assert all(np.isfinite(values).all() for values in stats.values())
    assert all((values >= 0).all() for values in stats.values())
