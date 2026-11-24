from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BUILD = load("exp152_build_test", "build_candidate.py")
SCORE = load("exp152_score_test", "score_candidate.py")
DM = load(
    "exp152_dual_moment_test",
    BUILD.ROOT / "experiments/147-dual-moment-emitter/dual_moment.py",
)


def test_preregistration_digest_and_frozen_constants() -> None:
    protocol = BUILD.load_preregistration()
    assert protocol["gene_wise_selector"]["ineligible_lambda"] == 0.5
    assert protocol["gene_wise_selector"]["hct_only_lambda"] == 1.0
    assert protocol["fixed_emission"]["seed_payload"] == (
        "14700:{target}:gene-ramp-v1"
    )
    assert protocol["common_surface"]["common_target_count"] == 120
    assert protocol["common_surface"]["eligible_gene_count_per_target"] == 17_942


def test_exact_midrank_and_consensus_lambda_formula_with_ties() -> None:
    hek = np.array([1.0, 2.0, 2.0, -4.0, 0.0, np.nan, 5.0, -6.0])
    hct = np.array([1.0, -3.0, 2.0, -1.0, 5.0, 1.0, 5.0, 6.0])
    lambdas, diagnostics = BUILD.consensus_gene_lambdas(hek, hct)
    expected = np.array(
        [
            0.5,
            0.5,
            7.0 / 12.0,
            0.5,
            0.5,
            0.5,
            5.0 / 6.0,
            0.5,
        ]
    )
    np.testing.assert_array_equal(lambdas, expected)
    assert diagnostics.eligible_gene_count == 6
    assert diagnostics.agreeing_gene_count == 4
    assert diagnostics.ramped_gene_count == 2


def test_coordinate_interpolation_then_exact_cpm_renormalization() -> None:
    incumbent = np.array([100_000.0, 300_000.0, 600_000.0])
    direction = np.array([300_000.0, 300_000.0, 400_000.0])
    lambdas = np.array([0.5, 0.75, 1.0])
    observed = BUILD.desired_gene_ramp_cpm(incumbent, direction, lambdas)
    raw = incumbent + lambdas * (direction - incumbent)
    expected = BUILD.normalize_exact_cpm(raw)
    np.testing.assert_array_equal(observed, expected)
    assert observed.sum(dtype=np.float64) == 1_000_000.0
    assert np.all(observed >= 0.0)


def test_gene_ramp_exact_margins_pseudobulk_and_determinism() -> None:
    rng = np.random.default_rng(152)
    incumbent = sp.csr_matrix(rng.poisson(4.0, size=(8, 12)).astype(np.int64) + 1)
    direction_dense = incumbent.toarray()
    direction_dense[:, :4] += np.arange(8, dtype=np.int64)[:, None]
    direction = sp.csr_matrix(direction_dense)
    lambdas = np.linspace(0.5, 1.0, incumbent.shape[1], dtype=np.float64)
    seed = BUILD.stable_gene_ramp_seed(14700, "SYNTH")
    first, diagnostics = BUILD.emit_gene_ramp(
        incumbent,
        direction,
        lambdas,
        dual_moment=DM,
        depth_contrast=0.5,
        quasi_spread=0.25,
        seed=seed,
    )
    second, _ = BUILD.emit_gene_ramp(
        incumbent,
        direction,
        lambdas,
        dual_moment=DM,
        depth_contrast=0.5,
        quasi_spread=0.25,
        seed=seed,
    )
    np.testing.assert_array_equal(first.toarray(), second.toarray())
    np.testing.assert_array_equal(
        np.asarray(first.sum(axis=0)).ravel(),
        np.asarray(incumbent.sum(axis=0)).ravel(),
    )
    np.testing.assert_array_equal(
        DM.scorer_pseudobulk(first),
        DM.scorer_pseudobulk(incumbent),
    )
    assert diagnostics.lambda_min == 0.5
    assert diagnostics.lambda_max == 1.0
    assert diagnostics.max_row_total <= 1_000_000


def test_uniform_lambda_one_matches_frozen_dual_moment_adapter() -> None:
    rng = np.random.default_rng(1521)
    incumbent = sp.csr_matrix(rng.poisson(5.0, size=(10, 9)).astype(np.int64) + 1)
    direction = sp.csr_matrix(rng.poisson(6.0, size=(10, 9)).astype(np.int64) + 1)
    seed = 812_152
    gene_ramp, _ = BUILD.emit_gene_ramp(
        incumbent,
        direction,
        np.ones(incumbent.shape[1], dtype=np.float64),
        dual_moment=DM,
        depth_contrast=0.5,
        quasi_spread=0.25,
        seed=seed,
    )
    uniform, _ = DM.emit_dual_moment(
        incumbent,
        direction,
        blend_lambda=1.0,
        depth_contrast=0.5,
        quasi_spread=0.25,
        seed=seed,
    )
    np.testing.assert_array_equal(gene_ramp.toarray(), uniform.toarray())


def test_builder_refuses_without_future_source_authorization(
    tmp_path: Path,
) -> None:
    protocol = BUILD.load_preregistration()
    missing = tmp_path / "SOURCE_BUILD_AUTHORIZATION.v1.json"
    with pytest.raises(PermissionError, match="future one-time"):
        BUILD.verify_build_authorization(protocol, missing)


def test_scorer_refuses_without_future_scoring_authorization(
    tmp_path: Path,
) -> None:
    missing = tmp_path / "SCORING_AUTHORIZATION.v1.json"
    with pytest.raises(PermissionError, match="future one-time"):
        SCORE.verify_scoring_authorization(missing)


def test_frozen_comparator_selection_and_gate_boundaries() -> None:
    assert (
        SCORE.choose_best_frozen_comparator(
            {
                "uniform_lambda_1": 0.1,
                "targetwise_selector": 0.2,
            }
        )
        == "targetwise_selector"
    )
    assert (
        SCORE.choose_best_frozen_comparator(
            {
                "uniform_lambda_1": 0.2,
                "targetwise_selector": 0.2,
            }
        )
        == "uniform_lambda_1"
    )

    protocol = SCORE.load_preregistration()
    a_delta = {
        "de_wilcoxon_lfc_nmae": 0.001,
        "de_wilcoxon_direction_reach_raw": 0.001,
        "de_wilcoxon_direction_fidelity_yield_raw": 0.0,
        "de_wilcoxon_sig_jaccard": -0.015,
        "expr_mse_unbiased_capped_norm": 0.0,
        "avg_score": 0.006,
    }
    versus_incumbent = {
        "scaled_delta": a_delta,
        "all_per_target_pds_exact": True,
        "paired_target_bootstrap": {"lower_2_5pct": 1e-12},
    }
    versus_best = {
        "scaled_delta": {"avg_score": 0.001},
        "paired_target_bootstrap": {"lower_2_5pct": -0.001},
    }
    result = SCORE.evaluate_primary_gates(
        protocol,
        versus_incumbent,
        versus_best,
    )
    assert result["passes"]
    versus_incumbent["paired_target_bootstrap"]["lower_2_5pct"] = 0.0
    assert not SCORE.evaluate_primary_gates(
        protocol,
        versus_incumbent,
        versus_best,
    )["passes"]
