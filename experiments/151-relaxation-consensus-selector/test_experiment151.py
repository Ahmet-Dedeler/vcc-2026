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


BUILD = load("exp151_build_test", "build_candidate.py")
SCORE = load("exp151_score_test", "score_candidate.py")
DM = load(
    "exp151_dual_moment_test",
    BUILD.ROOT / "experiments/147-dual-moment-emitter/dual_moment.py",
)


def test_preregistration_digest_and_frozen_selector_hashes() -> None:
    protocol = BUILD.load_preregistration()
    plan, _, _ = BUILD.compute_selector(protocol, verify_sources=True)
    selector = protocol["selector"]
    assert len(plan.top_targets) == 60
    assert len(plan.bottom_targets) == 60
    assert len(plan.hct_only_targets) == 12
    assert len(plan.primary_targets_hct_order) == 120
    assert len(plan.hct_targets_order) == 132
    assert set(plan.eligible_counts.values()) == {17_942}
    assert (
        plan.hashes["top_names_sha256_newline_terminated"]
        == selector["top_names_sha256_newline_terminated"]
    )
    assert (
        plan.hashes["bottom_names_sha256_newline_terminated"]
        == selector["bottom_names_sha256_newline_terminated"]
    )
    assert plan.top_targets[-1] == "IFNAR2"
    assert plan.bottom_targets[0] == "FZD2"


def test_lambda_assignment_covers_all_hct_targets() -> None:
    protocol = BUILD.load_preregistration()
    plan, _, hct = BUILD.compute_selector(protocol, verify_sources=False)
    targets = set(np.asarray(hct["targets"]).astype(str))
    assert set(plan.lambda_by_target) == targets
    assert sum(value == 1.0 for value in plan.lambda_by_target.values()) == 72
    assert sum(value == 0.5 for value in plan.lambda_by_target.values()) == 60
    assert all(plan.lambda_by_target[target] == 1.0 for target in plan.hct_only_targets)


def test_dual_moment_lambdas_preserve_columns_and_pseudobulk() -> None:
    rng = np.random.default_rng(151)
    incumbent = sp.csr_matrix(rng.poisson(4.0, size=(8, 12)).astype(np.int64) + 1)
    direction = incumbent.copy().tolil()
    direction[:, 0] = np.asarray(direction[:, 0].toarray()).ravel() + np.arange(8)
    direction = direction.tocsr()
    for blend_lambda in (0.5, 1.0):
        candidate, diagnostics = DM.emit_dual_moment(
            incumbent,
            direction,
            blend_lambda=blend_lambda,
            depth_contrast=0.5,
            quasi_spread=0.25,
            seed=BUILD.stable_seed(14700, "SYNTH", blend_lambda),
        )
        np.testing.assert_array_equal(
            np.asarray(candidate.sum(axis=0)).ravel(),
            np.asarray(incumbent.sum(axis=0)).ravel(),
        )
        np.testing.assert_array_equal(
            DM.scorer_pseudobulk(candidate),
            DM.scorer_pseudobulk(incumbent),
        )
        assert diagnostics.blend_lambda == blend_lambda
        assert diagnostics.max_row_total <= 1_000_000


def test_source_verification_does_not_hash_held_hct_expression(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protocol = BUILD.load_preregistration()
    original = BUILD.sha256
    held = BUILD.ROOT / protocol["source_artifacts"]["held_hct_harness"]["path"]

    def guarded(path: Path) -> str:
        assert path.resolve() != held.resolve()
        return original(path)

    monkeypatch.setattr(BUILD, "sha256", guarded)
    report = BUILD.verify_source_artifacts(protocol)
    assert report["held_hct_harness"]["verification"].startswith("existence_only")


def test_scorer_refuses_without_future_authorization(tmp_path: Path) -> None:
    missing = tmp_path / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
    with pytest.raises(PermissionError, match="future one-time"):
        SCORE.verify_authorization(missing)


def test_builder_refuses_without_source_build_authorization(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        BUILD,
        "BUILD_AUTHORIZATION_PATH",
        tmp_path / "SOURCE_BUILD_AUTHORIZATION.v1.json",
    )
    with pytest.raises(PermissionError, match="source build is not authorized"):
        BUILD.verify_build_authorization()


def test_scorer_rejects_unbound_authorization(tmp_path: Path) -> None:
    authorization = tmp_path / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
    authorization.write_text(
        json.dumps(
            {
                "experiment": "151-relaxation-consensus-selector",
                "status": "authorized_once",
                "scorer": {"sha256": "0" * 64},
                "bindings_sha256": {},
            }
        )
    )
    with pytest.raises(ValueError, match="scorer source"):
        SCORE.verify_authorization(authorization)


def test_scorer_rejects_partial_binding_set(tmp_path: Path) -> None:
    authorization = tmp_path / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
    authorization.write_text(
        json.dumps(
            {
                "experiment": "151-relaxation-consensus-selector",
                "status": "authorized_once",
                "scorer": {"sha256": SCORE.sha256(Path(SCORE.__file__))},
                "bindings_sha256": {
                    "experiments/151-relaxation-consensus-selector/"
                    "PREREGISTRATION.v1.json": SCORE.sha256(
                        SCORE.PREREGISTRATION_PATH
                    )
                },
            }
        )
    )
    with pytest.raises(ValueError, match="misses required bindings"):
        SCORE.verify_authorization(authorization)


def test_preregistered_gate_logic() -> None:
    protocol = SCORE.load_preregistration()
    de_metrics = (
        "de_wilcoxon_lfc_nmae",
        "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_direction_fidelity_yield_raw",
        "de_wilcoxon_sig_jaccard",
    )
    a_delta = {
        metric: 0.001 if index < 2 else -0.005
        for index, metric in enumerate(de_metrics)
    }
    a_delta.update(
        {
            "avg_score": 0.004,
            "expr_mse_unbiased_capped_norm": 0.0,
        }
    )
    b_delta = {metric: -0.005 for metric in de_metrics}
    b_delta.update(
        {
            "avg_score": 0.001,
            "expr_mse_unbiased_capped_norm": 0.0,
        }
    )
    versus_incumbent = {
        "scaled_delta": a_delta,
        "all_per_target_pds_exact": True,
        "paired_target_bootstrap": {"lower_2_5pct": 1e-12},
    }
    versus_uniform = {
        "scaled_delta": b_delta,
        "all_per_target_pds_exact": True,
        "paired_target_bootstrap": {"lower_2_5pct": -0.001},
    }
    result = SCORE.evaluate_primary_gates(
        protocol,
        versus_incumbent,
        versus_uniform,
    )
    assert result["passes"]

    versus_incumbent["paired_target_bootstrap"]["lower_2_5pct"] = 0.0
    assert not SCORE.evaluate_primary_gates(
        protocol,
        versus_incumbent,
        versus_uniform,
    )["passes"]
