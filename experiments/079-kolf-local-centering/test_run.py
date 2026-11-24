"""Synthetic tests for experiment 079.

No real donor, KOLF, or whitening matrix is opened by these tests.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPEC = importlib.util.spec_from_file_location("exp079_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def load_exp075():
    path = ROOT / "experiments/075-postwhiten-local-centering/run.py"
    spec = importlib.util.spec_from_file_location("exp075_reference", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_independent_transform_matches_frozen_exp075_reference() -> None:
    reference = load_exp075()
    rng = np.random.default_rng(7901)
    query = rng.normal(size=(300, 41))
    query[7] = 0.0
    observed, diagnostics = RUN.reciprocal_local_center(query)
    expected, reference_diagnostics = reference.apply_reciprocal_local_centering(
        query,
        k=64,
        beta=0.5,
        minimum_mutual=2,
    )
    np.testing.assert_allclose(observed, expected, rtol=0.0, atol=1e-12)
    assert diagnostics.fallback_rows == reference_diagnostics["fallback_rows"]
    assert diagnostics.max_row_norm_error < 1e-10


def test_transform_is_permutation_equivariant_and_norm_preserving() -> None:
    rng = np.random.default_rng(7902)
    query = rng.normal(size=(300, 29))
    transformed, _ = RUN.reciprocal_local_center(query)
    permutation = rng.permutation(len(query))
    inverse = np.argsort(permutation)
    permuted, _ = RUN.reciprocal_local_center(query[permutation])
    np.testing.assert_allclose(transformed, permuted[inverse], atol=1e-12)
    np.testing.assert_allclose(
        np.linalg.norm(transformed, axis=1),
        np.linalg.norm(query, axis=1),
        atol=1e-12,
    )


def test_full_300_query_is_frozen_before_truth_loader() -> None:
    official = [f"P{i:03d}" for i in range(300)] + [f"G{i:03d}" for i in range(17)]
    panel = official[:300]
    state: dict[str, object] = {}

    def build_query():
        state["query_built"] = True
        return np.arange(300 * len(official), dtype=np.float64).reshape(
            300, len(official)
        ), {"synthetic": True}

    def load_truth(score_positions: np.ndarray):
        assert state.get("query_built") is True
        assert len(score_positions) == 17
        state["truth_loaded"] = True
        return "truth"

    old_expected = RUN.EXPECTED_OFF_PANEL_GENES
    RUN.EXPECTED_OFF_PANEL_GENES = 17
    try:
        identity, centered, _, audit, truth = RUN.freeze_queries_then_open_truth(
            build_query,
            load_truth,
            official,
            panel,
        )
    finally:
        RUN.EXPECTED_OFF_PANEL_GENES = old_expected
    assert identity.shape == (300, 17)
    assert centered.shape == (300, 17)
    assert audit["transform"]["query_rows_seen_by_transform"] == 300
    assert truth == "truth"
    assert state["truth_loaded"] is True


def test_pds_uses_complete_shared_gallery() -> None:
    truth = np.eye(7, dtype=np.float64)
    perfect = RUN.pds_per_target(truth, truth)
    np.testing.assert_allclose(perfect, 1.0)
    swapped = truth.copy()
    swapped[[0, 1]] = swapped[[1, 0]]
    values = RUN.pds_per_target(swapped, truth)
    assert values[0] < 1.0
    assert values[1] < 1.0
    assert len(values) == 7


def test_target_cluster_bootstrap_is_deterministic() -> None:
    deltas = np.linspace(-0.01, 0.03, 273)
    first = RUN.target_cluster_bootstrap(deltas, replicates=500, seed=79)
    second = RUN.target_cluster_bootstrap(deltas, replicates=500, seed=79)
    assert first == second
    assert first["unit"] == "target_gene"
    assert abs(float(first["mean"]) - float(np.mean(deltas))) < 0.001


def test_predeclared_gate_requires_all_three_conditions() -> None:
    passed = RUN.promotion_gate(
        0.0015,
        {"lower_95": 0.00001},
    )
    assert passed["passed"] is True
    too_small = RUN.promotion_gate(
        0.00149,
        {"lower_95": 0.00001},
    )
    assert too_small["raw_pds_delta_positive"] is True
    assert too_small["raw_pds_delta_at_least_0_0015"] is False
    assert too_small["passed"] is False
    uncertain = RUN.promotion_gate(
        0.002,
        {"lower_95": 0.0},
    )
    assert uncertain["passed"] is False


def test_cli_has_no_kolf_tuning_arguments() -> None:
    parser = RUN.build_parser()
    option_strings = {
        option for action in parser._actions for option in action.option_strings
    }
    forbidden = {
        "--k",
        "--ks",
        "--beta",
        "--betas",
        "--minimum-mutual",
        "--bootstrap-seed",
        "--minimum-delta",
    }
    assert not (option_strings & forbidden)


def test_preregistration_matches_code_constants() -> None:
    path = HERE / "preregistered_protocol.json"
    protocol = json.loads(path.read_text())
    assert protocol["candidate"] == {
        "k": RUN.LOCAL_K,
        "beta": RUN.LOCAL_BETA,
        "minimum_mutual": RUN.MINIMUM_MUTUAL,
    }
    assert protocol["promotion_gate"]["minimum_raw_pds_delta"] == (
        RUN.MINIMUM_RAW_PDS_DELTA
    )
    assert protocol["bootstrap"] == {
        "unit": "target_gene",
        "replicates": RUN.BOOTSTRAP_REPLICATES,
        "seed": RUN.BOOTSTRAP_SEED,
        "interval": "percentile_95",
    }
    validated = RUN.validate_preregistration(path)
    assert validated["sha256"] == RUN.EXPECTED_PREREGISTRATION_SHA256


def test_output_guard_rejects_paths_outside_experiment() -> None:
    with tempfile.TemporaryDirectory() as directory:
        try:
            RUN.ensure_output_directory(Path(directory) / "outside")
        except ValueError as error:
            assert "inside experiment 079" in str(error)
        else:
            raise AssertionError("output path outside experiment 079 was accepted")


def test_aws_scripts_cannot_submit_or_create_compute() -> None:
    text = "\n".join(
        (HERE / "aws" / name).read_text() for name in ("launch.sh", "run.sh")
    )
    forbidden = (
        "vcc submit",
        "run-instances",
        "start-instances",
        "create-fleet",
        "request-spot-instances",
    )
    for marker in forbidden:
        assert marker not in text


def main() -> int:
    tests = [
        test_independent_transform_matches_frozen_exp075_reference,
        test_transform_is_permutation_equivariant_and_norm_preserving,
        test_full_300_query_is_frozen_before_truth_loader,
        test_pds_uses_complete_shared_gallery,
        test_target_cluster_bootstrap_is_deterministic,
        test_predeclared_gate_requires_all_three_conditions,
        test_cli_has_no_kolf_tuning_arguments,
        test_preregistration_matches_code_constants,
        test_output_guard_rejects_paths_outside_experiment,
        test_aws_scripts_cannot_submit_or_create_compute,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
