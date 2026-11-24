#!/usr/bin/env python3
"""Synthetic checks for experiment 085.

No donor, whitening, PerturbMatch, CD4, or KOLF matrix is opened.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp085_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def development_features() -> dict[str, np.ndarray]:
    base = (np.arange(RUN.EXP064_TARGETS, dtype=np.float64) + 0.5) / (
        RUN.EXP064_TARGETS
    )
    return {
        "orion_hct116": base.copy(),
        "orion_hek293t": base[::-1].copy(),
    }


def test_midrank_percentile_exact_ties() -> None:
    values = np.asarray([10.0, 0.0, 10.0, 5.0])
    observed = RUN.midrank_percentile(values)
    expected = np.asarray([0.75, 0.125, 0.75, 0.375])
    np.testing.assert_allclose(observed, expected, rtol=0.0, atol=0.0)
    assert observed[0] == observed[2]


def test_low_q20_gate_selected_with_two_donor_support() -> None:
    features = development_features()
    gains = {}
    for held in RUN.FOLDS:
        switched = features[held] <= 0.20
        values = np.full(RUN.EXP064_TARGETS, -0.02)
        values[switched] = 0.03
        gains[held] = values
    selection = RUN.select_gate(features, gains)
    selected = selection["selected"]
    assert selected["rule"] == RUN.GateRule("low", 0.20)
    assert min(selected["fold_deltas"].values()) > 0.0


def test_always_d_fallback_wins_when_one_donor_is_negative() -> None:
    features = development_features()
    gains = {
        "orion_hct116": np.ones(RUN.EXP064_TARGETS),
        "orion_hek293t": -np.ones(RUN.EXP064_TARGETS),
    }
    selected = RUN.select_gate(features, gains)["selected"]
    assert selected["identity"] is True
    assert selected["macro_delta"] == 0.0


def test_full_selection_permutation_is_deterministic() -> None:
    features = development_features()
    rng = np.random.default_rng(850)
    paired = rng.normal(size=RUN.EXP064_TARGETS)
    gains = {
        "orion_hct116": paired,
        "orion_hek293t": 0.5 * paired + rng.normal(scale=0.2, size=RUN.EXP064_TARGETS),
    }
    first = RUN.full_selection_permutation_null(
        features, gains, repetitions=100, seed=85
    )
    second = RUN.full_selection_permutation_null(
        features, gains, repetitions=100, seed=85
    )
    assert first == second
    assert 0.0 < first["one_sided_p"] <= 1.0


def test_training_gate_requires_every_condition() -> None:
    selection = {
        "selected": {
            "identity": False,
            "macro_delta": 0.002,
            "fold_deltas": {
                "orion_hct116": 0.001,
                "orion_hek293t": 0.003,
            },
        }
    }
    passed = RUN.training_promotion_gate(selection, {"one_sided_p": 0.04})
    assert passed["passed"] is True
    failed = RUN.training_promotion_gate(selection, {"one_sided_p": 0.06})
    assert failed["passed"] is False


def test_angular_displacement_and_pds_use_full_rows() -> None:
    left = np.eye(5)
    right = left.copy()
    right[[0, 1]] = right[[1, 0]]
    displacement = RUN.angular_displacement(left, right)
    assert displacement[0] == 1.0
    assert displacement[2] == 0.0
    perfect = RUN.pds_per_target(left, left)
    np.testing.assert_allclose(perfect, 1.0)
    swapped = RUN.pds_per_target(right, left)
    assert swapped[0] < 1.0 and swapped[1] < 1.0


def test_nonfinite_pds_fails_closed() -> None:
    prediction = np.eye(3)
    truth = np.eye(3)
    truth[0, 0] = np.nan
    try:
        RUN.pds_per_target(prediction, truth)
    except ValueError as error:
        assert "finite" in str(error)
    else:
        raise AssertionError("non-finite PDS input was accepted")


def test_source_clusters_are_deterministic_and_source_only() -> None:
    query = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [0.99, 0.01, 0.0],
            [0.0, 1.0, 0.0],
            [0.01, 0.99, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.01, 0.99],
        ]
    )
    first, audit_first = RUN.source_effect_clusters(
        query, reciprocal_k=1, cosine_threshold=0.8
    )
    second, audit_second = RUN.source_effect_clusters(
        query, reciprocal_k=1, cosine_threshold=0.8
    )
    np.testing.assert_array_equal(first, second)
    assert audit_first == audit_second
    assert len(np.unique(first)) == 3


def test_cluster_bootstrap_and_routing_null_are_deterministic() -> None:
    deltas = np.linspace(-0.01, 0.03, 24)
    labels = np.repeat(np.arange(8), 3)
    first = RUN.source_cluster_bootstrap(deltas, labels, repetitions=200, seed=85)
    second = RUN.source_cluster_bootstrap(deltas, labels, repetitions=200, seed=85)
    assert first == second
    supported = np.ones(24, dtype=bool)
    assignments = np.zeros(24, dtype=bool)
    assignments[:6] = True
    routing_first = RUN.matched_routing_permutation(
        deltas,
        assignments,
        supported,
        repetitions=200,
        seed=86,
    )
    routing_second = RUN.matched_routing_permutation(
        deltas,
        assignments,
        supported,
        repetitions=200,
        seed=86,
    )
    assert routing_first == routing_second
    assert routing_first["c_routes"] == 6


def test_too_few_source_clusters_rejects_without_crashing() -> None:
    result = RUN.source_cluster_bootstrap(
        np.linspace(-0.01, 0.02, 12),
        np.repeat(np.arange(4), 3),
        repetitions=10,
        minimum_clusters=8,
    )
    assert result["valid"] is False
    assert result["lower_95"] is None


def test_degenerate_routing_null_rejects_without_crashing() -> None:
    gains = np.arange(10, dtype=np.float64)
    supported = np.ones(10, dtype=bool)
    assignments = np.zeros(10, dtype=bool)
    result = RUN.matched_routing_permutation(
        gains, assignments, supported, repetitions=10
    )
    assert result["valid"] is False
    assert result["one_sided_p"] == 1.0


def test_unsupported_targets_are_rejected_by_routing_null() -> None:
    gains = np.arange(10, dtype=np.float64)
    supported = np.asarray([True] * 8 + [False] * 2)
    assignments = np.asarray([True, False] * 5)
    try:
        RUN.matched_routing_permutation(gains, assignments, supported, repetitions=10)
    except ValueError as error:
        assert "unsupported" in str(error)
    else:
        raise AssertionError("unsupported C route was accepted")


def test_production_ranks_only_pm_supported_targets() -> None:
    raw = np.asarray([0.0, 0.1, 0.2, 0.3, 0.4])
    supported = np.asarray([False, True, True, True, False])
    feature, assignments = RUN.production_routes(
        raw, supported, RUN.GateRule("low", 0.5)
    )
    np.testing.assert_array_equal(feature[~supported], [-1.0, -1.0])
    np.testing.assert_allclose(feature[supported], [1 / 6, 3 / 6, 5 / 6])
    np.testing.assert_array_equal(assignments, [False, True, True, False, False])


def test_execute_after_freeze_enforces_boundary_and_pre_gate() -> None:
    events: list[str] = []

    def freeze():
        events.append("freeze")
        return {"sealed_before_truth": True}

    def truth(_):
        assert events == ["freeze"]
        events.append("truth")
        return "held"

    frozen, held = RUN.execute_after_freeze(False, freeze, truth)
    assert frozen is None and held is None and events == []
    frozen, held = RUN.execute_after_freeze(True, freeze, truth)
    assert frozen == {"sealed_before_truth": True}
    assert held == "held"
    assert events == ["freeze", "truth"]


def test_access_marker_is_durable_before_truth_can_continue() -> None:
    events: list[str] = []
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        output = Path(directory)
        durable_freeze = {
            "freeze_manifest_sha256": "a" * 64,
            "uri": "s3://synthetic/output/exp085-kolf-pm-gate/run/pre_kolf/",
        }

        def uploader(path: Path, uri: str, digest: str) -> None:
            assert path.is_file()
            assert (
                uri == "s3://synthetic/output/exp085-kolf-pm-gate/"
                "KOLF_ACCESS_CONSUMED.json"
            )
            assert RUN.sha256_file(path) == digest
            events.append("durable")

        record = RUN.mark_kolf_access(
            output,
            {
                "freeze_manifest_sha256": "a" * 64,
            },
            durable_freeze,
            ("s3://synthetic/output/exp085-kolf-pm-gate/KOLF_ACCESS_CONSUMED.json"),
            uploader=uploader,
        )
        assert events == ["durable"]
        assert record["durable_sha256_verified_before_kolf_open"] is True


def test_pre_kolf_freeze_is_fully_persisted_before_access_marker() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        freeze_dir = Path(directory) / "pre_kolf"
        freeze_dir.mkdir()
        first = freeze_dir / "first.npy"
        second = freeze_dir / "freeze_manifest.json"
        first.write_bytes(b"frozen-array")
        second.write_bytes(b"frozen-manifest")
        frozen = {
            "freeze_manifest_path": str(second),
            "freeze_manifest_sha256": RUN.sha256_file(second),
            "files": {"first.npy": RUN.sha256_file(first)},
        }
        uploaded: list[tuple[str, str]] = []

        def uploader(path: Path, uri: str, digest: str) -> None:
            assert RUN.sha256_file(path) == digest
            uploaded.append((path.name, uri))

        result = RUN.persist_pre_kolf_freeze(
            frozen,
            "s3://synthetic/output/exp085-kolf-pm-gate/run/pre_kolf/",
            uploader=uploader,
        )
        assert [name for name, _ in uploaded] == [
            "first.npy",
            "freeze_manifest.json",
        ]
        assert result["files_verified"] == 2
        assert result["atomic_no_overwrite_uploads"] is True


def test_freeze_artifacts_are_written_and_hashed_before_truth() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        output = Path(directory)
        panel = tuple(f"T{i}" for i in range(RUN.PANEL_TARGETS))
        arm_c = np.arange(
            RUN.PANEL_TARGETS * RUN.OFFICIAL_GENES, dtype=np.float64
        ).reshape(RUN.PANEL_TARGETS, RUN.OFFICIAL_GENES)
        arm_d = arm_c + 1.0
        assignments = np.zeros(RUN.PANEL_TARGETS, dtype=bool)
        assignments[::3] = True
        gated = arm_d.copy()
        gated[assignments] = arm_c[assignments]
        production = {
            "arm_c": arm_c,
            "arm_d": arm_d,
            "gated": gated,
            "raw_feature": np.linspace(0.0, 1.0, RUN.PANEL_TARGETS),
            "feature": np.linspace(0.0, 1.0, RUN.PANEL_TARGETS),
            "assignments_c": assignments,
            "pm_supported": np.ones(RUN.PANEL_TARGETS, dtype=bool),
            "source_counts": np.resize(
                np.asarray([3, 2, 1, 0], dtype=np.uint8), RUN.PANEL_TARGETS
            ),
            "clusters": np.arange(RUN.PANEL_TARGETS) % 12,
            "audit": {"synthetic": True},
        }
        frozen = RUN.freeze_production_artifacts(
            output,
            panel,
            RUN.GateRule("low", 0.5),
            production,
        )
        assert frozen["sealed_before_truth"] is True
        manifest_path = Path(frozen["freeze_manifest_path"])
        assert manifest_path.is_file()
        assert RUN.sha256_file(manifest_path) == frozen["freeze_manifest_sha256"]
        assert (output / "pre_kolf" / "arm_C.npy").is_file()


def test_output_guard_and_one_time_marker_refuse_reuse() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        output = Path(directory) / "output"
        RUN.ensure_output_directory(output)
        (output / "occupied").write_text("x")
        try:
            RUN.ensure_output_directory(output)
        except FileExistsError:
            pass
        else:
            raise AssertionError("nonempty output directory was reused")
    with tempfile.TemporaryDirectory() as outside:
        try:
            RUN.ensure_output_directory(Path(outside) / "output")
        except ValueError:
            pass
        else:
            raise AssertionError("output path outside experiment 085 was accepted")


def test_protocol_and_manifest_match_code_contract() -> None:
    validated = RUN.validate_protocol(HERE / "preregistered_protocol.json")
    assert validated["protocol"]["held_truth"]["candidate_count"] == 1
    manifest = RUN.load_artifact_manifest(HERE / "artifact_manifest.json")
    assert manifest["fixed_axes"]["kolf_score_genes"] == RUN.KOLF_SCORE_GENES
    assert (
        manifest["fixed_axes"]["kolf_score_axis_sha256"]
        == RUN.EXPECTED_KOLF_SCORE_AXIS_SHA256
    )


def test_cli_exposes_no_gate_or_kolf_tuning_arguments() -> None:
    parser = RUN.build_parser()
    options = {option for action in parser._actions for option in action.option_strings}
    forbidden = {
        "--quantiles",
        "--direction",
        "--threshold",
        "--training-permutations",
        "--bootstrap-replicates",
        "--bootstrap-seed",
        "--routing-permutations",
        "--routing-seed",
        "--cluster-k",
        "--cluster-cosine",
        "--minimum-delta",
        "--expected-hct-whitener-sha256",
        "--expected-hek-whitener-sha256",
    }
    assert not (options & forbidden)


def test_kolf_control_is_never_loaded_and_submission_is_absent() -> None:
    text = (HERE / "run.py").read_text()
    assert 'np.load(args.kolf / "b_ctrl.npy"' not in text
    assert "vcc submit" not in text
    assert "cell-eval2" not in text


def test_aws_scripts_cannot_create_compute_or_submit() -> None:
    scripts = [HERE / "aws" / "launch.sh", HERE / "aws" / "run.sh"]
    if not all(path.is_file() for path in scripts):
        return
    text = "\n".join(path.read_text() for path in scripts)
    forbidden = (
        "run-instances",
        "start-instances",
        "create-fleet",
        "request-spot-instances",
        "vcc submit",
        "cell-eval2",
    )
    for marker in forbidden:
        assert marker not in text
    run_text = (HERE / "run.py").read_text()
    assert '"--if-none-match"' in run_text
    assert "EXP085_DURABLE_FREEZE_URI" in run_text
    worker = (HERE / "aws" / "run.sh").read_text()
    launcher = (HERE / "aws" / "launch.sh").read_text()
    assert "output/exp062-cd4-post-whiten/whiteners/" not in worker
    assert '${SOURCE_PREFIX}/held-operators/held_${held}.${suffix}' in worker
    assert "output/exp062-cd4-post-whiten/whiteners/" in launcher
    assert '${SOURCE_PREFIX}/held-operators/held_${held}.${suffix}' in launcher


def test_source_strata_cover_complete_partial_direct_missing() -> None:
    deltas = np.asarray([0.3, 0.2, -0.1, 0.0])
    counts = np.asarray([3, 2, 1, 0])
    strata = RUN.source_strata(deltas, counts)
    assert strata["complete"]["targets"] == 1
    assert strata["partial"]["targets"] == 2
    assert strata["direct"]["targets"] == 3
    assert strata["missing"]["targets"] == 1


def test_complete_evaluation_path_uses_one_fixed_candidate() -> None:
    rng = np.random.default_rng(855)
    old_targets = RUN.KOLF_TARGETS
    old_genes = RUN.KOLF_SCORE_GENES
    RUN.KOLF_TARGETS = 12
    RUN.KOLF_SCORE_GENES = 12
    try:
        truth, _ = np.linalg.qr(rng.normal(size=(12, 12)))
        arm_d_eval = truth.copy()
        arm_d_eval[[0, 1]] = arm_d_eval[[1, 0]]
        arm_d_eval[[2, 3]] = arm_d_eval[[3, 2]]
        arm_c_eval = arm_d_eval.copy()
        assignments = np.zeros(12, dtype=bool)
        assignments[:4] = True
        arm_c_eval[assignments] = truth[assignments]
        arm_c = np.zeros((12, 12), dtype=np.float64)
        arm_d = np.zeros_like(arm_c)
        arm_c[:, :12] = arm_c_eval
        arm_d[:, :12] = arm_d_eval
        gated = arm_d.copy()
        gated[assignments] = arm_c[assignments]
        held = RUN.HeldTruth(
            values=truth,
            targets=tuple(f"T{i}" for i in range(12)),
            score_genes=tuple(f"G{i}" for i in range(12)),
            panel_rows=np.arange(12),
            official_columns=np.arange(12),
            audit={"b_ctrl_loaded_or_used": False},
        )
        production = {
            "arm_c": arm_c,
            "arm_d": arm_d,
            "gated": gated,
            "assignments_c": assignments,
            "pm_supported": np.ones(12, dtype=bool),
            "clusters": np.arange(12),
            "source_counts": np.resize(np.asarray([3, 2, 1, 0], dtype=np.uint8), 12),
        }
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            result = RUN.evaluate_kolf_once(Path(directory), held, production)
            assert result["scores"]["gated_minus_D"] > 0.0
            assert result["routing"]["C_routes_on_supported"] == 4
            assert (Path(directory) / "per_target_pds.csv").is_file()
    finally:
        RUN.KOLF_TARGETS = old_targets
        RUN.KOLF_SCORE_GENES = old_genes


def main() -> int:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} synthetic checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
