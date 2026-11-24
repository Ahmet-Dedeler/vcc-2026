#!/usr/bin/env python3
"""Low-memory synthetic tests for make_submission_pm_cd4.py."""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "make_submission_pm_cd4.py"
ROOT = HERE.parent
VALIDATOR_PATH = ROOT / "experiments/071-production-pm-cd4/validate_authorization.py"
RUNNER_PATH = ROOT / "experiments/071-production-pm-cd4/aws/run.sh"
LAUNCHER_PATH = ROOT / "experiments/071-production-pm-cd4/aws/launch.sh"
README_PATH = ROOT / "experiments/071-production-pm-cd4/README.md"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "test_target_make_submission_pm_cd4",
        MODULE_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


M = load_module()


def load_validator():
    spec = importlib.util.spec_from_file_location(
        "test_target_validate_exp071_authorization",
        VALIDATOR_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


V = load_validator()


class LinearWhitener:
    def __init__(self, matrix: np.ndarray):
        self.matrix = np.asarray(matrix, dtype=np.float64)

    def apply(self, values, gamma, lam, renorm=False):
        assert gamma == M.WHITEN_GAMMA
        assert lam == M.WHITEN_LAM
        assert renorm is False
        return np.asarray(values, dtype=np.float64) @ self.matrix


class ParityModel:
    names = list(M.PRODUCTION_DONORS)
    target_genes = np.asarray(["G0", "G1", "G2", "G3"])
    effect_scale = 6.817
    pert_index = {
        "k562_gwps": {"T0": 0, "T1": 1},
        "orion_hct116": {"T0": 0, "T1": 1},
        "orion_hek293t": {"T0": 0, "T1": 1},
    }

    def effect(self, target, alpha, beta, *, weights):
        assert weights == {name: 1.0 for name in self.names}
        assert beta == M.BETA
        base = {
            "T0": np.asarray([1.0, 2.0, -1.0, 0.5]),
            "T1": np.asarray([-0.5, 1.5, 2.5, -2.0]),
        }[target]
        return alpha * base


class SupportModel:
    names = list(M.PRODUCTION_DONORS)
    target_genes = np.asarray(["G0", "G1", "G2", "G3"])
    effect_scale = 6.0
    _scale = {"k562_gwps": 2.0}
    tcols = {"k562_gwps": np.arange(4, dtype=np.int64)}
    dcols = {"k562_gwps": np.arange(4, dtype=np.int64)}
    pert_index = {
        "k562_gwps": {"both": 0, "k_only": 1},
        "orion_hct116": {"both": 0, "k_only": 1, "pm_only": 2},
        "orion_hek293t": {"both": 0, "pm_only": 1},
    }
    s = {
        "k562_gwps": np.asarray(
            [
                [2.0, 4.0, 6.0, 8.0],
                [1.0, 3.0, 5.0, 7.0],
            ]
        )
    }


class AxisModel:
    names = list(M.PRODUCTION_DONORS)
    target_genes = np.asarray(["PANEL", "A", "B", "C", "D"])
    tcols = {
        "k562_gwps": np.asarray([0, 1, 2], dtype=np.int64),
        "orion_hct116": np.asarray([0, 2, 3], dtype=np.int64),
        "orion_hek293t": np.asarray([0, 3, 4], dtype=np.int64),
    }


def test_disabled_path_matches_v9_effect_semantics() -> None:
    model = ParityModel()
    targets = ("T0", "T1")
    columns = np.asarray([0, 1, 2], dtype=np.int64)
    transform = LinearWhitener(
        np.asarray(
            [
                [1.2, 0.1, 0.0],
                [0.0, 0.8, 0.2],
                [0.1, 0.0, 1.1],
            ]
        )
    )

    baseline = np.vstack(
        [
            model.effect(
                target,
                M.ALPHA,
                M.BETA,
                weights={name: 1.0 for name in model.names},
            )
            for target in targets
        ]
    ).astype(np.float64)
    expected = baseline.copy()
    before_operator = float(np.median(np.linalg.norm(expected[:, columns], axis=1)))
    transformed = transform.apply(
        expected[:, columns],
        M.WHITEN_GAMMA,
        M.WHITEN_LAM,
        renorm=False,
    )
    after_operator = float(np.median(np.linalg.norm(transformed, axis=1)))
    expected[:, columns] = transformed * (before_operator / after_operator)
    actual, audit = M.build_production_effects(
        model,
        targets,
        transform,
        columns,
        np.asarray([3], dtype=np.int64),
        pm=None,
        cd4=None,
    )
    np.testing.assert_array_equal(actual, expected)
    assert audit["perturbmatch"] == {"enabled": False}
    assert audit["cd4"] == {"enabled": False}
    assert np.isclose(
        audit["whitening"]["operator_preserve_scale"],
        before_operator / after_operator,
    )
    assert audit["baseline_pre_sha256"] == audit["candidate_pre_sha256"]
    assert audit["post_whiten_sha256"] == audit["final_sha256"]


def test_perturbmatch_target_support_policy_and_source_share() -> None:
    model = SupportModel()
    targets = ("both", "k_only", "pm_only", "neither")
    baseline = np.arange(16, dtype=np.float64).reshape(4, 4)
    pm = {
        "scale": 2.0,
        "target_columns": np.asarray([0, 2], dtype=np.int64),
        "source_columns": np.asarray([0, 1], dtype=np.int64),
        "mean": np.asarray([1.0, 2.0]),
        "delta": np.asarray(
            [
                [5.0, 10.0],
                [99.0, 99.0],
            ]
        ),
        "pert_index": {"both": 0, "pm_only": 1},
    }
    actual, audit = M.apply_perturbmatch_k562_replacement(
        baseline,
        model,
        pm,
        targets,
        alpha=1.0,
    )

    # "both" has three incumbent source votes, so the K562 slot share is 1/3.
    incumbent_standardized = np.asarray([1.0, 2.0, 3.0, 4.0])
    pm_standardized = np.asarray([(5.0 - 1.0) / 2.0, (10.0 - 2.0) / 2.0])
    expected = baseline[0].copy()
    expected[[0, 2]] += (
        (1.0 / 3.0)
        * model.effect_scale
        * M.PM_BLEND
        * (pm_standardized - incumbent_standardized[[0, 2]])
    )
    np.testing.assert_allclose(actual[0], expected)
    np.testing.assert_array_equal(actual[0, [1, 3]], baseline[0, [1, 3]])
    np.testing.assert_array_equal(actual[1:], baseline[1:])
    assert audit["both_supported_modified_count"] == 1
    assert audit["k562_supported_pm_missing_preserved_count"] == 1
    assert audit["pm_supported_k562_missing_preserved_count"] == 1
    assert audit["neither_supported_preserved_count"] == 1
    assert audit["modified_target_source_vote_histogram"] == {"3": 1}


def test_panel_exclusion_for_operator_and_cd4_residual() -> None:
    operator_columns = M.production_operator_columns(AxisModel(), ("PANEL",))
    np.testing.assert_array_equal(operator_columns, np.asarray([1, 2, 3, 4]))

    baseline = np.asarray(
        [
            [4.0, 3.0, 2.0, 1.0],
            [1.0, 2.0, 3.0, 4.0],
        ]
    )
    cd4 = M.EXP062.CD4Prepared(
        mode=M.CD4_MODE,
        mean_effect=np.asarray(
            [
                [0.0, 2.0, 0.0, 1.0],
                [0.0, 0.0, 3.0, 0.0],
            ]
        ),
        gene_votes=np.asarray(
            [
                [0, 3, 0, 2],
                [0, 0, 3, 0],
            ],
            dtype=np.uint8,
        ),
        target_condition_count=np.asarray([3, 3], dtype=np.uint8),
        condition_details=[],
    )
    output, audit = M.apply_frozen_cd4_residual(
        baseline,
        cd4,
        np.asarray([1, 2, 3], dtype=np.int64),
        np.asarray([0], dtype=np.int64),
    )
    scale = audit["post_add_preserve"]["post_add_global_preserve_scale"]
    np.testing.assert_allclose(output[:, 0], baseline[:, 0] * scale, rtol=0, atol=0)
    assert audit["panel_residual_nonzero_entries"] == 0

    invalid = M.EXP062.CD4Prepared(
        mode=M.CD4_MODE,
        mean_effect=np.asarray([[1.0, 0.0, 0.0, 0.0]]),
        gene_votes=np.asarray([[1, 0, 0, 0]], dtype=np.uint8),
        target_condition_count=np.asarray([3], dtype=np.uint8),
        condition_details=[],
    )
    try:
        M.apply_frozen_cd4_residual(
            baseline[:1],
            invalid,
            np.asarray([1, 2, 3], dtype=np.int64),
            np.asarray([0], dtype=np.int64),
        )
    except AssertionError as error:
        assert "panel coordinate" in str(error)
    else:
        raise AssertionError("nonzero panel-coordinate CD4 residual was accepted")


def test_complete_mean3_support_semantics() -> None:
    prepared = M.EXP062.CD4Prepared(
        mode=M.CD4_MODE,
        mean_effect=np.asarray(
            [
                [1.0, 2.0],
                [3.0, 4.0],
                [5.0, 6.0],
                [7.0, 8.0],
            ]
        ),
        gene_votes=np.asarray(
            [
                [3, 3],
                [2, 2],
                [0, 0],
                [1, 3],
            ],
            dtype=np.uint8,
        ),
        target_condition_count=np.asarray([3, 2, 0, 3], dtype=np.uint8),
        condition_details=[],
    )
    frozen = M.freeze_complete_mean3(prepared)
    np.testing.assert_array_equal(frozen.mean_effect[0], [1.0, 2.0])
    np.testing.assert_array_equal(frozen.mean_effect[1], [0.0, 0.0])
    np.testing.assert_array_equal(frozen.mean_effect[2], [0.0, 0.0])
    np.testing.assert_array_equal(frozen.mean_effect[3], [7.0, 8.0])
    np.testing.assert_array_equal(frozen.gene_votes[1], [0, 0])
    np.testing.assert_array_equal(frozen.gene_votes[2], [0, 0])

    audit = M.cd4_support_audit(frozen)
    assert audit["fully_supported_targets"] == 2
    assert audit["incomplete_targets_preserved_before_global_scale"] == 2
    assert audit["target_condition_count_histogram"] == {
        "0": 1,
        "2": 1,
        "3": 2,
    }


def test_rng_streams_are_deterministically_separated() -> None:
    pick_a, emit_a = M.make_rngs(17)
    first_a = pick_a.choice(1000, 50, replace=False)
    emit_a.random(10)
    second_a = pick_a.choice(1000, 50, replace=False)

    pick_b, emit_b = M.make_rngs(17)
    first_b = pick_b.choice(1000, 50, replace=False)
    emit_b.random(10_000)
    second_b = pick_b.choice(1000, 50, replace=False)

    np.testing.assert_array_equal(first_a, first_b)
    np.testing.assert_array_equal(second_a, second_b)

    # The emission stream is deterministic too, but distinct from the pick stream.
    _, emit_c = M.make_rngs(17)
    np.testing.assert_array_equal(
        M.make_rngs(17)[1].integers(0, 2**31, 20),
        emit_c.integers(0, 2**31, 20),
    )
    assert not np.array_equal(
        M.make_rngs(17)[0].integers(0, 2**31, 20),
        M.make_rngs(17)[1].integers(0, 2**31, 20),
    )


def test_profile_matches_frozen_v9_gate() -> None:
    b_ctrl = np.log1p(np.asarray([0.0, 0.4, 8.0, 80.0]) * M.TS)
    effect = np.asarray([5.0, 4.0, -0.5, 0.25])
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl),
        M.EXP061.profile_from_effect(effect, b_ctrl),
    )


def test_effect_scale_is_off_by_default_and_shrinks_toward_the_control() -> None:
    """`--effect-scale` must be a true no-op at 1.0 and reach the control at 0.0.

    The default has to be bit-identical, not merely close: every build before 096 was
    produced without this parameter, and the v11 byte-reproduction test compares against
    one of them. And at 0.0 the emitted profile must be exactly the context's control
    profile -- that is the board's measured expression-error floor (~1.042 raw, team
    "Kun Lab"'s zero-effect submission), so it is the anchor the knob interpolates from.
    """
    b_ctrl = np.log1p(np.asarray([0.0, 0.4, 8.0, 80.0]) * M.TS)
    effect = np.asarray([5.0, 4.0, -0.5, 0.25])

    # 1.0 is bit-identical to omitting the argument entirely.
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, effect_scale=1.0),
        M.profile_from_effect(effect, b_ctrl),
    )

    # 0.0 reproduces the control profile the gate would emit for a zero effect.
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, effect_scale=0.0),
        M.profile_from_effect(np.zeros_like(effect), b_ctrl),
    )

    # Intermediate values interpolate. NOTE: this is asserted on the DISTANCE to the
    # control, not per gene. `profile_from_effect` renormalises to a composition, so
    # shrinking one gene's effect changes the total and shifts every other gene --
    # a low-effect gene can move away from its full-effect value while the profile as a
    # whole moves toward the control. Asserting per-gene monotonicity fails here, and
    # that failure is the compositional constraint, not a bug.
    none = M.profile_from_effect(effect, b_ctrl, effect_scale=0.0)
    d = [float(np.linalg.norm(M.profile_from_effect(effect, b_ctrl, effect_scale=a) - none))
         for a in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert d[0] == 0.0, "effect_scale=0 must sit exactly on the control profile"
    assert all(x < y for x, y in zip(d, d[1:])), f"distance to control not monotone: {d}"


def test_emission_spread_is_the_shipped_scalar_unless_a_tail_is_asked_for() -> None:
    """An unflagged build must stay bit-identical to v11.

    The emission spread is the only thing experiment 093 added to the count path, so the
    guard that matters is that the DEFAULT vector is exactly what the shipped call built:
    `np.full(profile.size, SPREAD)`, for every argument combination that leaves
    `--tail-spread` unset.
    """
    b_ctrl = np.log1p(np.asarray([0.0, 0.4, 8.0, 80.0, 300.0]) * M.TS)
    effect = np.asarray([5.0, 4.0, -0.5, 0.25, 2.0])
    expected = np.full(effect.size, M.SPREAD)
    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl), expected
    )
    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl, keep_top=3), expected
    )
    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl, keep_top_frac=0.5), expected
    )

    parsed = M.build_parser().parse_args(["--output", "unused.h5ad"])
    assert parsed.tail_spread is None
    assert parsed.keep_top == 0
    assert parsed.keep_top_frac is None


def test_emission_spread_tail_policy_keeps_the_top_k_gated_effects() -> None:
    # Gene 0 is unexpressed and gene 1 sits under the 5 CPM gate, so neither can be
    # called however large its predicted effect is; both belong in the tail. The three
    # gated genes rank 4 (|2.0|) > 2 (|-0.5|) > 3 (|0.25|).
    b_ctrl = np.log1p(np.asarray([0.0, 1e-6, 8.0, 80.0, 300.0]) * M.TS)
    effect = np.asarray([9.0, 8.0, -0.5, 0.25, 2.0])
    assert (1e6 * np.expm1(b_ctrl) / np.expm1(b_ctrl).sum() > M.MIN_CPM).tolist() == [
        False, False, True, True, True
    ]

    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl, tail_spread=1.0),
        np.ones(effect.size),
    )
    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl, tail_spread=0.8, keep_top=2),
        np.asarray([0.8, 0.8, M.SPREAD, 0.8, M.SPREAD]),
    )
    # 3 gated genes * 0.5 -> 2 kept, the same pair.
    np.testing.assert_array_equal(
        M.emission_spread_vector(effect, b_ctrl, tail_spread=0.6, keep_top_frac=0.5),
        np.asarray([0.6, 0.6, M.SPREAD, 0.6, M.SPREAD]),
    )
    # keep_top past the gate cannot promote an ungated gene above a gated one.
    kept = M.emission_spread_vector(effect, b_ctrl, tail_spread=1.0, keep_top=4)
    np.testing.assert_array_equal(kept, np.asarray([1.0, 1.0, M.SPREAD, M.SPREAD,
                                                    M.SPREAD]))

    for bad in ({"tail_spread": 1.5}, {"tail_spread": -0.1},
                {"tail_spread": 0.8, "keep_top": 2, "keep_top_frac": 0.1},
                {"tail_spread": 0.8, "keep_top": -1},
                {"tail_spread": 0.8, "keep_top_frac": 1.5}):
        try:
            M.emission_spread_vector(effect, b_ctrl, **bad)
        except ValueError:
            continue
        raise AssertionError(f"emission_spread_vector accepted {bad}")


def test_atomic_h5_output_publishes_or_cleans_up() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        output = root / "candidate.h5ad"
        with M.atomic_h5_output(output) as handle:
            handle.create_dataset("ok", data=np.asarray([1, 2, 3]))
        assert output.is_file()
        assert not list(root.glob("*.partial-*"))

        failed = root / "failed.h5ad"
        try:
            with M.atomic_h5_output(failed) as handle:
                handle.create_dataset("partial", data=np.asarray([1]))
                raise RuntimeError("intentional")
        except RuntimeError as error:
            assert str(error) == "intentional"
        else:
            raise AssertionError("intentional atomic-write failure did not propagate")
        assert not failed.exists()
        assert not list(root.glob("*.partial-*"))


def test_frozen_constants_match_reused_helpers() -> None:
    M.assert_frozen_contract()
    # Frozen unless deliberately overridden for a board probe; the default is unchanged.
    assert M.ALPHA == float(os.environ.get("VCC_ALPHA", "1.0796"))
    assert M.CD4_FRACTION == float(os.environ.get("VCC_CD4_FRACTION", "0.05"))
    assert M.PM_BLEND == 0.50
    assert M.CD4_MODE == "gene_norm"
    assert M.WHITEN_GAMMA == 0.65
    assert M.SPREAD == 0.25
    assert M.SEED == 0
    assert M.EXPECTED_PM_REPLACEMENT_TARGETS == 248
    assert M.EXPECTED_CD4_COMPLETE_TARGETS == 276


def write_agg(path: Path, values: dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["statistic", *V.METRICS])
        writer.writeheader()
        writer.writerow({"statistic": "mean", **values})


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def make_results(
    baseline_path: Path,
    candidate_path: Path,
    baseline: dict[str, float],
    candidate: dict[str, float],
) -> dict[str, object]:
    metrics = {}
    baseline_scaled = []
    candidate_scaled = []
    for metric in V.METRICS:
        base_scaled = V.scaled(metric, baseline[metric])
        cand_scaled = V.scaled(metric, candidate[metric])
        raw_delta = candidate[metric] - baseline[metric]
        metrics[metric] = {
            "preferred": "lower" if metric in V.LOWER_IS_BETTER else "higher",
            "baseline_raw": baseline[metric],
            "candidate_raw": candidate[metric],
            "raw_candidate_minus_baseline": raw_delta,
            "raw_favorable_delta": (
                -raw_delta if metric in V.LOWER_IS_BETTER else raw_delta
            ),
            "baseline_scaled": base_scaled,
            "candidate_scaled": cand_scaled,
            "scaled_candidate_minus_baseline": cand_scaled - base_scaled,
        }
        baseline_scaled.append(base_scaled)
        candidate_scaled.append(cand_scaled)
    base_avg = sum(baseline_scaled) / len(V.METRICS)
    cand_avg = sum(candidate_scaled) / len(V.METRICS)
    return {
        "experiment": V.EXPERIMENT,
        "status": "scored",
        "comparison": (
            "joint frozen candidate minus exact clean gamma-0.65 HCT baseline"
        ),
        "metrics": metrics,
        "avg6": {
            "baseline": base_avg,
            "candidate": cand_avg,
            "candidate_minus_baseline": cand_avg - base_avg,
        },
        "score_files": {
            "baseline": {
                "path": "/work/exp065/output/baseline-agg_results.csv",
                "sha256": V.file_sha256(baseline_path),
            },
            "candidate": {
                "path": "/work/exp065/output/candidate-run/agg_results.csv",
                "sha256": V.file_sha256(candidate_path),
            },
        },
    }


def refresh_internal_manifest(root: Path, candidate_digest: str) -> None:
    entries = {}
    for relative in V.REQUIRED_INTERNAL_MANIFEST_ENTRIES:
        if relative == "candidate.h5ad":
            entries[relative] = candidate_digest
        else:
            entries[relative] = V.file_sha256(root / relative)
    lines = [
        f"{digest}  /work/exp065/output/{relative}"
        for relative, digest in sorted(entries.items())
    ]
    (root / "artifact-sha256.txt").write_text("\n".join(lines) + "\n")


def refresh_terminal_manifest(root: Path, run_id: str) -> Path:
    prefix = f"output/{V.S3_EXPERIMENT}/{run_id}"
    artifacts = []
    for relative in V.REQUIRED_TERMINAL_ARTIFACTS:
        path = root / relative
        artifacts.append(
            {
                "path": relative,
                "key": f"{prefix}/{relative}",
                "size": path.stat().st_size,
                "sha256": V.file_sha256(path),
                "etag": '"synthetic"',
                "last_modified": "2026-08-23T00:00:00+00:00",
            }
        )
    path = root / "terminal-s3-manifest.json"
    write_json(
        path,
        {
            "schema": 1,
            "experiment": V.EXPERIMENT,
            "run_id": run_id,
            "s3": {"bucket": "synthetic-test", "prefix": prefix},
            "artifacts": artifacts,
        },
    )
    return path


def make_authorization_fixture(root: Path) -> tuple[str, str, Path]:
    run_id = "20260823T040540Z"
    prefix = f"output/{V.S3_EXPERIMENT}/{run_id}"
    baseline = {metric: V.ANCHORS[metric][0] for metric in V.METRICS}
    candidate = {
        metric: baseline[metric] + 0.20 * (V.ANCHORS[metric][1] - baseline[metric])
        for metric in V.METRICS
    }
    baseline_path = root / "baseline-agg_results.csv"
    candidate_path = root / "candidate-run/agg_results.csv"
    write_agg(baseline_path, baseline)
    write_agg(candidate_path, candidate)
    write_json(
        root / "results.json",
        make_results(baseline_path, candidate_path, baseline, candidate),
    )

    candidate_digest = "a" * 64
    (root / "candidate.h5ad.sha256").write_text(
        f"{candidate_digest}  /work/exp065/output/candidate.h5ad\n"
    )
    write_json(
        root / "build-report.json",
        {
            "experiment": V.EXPERIMENT,
            "status": "candidate_built",
            "protocol": dict(V.EXPECTED_PROTOCOL),
            "inputs_sha256": dict(V.EXPECTED_BUILD_INPUT_HASHES),
            "cd4_inputs_sha256": V.EXPECTED_CD4_HASHES,
            "effect_audit": {
                **V.EXPECTED_EFFECT_HASHES,
                "panel_cd4_residual_nonzero_entries": 0,
                "missing_pm_rows_equal_to_baseline_before_cd4": True,
            },
            "axis": {
                "hct_targets": 132,
                "official_genes": 18533,
                "operator_genes": 17942,
            },
            "cd4_audit": {
                "mode": "gene_norm",
                "fully_supported_targets": 125,
            },
            "candidate_h5ad": {
                "path": "/work/exp065/output/candidate.h5ad",
                "bytes": 123,
                "sha256": candidate_digest,
            },
        },
    )
    write_json(
        root / "score-memory-guard.json",
        {
            "max_rss_gib": 55.0,
            "peak_process_tree_rss_gib": 20.0,
            "return_code": 0,
            "killed_for_memory": False,
            "command": ["cell-eval2", "run", "--preset", "vcc2026"],
        },
    )
    write_json(root / "runtime.json", {"python": "3.12.11", "cell_eval2": "0.16.0"})
    (root / "exit-code.txt").write_text("0\n")
    refresh_internal_manifest(root, candidate_digest)
    manifest = refresh_terminal_manifest(root, run_id)
    return run_id, prefix, manifest


def assert_authorization_rejected(
    root: Path,
    manifest: Path,
    run_id: str,
    prefix: str,
    expected_text: str,
) -> None:
    try:
        V.validate_authorization(
            root,
            manifest,
            exp065_run_id=run_id,
            expected_prefix=prefix,
        )
    except (ValueError, FileNotFoundError) as error:
        assert expected_text in str(error), str(error)
    else:
        raise AssertionError(f"authorization unexpectedly accepted: {expected_text}")


def test_authorization_accepts_only_exploratory_scope() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        result = V.validate_authorization(
            root,
            manifest,
            exp065_run_id=run_id,
            expected_prefix=prefix,
        )
        assert result["status"] == "exploratory_packaging_gate_passed"
        assert result["authorization_scope"] == "exploratory_packaging_only"
        assert result["challenge_submission_authorized"] is False
        assert result["semantic_limitations"]["exact_v9_equivalence"] is False


def test_authorization_rejects_run_prefix_and_gate_spoofing() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        assert_authorization_rejected(
            root,
            manifest,
            "20260823T040541Z",
            prefix,
            "expected prefix",
        )
        assert_authorization_rejected(
            root,
            manifest,
            run_id,
            f"{prefix}-spoofed",
            "expected prefix",
        )

    mutations = (
        ("experiment", "spoofed", "build-report experiment"),
        ("status", "completed", "build-report status"),
    )
    for key, value, expected_text in mutations:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            root = Path(directory)
            run_id, prefix, manifest = make_authorization_fixture(root)
            report_path = root / "build-report.json"
            report = json.loads(report_path.read_text())
            report[key] = value
            write_json(report_path, report)
            refresh_internal_manifest(root, "a" * 64)
            manifest = refresh_terminal_manifest(root, run_id)
            assert_authorization_rejected(root, manifest, run_id, prefix, expected_text)

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        report_path = root / "build-report.json"
        report = json.loads(report_path.read_text())
        report["protocol"]["alpha"] = 1.0796
        write_json(report_path, report)
        refresh_internal_manifest(root, "a" * 64)
        manifest = refresh_terminal_manifest(root, run_id)
        assert_authorization_rejected(root, manifest, run_id, prefix, "protocol.alpha")


def test_authorization_rejects_nonfinite_and_hash_relationship_drift() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        result_path = root / "results.json"
        result = json.loads(result_path.read_text())
        result["metrics"]["pds_cosine"]["candidate_raw"] = float("nan")
        write_json(result_path, result)
        refresh_internal_manifest(root, "a" * 64)
        manifest = refresh_terminal_manifest(root, run_id)
        assert_authorization_rejected(
            root, manifest, run_id, prefix, "non-finite numeric value"
        )

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        report_path = root / "build-report.json"
        report = json.loads(report_path.read_text())
        report["candidate_h5ad"]["sha256"] = "b" * 64
        write_json(report_path, report)
        refresh_internal_manifest(root, "a" * 64)
        manifest = refresh_terminal_manifest(root, run_id)
        assert_authorization_rejected(
            root, manifest, run_id, prefix, "candidate build hash"
        )


def test_authorization_rejects_unsafe_or_duplicate_manifest_paths() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        artifact_manifest = root / "artifact-sha256.txt"
        artifact_manifest.write_text(
            artifact_manifest.read_text()
            + f"{'0' * 64}  /work/exp065/output/../escape\n"
        )
        manifest = refresh_terminal_manifest(root, run_id)
        assert_authorization_rejected(
            root, manifest, run_id, prefix, "unsafe exp065 output path"
        )

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        run_id, prefix, manifest = make_authorization_fixture(root)
        terminal = json.loads(manifest.read_text())
        terminal["artifacts"].append(dict(terminal["artifacts"][0]))
        write_json(manifest, terminal)
        assert_authorization_rejected(
            root, manifest, run_id, prefix, "duplicate terminal artifact path"
        )


def test_frozen_input_and_helper_hash_drift_fails_before_use() -> None:
    assert "donors/orion_hct116/delta.npy" in M.FROZEN_INPUT_SHA256
    assert "perturbmatch/delta.npy" in M.FROZEN_INPUT_SHA256
    assert "cd4/cd4_stim48hr/meta.json" in M.FROZEN_INPUT_SHA256
    assert set(M.verify_frozen_implementation_hashes()) == set(
        M.FROZEN_IMPLEMENTATION_SHA256
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "input.bin"
        path.write_bytes(b"frozen")
        expected = {"input": M.file_sha256(path)}
        assert (
            M.verify_expected_hashes({"input": path}, expected, label="synthetic")
            == expected
        )
        path.write_bytes(b"drifted")
        try:
            M.verify_expected_hashes({"input": path}, expected, label="synthetic")
        except ValueError as error:
            assert "hash drift" in str(error)
        else:
            raise AssertionError("drifted frozen input was accepted")


def test_partial_vcc_failure_never_publishes_final_package() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        fake = root / "fake-prep.sh"
        fake.write_text(
            "#!/usr/bin/env bash\n"
            "set -eu\n"
            "out=''\n"
            "while (($#)); do\n"
            '  if [[ "$1" == \'-o\' ]]; then out="$2"; shift 2; else shift; fi\n'
            "done\n"
            'printf partial > "$out"\n'
            "exit 17\n"
        )
        fake.chmod(0o755)
        final = root / "candidate.vcc"
        command = (
            'EXP071_LIBRARY_ONLY=1 source "$1"; '
            'set +e; RUN_ID=synthetic atomic_vcc_prep "$2" "$3"; '
            "status=$?; set -e; "
            'test "$status" = 17; '
            'test ! -e "$2"; '
            'test -z "$(find "$4" -name \'*.partial-*\' -print -quit)"'
        )
        subprocess.run(
            [
                "bash",
                "-c",
                command,
                "bash",
                str(RUNNER_PATH),
                str(final),
                str(fake),
                str(root),
            ],
            check=True,
        )


def test_cloud_contract_has_dual_locks_pinned_sources_and_success_only_vcc() -> None:
    runner = RUNNER_PATH.read_text()
    launcher = LAUNCHER_PATH.read_text()
    readme = README_PATH.read_text()

    assert "exec 8>/work/exp061.lock" in runner
    assert "exec 9>/work/vcc-heavy.lock" in runner
    assert "SOURCE_ARCHIVE_SHA256" in runner
    assert "RUNNER_SHA256" in runner
    assert "SOURCE_MANIFEST_SHA256" in runner
    assert 'PREP_MAX_RSS_GIB="${PREP_MAX_RSS_GIB:-110}"' in runner
    assert "mem_kib < 120000000" in runner
    assert '--max-rss-gib "$PREP_MAX_RSS_GIB"' in runner
    dry_run_block = runner[
        runner.index('rm -f "$DRY_RUN_VCC"') : runner.index(
            'atomic_vcc_prep "$VCC"'
        )
    ]
    assert '--max-rss-gib "$PREP_MAX_RSS_GIB"' in dry_run_block
    build_block = runner[
        runner.index('H5AD="$OUT/v10-pm-cd4.h5ad"') : runner.index(
            'sha256sum "$H5AD"'
        )
    ]
    assert '--max-rss-gib "$MAX_RSS_GIB"' in build_block
    assert "validate_authorization.py" in runner
    assert "atomic_vcc_prep" in runner
    assert ".partial-${run_token}.vcc" in runner
    assert "--exclude '*.vcc'" in runner
    assert "--exclude '*.vcc.sha256'" in runner
    assert 'if [[ "$status" == "0" && -s "$SUCCESS_MARKER" ]]' in runner
    assert "vcc submit" not in runner

    assert 'EXP065_RUN_ID="${EXP065_RUN_ID:?' in launcher
    assert "terminal-s3-manifest.json" in launcher
    assert "source-manifest.json" in launcher
    assert "head-object" in launcher
    assert "sha256sum --check -" in launcher
    assert "challenge_submission_authorized" in launcher

    assert "Experiment 065 completed successfully" in readme
    assert "20260823T040540Z" in readme
    assert "not mathematically identical" in readme
    assert "exploratory packaging" in readme
    assert "does not authorize challenge submission" in readme
    assert "one-day lifecycle" in readme
    assert "not durable" in readme


def test_keep_spread_lowers_the_kept_genes_instead_of_raising_the_tail() -> None:
    """Experiment 094: `keep_spread` decouples the two ends of the spread channel.

    v14 and v15 both raised the tail to 1.0 and left the kept genes at M.SPREAD, which
    only ever shrinks `n_pred`. Both lost `fid` on the board to the coverage term
    min(1, n_pred/n_real). Holding the tail at the shipped 0.25 and driving the kept
    genes to 0.0 instead is an ordering-only intervention: the tail is bit-identical to
    v11, so the call set barely moves, while the genes carrying real predicted effect
    gain rank in the confidence order that `reach` walks.
    """
    effect = np.array([0.0, 0.9, -0.5, 0.2, 0.0])
    b_ctrl = np.log1p(np.array([5e4, 4e4, 3e4, 2e4, 1e4]))

    v11 = M.emission_spread_vector(effect, b_ctrl, min_cpm=0.0)
    assert np.array_equal(v11, np.full(5, M.SPREAD))

    got = M.emission_spread_vector(
        effect, b_ctrl, min_cpm=0.0, tail_spread=M.SPREAD, keep_top=2, keep_spread=0.0
    )
    # top two by |effect| are indices 1 (0.9) and 2 (-0.5)
    assert np.array_equal(got, np.array([M.SPREAD, 0.0, 0.0, M.SPREAD, M.SPREAD]))
    # the tail is bit-identical to v11 -- the ONLY difference from v11 is the kept set
    tail = np.array([0, 3, 4])
    assert np.array_equal(got[tail], v11[tail])

    # keep_spread defaults to M.SPREAD, reproducing the v14/v15 shape exactly
    v15 = M.emission_spread_vector(
        effect, b_ctrl, min_cpm=0.0, tail_spread=1.0, keep_top=2
    )
    assert np.array_equal(v15, np.array([1.0, M.SPREAD, M.SPREAD, 1.0, 1.0]))


def test_keep_spread_is_validated_and_fails_closed() -> None:
    effect = np.array([0.0, 0.9, -0.5])
    b_ctrl = np.log1p(np.array([5e4, 4e4, 3e4]))
    for bad in (-0.1, 1.5, float("nan")):
        try:
            M.emission_spread_vector(
                effect, b_ctrl, min_cpm=0.0, tail_spread=0.25, keep_top=1,
                keep_spread=bad,
            )
        except ValueError:
            pass
        else:
            raise AssertionError(f"--keep-spread {bad} was accepted")
    # a keep_spread with nothing to apply it to is an error, not a silent no-op
    try:
        M.emission_spread_vector(
            effect, b_ctrl, min_cpm=0.0, tail_spread=0.25, keep_spread=0.0
        )
    except ValueError:
        pass
    else:
        raise AssertionError("--keep-spread without --keep-top was accepted")


def test_whitener_registry_selects_a_pair_and_still_fails_closed() -> None:
    """Experiment 091 needs a second operator; the guard must not be weakened to allow it.

    `FROZEN_INPUT_SHA256` pins exactly one whitener pair, which is correct for a single
    production artefact and blocks any A/B against it. The registry swaps the pair by
    NAME, so an operator must be both registered and hash-matched; an unknown id, or a
    registered id whose file drifted by one byte, still aborts before any model input is
    read.
    """
    # The table is now assembled from two halves -- the panel half in
    # scripts/frozen_inputs.json, the rest in source -- so "the default id reproduces the
    # frozen table" means it reproduces their union, with nothing added or dropped.
    base = M.whitener_expected_hashes("union")
    frozen = M.load_frozen_bundle()
    assert base == {**dict(frozen["files"]), **M.FROZEN_INPUT_SHA256}, (
        "the default id must reproduce the frozen table"
    )
    assert not any(key.startswith("bundle/") for key in M.FROZEN_INPUT_SHA256), (
        "the bundle half belongs in scripts/frozen_inputs.json, not in source"
    )

    pruned = M.whitener_expected_hashes("union_pruned")
    npz, js = M.WHITENER_REGISTRY["union_pruned"]
    assert pruned["whitener/whitener_union.npz"] == npz
    assert pruned["whitener/whitener_union.json"] == js
    assert npz != base["whitener/whitener_union.npz"]
    assert js != base["whitener/whitener_union.json"]

    # every other frozen input is untouched by the swap
    for key in base:
        if not key.startswith("whitener/"):
            assert pruned[key] == base[key], f"{key} changed with the operator"

    for bad in ("nope", "", "UNION", "union_pruned "):
        try:
            M.whitener_expected_hashes(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"--whitener-id {bad!r} was accepted")

    # both registered operators must carry two distinct 64-hex digests
    for name, pair in M.WHITENER_REGISTRY.items():
        assert len(pair) == 2, name
        for digest in pair:
            assert len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), name
        assert pair[0] != pair[1], name


def test_graded_spread_map_follows_the_rank_key_and_leaves_the_default_alone() -> None:
    """Experiment 095: the best ordering key needs a graded map, not a two-level block.

    Cross-donor agreement beats |effect| as a FULL ordering (+0.00678 vs +0.00448) and
    loses as a block promotion (+0.00341 vs +0.00583), because a block preserves the
    emitted p-order inside itself. The graded map assigns spread by rank across the whole
    gated axis so the realised p-order can follow the key.
    """
    effect = np.asarray([0.0, 0.9, -0.5, 0.2, 0.05])
    b_ctrl = np.log1p(np.asarray([5e4, 4e4, 3e4, 2e4, 1e4]))

    # the default is untouched by the new parameters
    np.testing.assert_array_equal(
        M.emission_spread_vector(
            effect, b_ctrl, min_cpm=0.0, spread_map="graded", spread_lo=0.0, spread_hi=0.9
        ),
        np.full(effect.size, M.SPREAD),
    )

    got = M.emission_spread_vector(
        effect, b_ctrl, min_cpm=0.0, tail_spread=0.25,
        spread_map="graded", spread_lo=0.0, spread_hi=0.6,
    )
    # monotone in |effect| rank: the biggest effect gets spread_lo, the smallest spread_hi
    order = np.argsort(-np.abs(effect), kind="stable")
    assert got[order[0]] == 0.0
    assert got[order[-1]] == 0.6
    assert np.all(np.diff(got[order]) > 0), "graded spread must be strictly ranked"

    # an explicit rank_score overrides |effect| entirely
    rank_score = np.asarray([0.1, 0.2, 0.9, 0.0, 0.5])
    keyed = M.emission_spread_vector(
        effect, b_ctrl, min_cpm=0.0, tail_spread=0.25, rank_score=rank_score,
        spread_map="graded", spread_lo=0.0, spread_hi=0.6,
    )
    assert keyed[2] == 0.0, "highest rank_score must get the lowest spread"
    assert keyed[3] == 0.6, "lowest rank_score must get the highest spread"
    ko = np.argsort(-rank_score, kind="stable")
    assert np.all(np.diff(keyed[ko]) > 0)


def test_graded_spread_map_fails_closed() -> None:
    effect = np.asarray([0.0, 0.9, -0.5])
    b_ctrl = np.log1p(np.asarray([5e4, 4e4, 3e4]))
    for kwargs in (
        {"spread_map": "nope"},
        {"spread_map": "graded", "spread_lo": -0.1},
        {"spread_map": "graded", "spread_hi": 1.5},
        {"spread_map": "graded", "spread_lo": 0.7, "spread_hi": 0.2},
        {"spread_map": "graded", "rank_score": np.asarray([1.0, 2.0])},
    ):
        try:
            M.emission_spread_vector(effect, b_ctrl, min_cpm=0.0, tail_spread=0.25, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {kwargs}")


def test_donor_agreement_is_geometric_magnitude_gated_on_sign_unanimity() -> None:
    """The key experiment 095 measured at +0.00678, the best buildable ordering found."""
    import types

    model = types.SimpleNamespace(
        names=["a", "b", "c"],
        target_genes=np.asarray(["g0", "g1", "g2", "g3"]),
        tcols={n: np.arange(4) for n in "abc"},
        dcols={n: np.arange(4) for n in "abc"},
        pert_index={"a": {"T": 0}, "b": {"T": 0}, "c": {}},
        _renorm={"a": 1.0, "b": 1.0, "c": 1.0},
        s={
            "a": np.asarray([[2.0, -3.0, 4.0, 0.0]]),
            "b": np.asarray([[8.0, -12.0, -1.0, 5.0]]),
            "c": np.asarray([[1.0, 1.0, 1.0, 1.0]]),
        },
    )
    got = M.MultiDonorModel.donor_agreement(model, "T")
    # g0: both +, geometric mean sqrt(2*8) = 4
    # g1: both -, sqrt(3*12) = 6
    # g2: signs disagree -> 0.  g3: a is exactly 0 -> not unanimous -> 0
    np.testing.assert_allclose(got, [4.0, 6.0, 0.0, 0.0], rtol=1e-9, atol=1e-9)

    # a perturbation only one donor carries has nothing to agree about
    model.pert_index = {"a": {"T": 0}, "b": {}, "c": {}}
    np.testing.assert_array_equal(
        M.MultiDonorModel.donor_agreement(model, "T"), np.zeros(4)
    )


def test_overdispersion_preserves_panel_totals_and_is_off_by_default() -> None:
    """The exact-total repair is what makes the emission channel free on pds/mse/nmae.

    Over-dispersion widens the scatter by orders of magnitude, which is the whole point --
    experiment 095 measured it as the only emission lever that moves the DE ordering -- so
    the guarantee that each gene's panel total survives it is the one that must not break.
    """
    profile = np.asarray([0.5, 0.3, 0.15, 0.05])
    libs = np.full(400, 100_000)
    totals = np.floor(400 * 100_000 * profile).astype(np.int64)

    def emit(phi, seed=0):
        v, c, ip = M.emit_counts(
            profile, libs, 0.1, np.random.default_rng(seed),
            mode="quasi", spread=np.ones(4), overdispersion=phi,
        )
        blk = np.zeros((400, 4))
        for r in range(400):
            sl = slice(ip[r], ip[r + 1])
            blk[r, c[sl]] = v[sl]
        return blk

    base = emit(0.0)
    np.testing.assert_array_equal(base.sum(0).astype(np.int64), totals)

    for phi in (0.1, 0.45, 2.0):
        blk = emit(phi)
        np.testing.assert_array_equal(
            blk.sum(0).astype(np.int64), totals,
            err_msg=f"overdispersion {phi} broke the exact-total repair",
        )
        assert (blk >= 0).all(), f"overdispersion {phi} produced a negative count"
        # and it must actually widen the scatter, or it is not doing anything
        assert blk.std(0).mean() > 5 * base.std(0).mean(), f"phi={phi} did not disperse"

    # per-gene phi: each gene gets its own dispersion, and phi=0 genes stay at the floor
    per_gene = M.emit_counts(
        profile, libs, 0.1, np.random.default_rng(5),
        mode="quasi", spread=np.ones(4), overdispersion=np.asarray([0.0, 1.0, 4.0, 8.0]),
    )
    blk = np.zeros((400, 4))
    for r in range(400):
        sl = slice(per_gene[2][r], per_gene[2][r + 1])
        blk[r, per_gene[1][sl]] = per_gene[0][sl]
    np.testing.assert_array_equal(blk.sum(0).astype(np.int64), totals)
    # compare the coefficient of variation, not the raw sd: these genes carry very
    # different panel totals, so sd is dominated by the mean rather than by phi
    cv = blk.std(0) / np.maximum(blk.mean(0), 1e-9)
    assert cv[0] < cv[1] < cv[2] < cv[3], f"per-gene phi is not monotone in CV: {cv}"
    try:
        M.emit_counts(profile, libs, 0.1, np.random.default_rng(0), mode="quasi",
                      spread=np.ones(4), overdispersion=np.asarray([1.0, 2.0]))
    except ValueError:
        pass
    else:
        raise AssertionError("a wrong-length per-gene overdispersion was accepted")

    # default is off: bit-identical to the shipped Poisson branch on the same seed
    np.testing.assert_array_equal(emit(0.0, seed=7), emit(0.0, seed=7))
    v1, c1, i1 = M.emit_counts(profile, libs, 0.1, np.random.default_rng(11),
                               mode="quasi", spread=np.ones(4))
    v2, c2, i2 = M.emit_counts(profile, libs, 0.1, np.random.default_rng(11),
                               mode="quasi", spread=np.ones(4), overdispersion=0.0)
    np.testing.assert_array_equal(v1, v2)
    np.testing.assert_array_equal(c1, c2)
    np.testing.assert_array_equal(i1, i2)


def test_agreement_shrink_vector_is_the_measured_soft_alpha() -> None:
    """`alpha = exp(-gamma * rank / 2000)` -- experiment 096's form, not a re-derivation.

    `mech_miller.py` swept five alphas on both harnesses (hard top-K at 500/2000/8000 and
    this soft exponential at gamma 0.5/1.0) and the soft form at gamma=1.0 won on both.
    Reproducing its rank construction exactly is what lets the 1.17-1.38x it measured be
    quoted for the code that ships.
    """
    rng = np.random.default_rng(0)
    scores = rng.random(500) + 0.5                      # tie-free, so argsort is unambiguous

    # mech_miller.py's own two lines, verbatim
    rank = np.argsort(np.argsort(-scores))
    np.testing.assert_array_equal(
        M.agreement_shrink_vector(scores, 1.0),
        np.exp(-1.0 * rank / M.AGREEMENT_SHRINK_RANK_SCALE),
    )
    np.testing.assert_array_equal(
        M.agreement_shrink_vector(scores, 0.5),
        np.exp(-0.5 * rank / M.AGREEMENT_SHRINK_RANK_SCALE),
    )
    assert M.AGREEMENT_SHRINK_RANK_SCALE == 2000.0
    assert M.AGREEMENT_SHRINK_GAMMA == 0.25
    # 'zero', not 'mean'. The caller specified 'mean' on the reasoning that the metric's
    # scored 0-anchor is the context's mean perturbation response; measurement then showed
    # our own panel mean is uncorrelated with the reference's (cos +0.016 / -0.066), so the
    # premise fails and 'zero' wins on both harnesses at every gamma. See
    # experiments/096-mse-floor/shrink_target.log.
    assert M.AGREEMENT_SHRINK_TOWARD == "zero"

    # monotone increasing in agreement, in (0, 1], and the most-agreed gene is untouched
    got = M.agreement_shrink_vector(np.asarray([1.0, 9.0, 3.0, 7.0]), 1.0, rank_scale=2.0)
    assert got.argmax() == 1 and got.argmin() == 0
    order = np.argsort(-np.asarray([1.0, 9.0, 3.0, 7.0]))
    assert list(np.argsort(-got)) == list(order), "alpha order must follow agreement order"
    assert got.max() == 1.0 and (got > 0).all() and (got <= 1.0).all()

    # TIES SHARE A MIDRANK. The key is exactly 0 on every gene the donors do not agree about
    # -- a median of 16,545 of 18,533 per target -- so if ties took consecutive ranks, the
    # great majority of a row's alphas would be set by position in the gene axis rather than by
    # the key. Equal keys must give equal alphas.
    tied = M.agreement_shrink_vector(np.asarray([0.0, 5.0, 0.0, 0.0]), 1.0, rank_scale=2.0)
    assert tied[1] == 1.0
    assert tied[0] == tied[2] == tied[3], f"tied genes got different alphas: {tied}"
    # the shared rank is the MEAN of the positions the tie spans, here 1, 2, 3 -> 2
    np.testing.assert_allclose(tied[0], np.exp(-1.0 * 2.0 / 2.0), rtol=0, atol=1e-15)
    np.testing.assert_array_equal(
        tied, M.agreement_shrink_vector(np.asarray([0.0, 5.0, 0.0, 0.0]), 1.0, rank_scale=2.0)
    )
    # and permuting a tied block cannot move any alpha
    perm = M.agreement_shrink_vector(np.asarray([0.0, 0.0, 5.0, 0.0]), 1.0, rank_scale=2.0)
    assert perm[2] == 1.0 and perm[0] == perm[1] == perm[3] == tied[0]

    # a target fewer than two donors carry has nothing to agree about: no evidence, no
    # shrink. Falling back to a self-referential key here is what experiment 096 measured
    # as useless-to-harmful, so the all-ones is deliberate, not a degenerate case.
    np.testing.assert_array_equal(
        M.agreement_shrink_vector(np.zeros(6), 1.0), np.ones(6)
    )
    # gamma = 0 is an exact no-op at any key
    np.testing.assert_array_equal(
        M.agreement_shrink_vector(np.asarray([3.0, 1.0, 2.0]), 0.0), np.ones(3)
    )

    for bad in ({"gamma": -0.5}, {"gamma": float("nan")}, {"gamma": float("inf")}):
        try:
            M.agreement_shrink_vector(np.asarray([1.0, 2.0]), **bad)
        except ValueError:
            continue
        raise AssertionError(f"agreement_shrink_vector accepted {bad}")
    for bad in (0.0, -1.0, float("nan")):
        try:
            M.agreement_shrink_vector(np.asarray([1.0, 2.0]), 1.0, rank_scale=bad)
        except ValueError:
            continue
        raise AssertionError(f"agreement_shrink_vector accepted rank_scale={bad}")
    try:
        M.agreement_shrink_vector(np.ones((2, 3)), 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("agreement_shrink_vector accepted a 2-D agreement matrix")


def test_agreement_shrink_is_off_by_default_and_shrinks_the_least_agreed_genes() -> None:
    """OFF must be bit-identical; ON must move the low-agreement genes and only those.

    Gene 0 is the ANCHOR: zero effect and zero centre, so it is invariant under every
    setting and its realised log-ratio to the control is exactly the renormalisation
    constant. Subtracting it recovers the PRE-normalisation log effect from the emitted
    composition, which is the only way to assert direction here. `profile_from_effect`
    renormalises to a composition, so naive per-gene monotonicity does NOT hold: holding
    the top-agreement gene at alpha=1 while every other gene shrinks moves the total, and
    the L2 distance from the shrunk profile to the control is NOT monotone in gamma. That
    is the compositional constraint, not a bug, and asserting on it would fail.
    """
    b_ctrl = np.log1p(M.TS * np.asarray([0.20, 0.22, 0.18, 0.15, 0.12, 0.08, 0.05]))
    effect = np.asarray([0.0, 1.3, -1.0, 0.8, -0.6, 0.45, -0.3])
    agreement = np.asarray([0.5, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0])     # anchor ranks LAST
    center = np.zeros_like(effect)                                  # anchor centre is 0 too

    # (a) OFF is bit-identical to omitting the arguments entirely.
    base = M.profile_from_effect(effect, b_ctrl)
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, gene_scale=None), base
    )
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, gene_scale=np.ones(effect.size)), base
    )
    # and alpha == 1 is bit-identical even WITH a centre -- the interpolation is written
    # as `a*e + (1-a)*c` precisely so the untouched genes reproduce the unflagged build
    # exactly rather than merely closely.
    np.testing.assert_array_equal(
        M.profile_from_effect(
            effect, b_ctrl, gene_scale=np.ones(effect.size),
            shrink_center=np.asarray([0.3, -0.2, 0.1, 0.4, -0.5, 0.2, 0.0]),
        ),
        base,
    )
    parsed = M.build_parser().parse_args(["--output", "unused.h5ad"])
    assert parsed.shrink_key == "none"
    assert parsed.shrink_gamma is None
    assert parsed.shrink_toward is None

    # (b) ON changes the profile, and in the expected direction.
    ctrl = M.profile_from_effect(np.zeros_like(effect), b_ctrl)

    def realised(gamma, toward_mean):
        alpha = M.agreement_shrink_vector(agreement, gamma, rank_scale=2.0)
        profile = M.profile_from_effect(
            effect, b_ctrl, gene_scale=alpha,
            shrink_center=center if toward_mean else None,
        )
        lfc = np.log(profile / ctrl)
        return profile, lfc - lfc[0]        # the anchor removes log(Zc/Zp) exactly

    for toward_mean in (True, False):
        on, _ = realised(1.0, toward_mean)
        assert not np.array_equal(on, base), "the shrinkage did not change the profile"
        mags = np.vstack([np.abs(realised(g, toward_mean)[1])
                          for g in (0.0, 0.25, 0.5, 1.0, 2.0, 8.0)])
        # gene 1 has rank 0, so alpha == 1 at every gamma and it must not move at all
        np.testing.assert_allclose(mags[:, 1], mags[0, 1], rtol=0, atol=1e-12)
        # every other gene's realised log effect shrinks strictly and monotonically
        for j in (2, 3, 4, 5, 6):
            column = mags[:, j]
            assert all(x > y for x, y in zip(column, column[1:])), (
                f"gene {j} not monotone under gamma (toward_mean={toward_mean}): {column}"
            )
        # and never past the control or beyond the unshrunk value
        assert (mags <= np.abs(realised(0.0, toward_mean)[1])[None, :] + 1e-12).all()
        assert (mags >= -1e-12).all()

    # (c) a low-agreement gene shrinks MORE than a high-agreement one. Asserted on a block
    # of genes with identical control abundance and identical effect, so the shared
    # renormalisation constant cannot reorder them.
    flat = np.log1p(M.TS * np.full(6, 1 / 6))
    block = M.agreement_shrink_vector(np.arange(6.0, 0.0, -1.0), 1.0, rank_scale=2.0)
    assert all(x > y for x, y in zip(block, block[1:]))
    up = M.profile_from_effect(np.full(6, 0.9), flat, gene_scale=block)
    assert all(x > y for x, y in zip(up, up[1:])), f"positive effects not ordered: {up}"
    down = M.profile_from_effect(np.full(6, -0.9), flat, gene_scale=block)
    assert all(x < y for x, y in zip(down, down[1:])), f"negative effects not ordered: {down}"

    # a gene the CPM gate zeroed must STAY at exactly zero, centre or no centre: the centre
    # is re-gated, or `(1 - alpha) * mbar_g` would resurrect a gene no metric can see.
    gated_off = np.log1p(M.TS * np.asarray([0.5, 0.5 - 1e-9, 1e-9]))
    assert (1e6 * np.expm1(gated_off) / np.expm1(gated_off).sum() > M.MIN_CPM).tolist() == [
        True, True, False
    ]
    with_center = M.profile_from_effect(
        np.asarray([0.4, -0.3, 5.0]), gated_off,
        gene_scale=np.asarray([1.0, 0.5, 0.0]),
        shrink_center=np.asarray([0.1, 0.2, 9.0]),
    )
    zeroed = M.profile_from_effect(
        np.asarray([0.4, -0.3, 0.0]), gated_off,
        gene_scale=np.asarray([1.0, 0.5, 0.0]),
        shrink_center=np.asarray([0.1, 0.2, 0.0]),
    )
    np.testing.assert_array_equal(with_center, zeroed)

    for bad in ({"gene_scale": np.ones(3)},
                {"gene_scale": np.ones(effect.size), "shrink_center": np.zeros(3)},
                {"shrink_center": np.zeros(effect.size)}):
        try:
            M.profile_from_effect(effect, b_ctrl, **bad)
        except ValueError:
            continue
        raise AssertionError(f"profile_from_effect accepted {list(bad)}")


def test_shrink_toward_mean_lands_on_the_baseline_and_zero_lands_on_the_control() -> None:
    """The two targets are different estimators, and the metric prices them differently.

    `expr_mse_unbiased_capped_norm`'s scored 0-anchor is the context's mean perturbation
    response assigned identically to every perturbation -- NOT the control. The board has
    both ends: team "Kun Lab" submitted an explicit zero-effect prediction and scored raw
    1.0424, against a scored baseline of raw ~0.9879. So a fully shrunk gene must land on
    `mbar`, not on 0, and this test pins which anchor each option reaches.
    """
    b_ctrl = np.log1p(M.TS * np.asarray([0.30, 0.24, 0.20, 0.14, 0.08, 0.04]))
    effect = np.asarray([1.1, -0.8, 0.6, -0.45, 0.3, -0.2])
    mbar = np.asarray([0.25, -0.10, 0.15, 0.05, -0.20, 0.02])
    alpha = M.agreement_shrink_vector(np.arange(6.0, 0.0, -1.0), 1.0, rank_scale=2.0)

    # the James-Stein form, against an independent computation of it
    np.testing.assert_allclose(
        M.profile_from_effect(effect, b_ctrl, gene_scale=alpha, shrink_center=mbar),
        M.profile_from_effect(mbar + alpha * (effect - mbar), b_ctrl),
        rtol=0, atol=1e-15,
    )
    # ...and toward zero is the plain multiply
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, gene_scale=alpha),
        M.profile_from_effect(alpha * effect, b_ctrl),
    )

    # fully shrunk: 'mean' reaches the baseline profile, 'zero' reaches the control. These
    # are different profiles, which is the whole point of the flag.
    zeros = np.zeros(6)
    np.testing.assert_allclose(
        M.profile_from_effect(effect, b_ctrl, gene_scale=zeros, shrink_center=mbar),
        M.profile_from_effect(mbar, b_ctrl), rtol=0, atol=1e-15,
    )
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, gene_scale=zeros),
        M.profile_from_effect(zeros, b_ctrl),
    )
    assert not np.allclose(
        M.profile_from_effect(mbar, b_ctrl), M.profile_from_effect(zeros, b_ctrl)
    ), "the fixture's mbar is degenerate; the two anchors must differ"

    # --effect-scale, --shrink-key and VCC_ALPHA are three knobs on ONE quantity. The centre
    # is scaled with the effect so the two post-gate knobs commute rather than fighting.
    np.testing.assert_allclose(
        M.profile_from_effect(
            effect, b_ctrl, effect_scale=0.6, gene_scale=alpha, shrink_center=mbar
        ),
        M.profile_from_effect(0.6 * effect, b_ctrl, gene_scale=alpha,
                              shrink_center=0.6 * mbar),
        rtol=0, atol=1e-15,
    )


def test_shrink_flags_are_rejected_without_the_agreement_key() -> None:
    """A modifier with no key to modify is an error, not a silent no-op.

    Same contract as `--keep-spread` without `--keep-top`: the flags that only mean
    something under `--shrink-key agreement` must abort the build rather than let an
    operator believe a shrinkage was applied when none was. `run()` reaches this check
    before it opens any model input, so this costs one pass over the frozen source hashes.
    """
    # --help must RENDER. argparse runs every help string through %-formatting, so a bare "%"
    # in one (e.g. "keeps 0.17% of its effect") raises TypeError and takes down --help for the
    # whole builder -- which is what an operator runs first on a fresh worker. This caught a
    # real one; it costs a millisecond to keep.
    rendered = M.build_parser().format_help()
    for flag in ("--shrink-key", "--shrink-gamma", "--shrink-toward",
                 "--shrink-pert-key", "--shrink-pert-gamma", "--shrink-pert-floor"):
        assert flag in rendered, f"{flag} is missing from --help"

    cases = [
        (["--shrink-gamma", "0.5"], "without --shrink-key agreement"),
        (["--shrink-toward", "mean"], "--shrink-toward does nothing"),
        (["--shrink-toward", "zero"], "--shrink-toward does nothing"),
        (["--shrink-pert-gamma", "2.0"], "without --shrink-pert-key donor-cos"),
        (["--shrink-pert-floor", "0.2"], "without --shrink-pert-key donor-cos"),
        # each modifier is rejected against ITS OWN key, so the other key being on is not a
        # licence for it -- this is the case a single "any shrinkage is on" check would miss
        (["--shrink-pert-key", "donor-cos", "--shrink-gamma", "0.5"],
         "without --shrink-key agreement"),
        (["--shrink-key", "agreement", "--shrink-pert-floor", "0.2"],
         "without --shrink-pert-key donor-cos"),
    ]
    for extra, expected in cases:
        args = M.build_parser().parse_args(
            ["--output", str(HERE / "never-written.h5ad"), *extra]
        )
        try:
            M.run(args)
        except ValueError as error:
            assert expected in str(error), f"{extra}: {error}"
        else:
            raise AssertionError(f"run() accepted {extra}")
    assert not (HERE / "never-written.h5ad").exists()


def test_shrink_threading_is_end_to_end_and_off_is_byte_identical() -> None:
    """The requirement is on the ARTIFACT, not on one helper: OFF must reproduce the bytes.

    `profile_from_effect` being bit-identical is necessary but not sufficient -- the flags
    also have to reach it, and a mis-threaded keyword would be invisible to a unit test of
    the helper. This drives the real writer on a synthetic three-context bundle and compares
    sha256 of the emitted .h5ad, which is the guarantee the v11 byte-reproduction check
    actually depends on. It also pins that alpha == 1 with a centre reproduces those same
    bytes, which is why the interpolation is written `a*e + (1-a)*c` rather than
    `c + a*(e-c)`.
    """
    rng = np.random.default_rng(3)
    n_genes = M.N_GENES
    genes = np.asarray([f"G{i}" for i in range(n_genes)])
    spec = M.ChallengeSpec(
        genes=genes, perts=np.asarray(["G7", "G11"]),
        contexts=["A", "B", "C"], cells_per_pert=4, manifest={},
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle = root / "bundle"
        bundle.mkdir()
        for context in spec.contexts:
            # 40 dense columns per cell keeps every context's CPM gate non-degenerate
            dense = M.sp.csr_matrix(
                (np.ones(30 * 40), (np.repeat(np.arange(30), 40), np.tile(np.arange(40), 30))),
                shape=(30, n_genes),
            )
            noise = M.sp.random(
                30, n_genes, density=0.02, format="csr",
                random_state=ord(context), data_rvs=lambda n: rng.integers(1, 40, n),
            )
            controls = M.ad.AnnData(X=(dense + noise).tocsr())
            controls.var_names = genes
            controls.write_h5ad(bundle / f"context_{context}.h5ad")

        effects = rng.normal(0.0, 0.5, (2, n_genes))
        alpha = np.vstack([
            M.agreement_shrink_vector(rng.random(n_genes), 1.0) for _ in range(2)
        ])
        mbar = effects.mean(axis=0)

        def build(tag, **kwargs):
            sink = io.StringIO()
            with contextlib.redirect_stdout(sink):
                report = M.write_submission_from_effects(
                    spec, effects, bundle, root / f"{tag}.h5ad",
                    max_rss_gib=0.0, **kwargs,
                )
            return report["sha256"]

        off = build("off")
        assert off == build("off_explicit", shrink_scales=None, shrink_center=None)
        # alpha == 1 everywhere reproduces the unflagged bytes even WITH a centre
        assert off == build("alpha_one", shrink_scales=np.ones_like(alpha),
                            shrink_center=mbar)
        toward_zero = build("zero", shrink_scales=alpha)
        toward_mean = build("mean", shrink_scales=alpha, shrink_center=mbar)
        assert len({off, toward_zero, toward_mean}) == 3, (
            "the shrinkage did not reach the emitted artifact, or the two targets agree"
        )

        for bad in ({"shrink_scales": np.ones((3, n_genes))},
                    {"shrink_center": mbar},
                    {"shrink_scales": np.ones_like(alpha), "shrink_center": np.ones(5)}):
            try:
                build("rejected", **bad)
            except ValueError:
                continue
            raise AssertionError(f"write_submission_from_effects accepted {list(bad)}")


def test_a_key_with_no_ordering_is_a_no_op_and_the_on_target_gene_is_protected() -> None:
    """The shipped-path bug this family was measured through, and its structural fix.

    `donor_agreement` returns an ALL-ZERO row whenever fewer than two donors carry the target.
    Under the original descending sort with ties broken by gene index, such a row produced
    `exp(-gamma * g / 2000)` indexed by POSITION IN THE GENE AXIS -- gene 0 kept full amplitude
    and gene 18,000 was shrunk to nothing, for no reason but where it sat in the array. On the
    real challenge panel that was 111 of 900 rows, 12% of what we would ship, and those rows
    were scored as though the shrinkage had helped, which inflated the family's estimated gain
    from +0.0055 to +0.0102.

    The invariant asserted here is the general one, not the special case: a row whose key
    carries NO ORDERING must be bit-identical to that row with the shrinkage off. That covers
    the all-zero row, the constant non-zero row (same defect, never counted), and the
    single-gene axis, and it is what a special case for `not scores.any()` would have missed.
    """
    b_ctrl = np.log1p(M.TS * np.asarray([0.30, 0.24, 0.20, 0.14, 0.08, 0.04]))
    effect = np.asarray([1.1, -0.8, 0.6, -0.45, 0.3, -0.2])

    for label, key in (("all zero (fewer than two donors)", np.zeros(6)),
                       ("constant non-zero", np.full(6, 0.7)),
                       ("constant negative", np.full(6, -1.5))):
        for gamma in (0.5, 1.0, 4.0):
            alpha = M.agreement_shrink_vector(key, gamma)
            np.testing.assert_array_equal(
                alpha, np.ones(6), err_msg=f"{label} at gamma={gamma} was not a no-op"
            )
            # the invariant that matters is on the PROFILE, and it is bit-identity
            for center in (None, np.asarray([0.2, -0.1, 0.3, 0.0, -0.4, 0.1])):
                np.testing.assert_array_equal(
                    M.profile_from_effect(effect, b_ctrl, gene_scale=alpha,
                                          shrink_center=center),
                    M.profile_from_effect(effect, b_ctrl),
                    err_msg=f"{label} at gamma={gamma} moved the profile",
                )

    # one distinct value is no ordering; two is an ordering and must still shrink
    assert not np.array_equal(
        M.agreement_shrink_vector(np.asarray([0.0, 0.0, 1.0, 0.0, 0.0, 0.0]), 1.0), np.ones(6)
    )

    # `protect` pins one gene at exactly 1.0. cell_eval2 excludes each row's own target
    # coordinate from all six scored members, so this is free -- but without it a knockdown
    # whose gene ranks poorly is shrunk 9x at gamma=0.5 and 99x at gamma=1.0, a defect no
    # number we check would catch.
    key = np.asarray([9.0, 0.1, 8.0, 0.2, 7.0, 0.3])
    bare = M.agreement_shrink_vector(key, 1.0, rank_scale=1.0)
    assert bare[1] < 0.02, "fixture must put the protected gene deep in the tail"
    for position in range(6):
        guarded = M.agreement_shrink_vector(key, 1.0, rank_scale=1.0, protect=position)
        assert guarded[position] == 1.0, f"gene {position} was not protected"
        others = [i for i in range(6) if i != position]
        np.testing.assert_array_equal(guarded[others], bare[others])
    for bad in (-1, 6, 99):
        try:
            M.agreement_shrink_vector(key, 1.0, protect=bad)
        except ValueError:
            continue
        raise AssertionError(f"protect={bad} was accepted")


def test_per_perturbation_weight_is_rank01_of_the_donor_cosine() -> None:
    """Experiment 097's per-perturbation half: `z_p = rank01(mean pairwise donor cosine)^g`.

    `rank01(x) = (argsort(argsort(x)) + 0.5) / n`, ascending, so the best-agreed target keeps
    ~all of its amplitude and the worst keeps 0.5/n of it. Taken over the CARRIED targets only,
    with uncarried targets left at z = 1 -- the same "no cross-donor evidence, no shrink" policy
    `agreement_shrink_vector` applies per gene, and deliberately NOT 097's own choice of
    dropping them to the bottom of the rank scale, which is a self-referential penalty on
    targets we have no information about.
    """
    score = np.asarray([0.30, -0.10, 0.05, 0.20, 0.40])
    carried = np.ones(5, dtype=bool)

    # against an independent computation of 097's rank01, unfloored
    rank01 = (np.argsort(np.argsort(score)).astype(float) + 0.5) / score.size
    for gamma in (0.5, 1.0, 2.0, 3.0):
        np.testing.assert_allclose(
            M.pert_shrink_vector(score, carried, gamma, floor=0.0), rank01 ** gamma,
            rtol=0, atol=1e-15,
        )
    assert M.PERT_SHRINK_GAMMA == 1.0
    assert M.PERT_SHRINK_FLOOR == 0.3
    assert M.PERT_SHRINK_MISSING == "damp"

    # THE FLOOR GOES INSIDE THE EXPONENT: z = max(rank01, floor)**gamma. 097's grid is reported
    # in those units, so clamping after the exponent instead would put the builder at 0.3 where
    # their table says 0.09 and every number they quote would silently disagree.
    np.testing.assert_allclose(
        M.pert_shrink_vector(score, carried, 2.0, floor=0.3),
        np.maximum(rank01, 0.3) ** 2, rtol=0, atol=1e-15,
    )
    assert M.pert_shrink_vector(score, carried, 2.0, floor=0.3).min() == 0.3 ** 2
    assert not np.allclose(
        M.pert_shrink_vector(score, carried, 2.0, floor=0.3),
        np.maximum(M.pert_shrink_vector(score, carried, 2.0, floor=0.0), 0.3),
    ), "the floor must be inside the exponent, not applied after it"

    # monotone in the key, in (0, 1], best target highest
    weight = M.pert_shrink_vector(score, carried, 1.0, floor=0.0)
    assert list(np.argsort(weight)) == list(np.argsort(score))
    assert (weight > 0).all() and (weight <= 1.0).all()
    # gamma = 0 is an exact no-op under every policy
    for policy in ("damp", "neutral", "keep"):
        np.testing.assert_array_equal(
            M.pert_shrink_vector(score, carried, 0.0, missing=policy), np.ones(5)
        )

    # the three missing-target policies. `keep` gives an evidence-free target the TOP weight,
    # which is why it is not the default: z is relative, so the rows we know least about would
    # outrank every row we know something about.
    mask = np.asarray([True, True, True, False, False])
    keep_z = M.pert_shrink_vector(score, mask, 1.0, floor=0.0, missing="keep")
    assert keep_z[3] == 1.0 and keep_z[4] == 1.0
    sub_rank = (np.argsort(np.argsort(score[:3])).astype(float) + 0.5) / 3
    np.testing.assert_allclose(keep_z[:3], sub_rank, rtol=0, atol=1e-15)
    assert keep_z[3] > keep_z[:3].max(), "keep makes evidence-free targets rank top"

    neutral_z = M.pert_shrink_vector(score, mask, 1.0, floor=0.0, missing="neutral")
    np.testing.assert_allclose(neutral_z[:3], sub_rank, rtol=0, atol=1e-15)
    assert neutral_z[3] == neutral_z[4] == float(np.median(sub_rank))

    # `damp` scores them 0.0 and ranks the whole panel together -- note that is NOT the bottom
    # when real cosines go negative, which they do on the challenge panel
    damp_z = M.pert_shrink_vector(score, mask, 1.0, floor=0.0, missing="damp")
    # both uncarried targets score exactly 0.0, so they TIE and share a midrank -- ordinal
    # ranks would hand them 0.3 and 0.5 purely by target index, which is the same defect the
    # per-gene half had. Reference: scores [0.30, -0.10, 0.05, 0.0, 0.0] sort ascending to
    # -0.10 | 0.0, 0.0 | 0.05 | 0.30 at positions 0 | 1, 2 | 3 | 4, so the tied pair takes 1.5.
    np.testing.assert_allclose(
        damp_z, (np.asarray([4.0, 0.0, 3.0, 1.5, 1.5]) + 0.5) / 5, rtol=0, atol=1e-15
    )
    assert damp_z[3] == damp_z[4], "tied uncarried targets must get the same weight"
    assert damp_z[3] > damp_z[1], "a 0.0 score outranks a genuinely negative cosine"

    # the same structural guard as the per-gene half: no distinct values, no ordering
    for policy in ("damp", "neutral", "keep"):
        np.testing.assert_array_equal(
            M.pert_shrink_vector(np.full(5, 0.3), carried, 1.0, missing=policy), np.ones(5)
        )
        np.testing.assert_array_equal(
            M.pert_shrink_vector(score, np.zeros(5, dtype=bool), 1.0, missing=policy),
            np.ones(5),
        )
    # ties share a midrank, so equal keys give equal weights
    tied = M.pert_shrink_vector(np.asarray([0.1, 0.1, 0.9, 0.1]), np.ones(4, bool), 1.0,
                                floor=0.0)
    assert tied[0] == tied[1] == tied[3] < tied[2]

    try:
        M.pert_shrink_vector(score, carried, 1.0, missing="nope")
    except ValueError:
        pass
    else:
        raise AssertionError("an unknown --shrink-pert-missing policy was accepted")

    for bad in ({"gamma": -1.0}, {"gamma": float("nan")}):
        try:
            M.pert_shrink_vector(score, carried, **bad)  # noqa
        except ValueError:
            continue
        raise AssertionError(f"pert_shrink_vector accepted {bad}")
    for bad in (-0.1, 1.5, float("nan")):
        try:
            M.pert_shrink_vector(score, carried, 1.0, floor=bad)
        except ValueError:
            continue
        raise AssertionError(f"pert_shrink_vector accepted floor={bad}")
    try:
        M.pert_shrink_vector(score, np.ones(3, dtype=bool), 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("a mismatched carried mask was accepted")


def test_the_two_shrink_halves_compose_and_each_is_reachable_alone() -> None:
    """Per-gene only, per-perturbation only, and both together must each behave.

    The composition is `z[:, None] * alpha`, so a per-row scalar cannot reorder genes within a
    row -- which is why the per-perturbation half leaves a row's DE membership untouched while
    the per-gene half can move it.
    """
    b_ctrl = np.log1p(M.TS * np.asarray([0.30, 0.24, 0.20, 0.14, 0.08, 0.04]))
    effect = np.asarray([1.1, -0.8, 0.6, -0.45, 0.3, -0.2])
    alpha = M.agreement_shrink_vector(np.asarray([6.0, 1.0, 5.0, 2.0, 4.0, 3.0]), 1.0,
                                      rank_scale=2.0)
    z = 0.4

    gene_only = M.profile_from_effect(effect, b_ctrl, gene_scale=alpha)
    pert_only = M.profile_from_effect(effect, b_ctrl, gene_scale=np.full(6, z))
    both = M.profile_from_effect(effect, b_ctrl, gene_scale=z * alpha)
    base = M.profile_from_effect(effect, b_ctrl)
    assert len({base.tobytes(), gene_only.tobytes(), pert_only.tobytes(),
                both.tobytes()}) == 4, "the three arms must be distinguishable"

    # per-perturbation alone is exactly a global amplitude cut on that row, so it is the same
    # thing --effect-scale does -- which is the check that it cannot reorder genes
    np.testing.assert_allclose(
        pert_only, M.profile_from_effect(effect, b_ctrl, effect_scale=z), rtol=0, atol=1e-15
    )
    # and z = 1 on every target reduces to the per-gene arm exactly
    np.testing.assert_array_equal(
        M.profile_from_effect(effect, b_ctrl, gene_scale=1.0 * alpha), gene_only
    )


def test_pairwise_agreement_matches_unanimity_at_two_donors_and_grades_more_at_three() -> None:
    """The property that makes experiment 097's harness numbers apply to the production key.

    Both harnesses use two donors. With exactly two donors, "mean over unordered pairs" has one
    pair and "unanimity across every carrying donor" IS that pair's sign test, so the two
    definitions are the same function and no harness measurement can distinguish them. They
    separate only at production's third donor, and there `donor_agreement` collapses -- it
    demands unanimity across all three and zeroes any gene missing from any axis, grading a
    median 10.6% of panel genes on three-donor targets against 54.5% on two-donor ones, so the
    targets with LESS evidence are damped LESS. The pairwise form removes that inversion.

    `src/vcc2026` is frozen and `--keep-key agreement` is board-validated on the unanimity
    method, so this is an alternative key for --shrink-key only and the shipped one is untouched.
    """
    import types

    def make(names, rows):
        return types.SimpleNamespace(
            names=list(names),
            target_genes=np.asarray(["g0", "g1", "g2", "g3"]),
            tcols={n: np.arange(4) for n in names},
            dcols={n: np.arange(4) for n in names},
            pert_index={n: {"T": 0} for n in names},
            _renorm={n: 1.0 for n in names},
            s={n: np.asarray([rows[n]]) for n in names},
        )

    # --- two donors: bit-identical, which is what makes 097's numbers transferable
    two = make(["a", "b"], {"a": [2.0, -3.0, 4.0, 0.0], "b": [8.0, -12.0, -1.0, 5.0]})
    np.testing.assert_array_equal(
        M.donor_agreement_pairwise(two, "T"),
        M.MultiDonorModel.donor_agreement(two, "T"),
    )
    np.testing.assert_allclose(
        M.donor_agreement_pairwise(two, "T"), [4.0, 6.0, 0.0, 0.0], rtol=1e-9, atol=1e-9
    )

    # --- three donors: g2 is graded by the a/b pair even though c disagrees there, which is
    # exactly the gene unanimity throws away
    three = make(["a", "b", "c"], {"a": [2.0, -3.0, 4.0, 0.0],
                                   "b": [8.0, -12.0, 9.0, 5.0],
                                   "c": [1.0, -1.0, -1.0, 1.0]})
    unanimous = M.MultiDonorModel.donor_agreement(three, "T")
    pairwise = M.donor_agreement_pairwise(three, "T")
    assert unanimous[2] == 0.0, "unanimity must discard a gene one donor dissents on"
    assert pairwise[2] > 0.0, "the pairwise key must keep the pair that agrees"
    # and it grades strictly more genes, which is the whole point
    assert int((pairwise > 0).sum()) > int((unanimous > 0).sum())
    # g0 and g1 are unanimous, so every pair contributes and the mean is a plain average
    np.testing.assert_allclose(
        pairwise[0], np.mean([np.sqrt(2 * 8), np.sqrt(2 * 1), np.sqrt(8 * 1)]),
        rtol=1e-9, atol=1e-9,
    )

    # --- fewer than two donors: all zero, same contract as donor_agreement
    lone = make(["a", "b"], {"a": [1.0, 2.0, 3.0, 4.0], "b": [1.0, 1.0, 1.0, 1.0]})
    lone.pert_index = {"a": {"T": 0}, "b": {}}
    np.testing.assert_array_equal(M.donor_agreement_pairwise(lone, "T"), np.zeros(4))
    np.testing.assert_array_equal(
        M.donor_agreement_pairwise(lone, "T"),
        M.MultiDonorModel.donor_agreement(lone, "T"),
    )

    # --- the flag exposes it, and --keep-key deliberately does NOT get the choice
    parser = M.build_parser()
    assert parser.parse_args(
        ["--output", "x.h5ad", "--shrink-key", "agreement-pairwise"]
    ).shrink_key == "agreement-pairwise"
    for action in parser._actions:
        if action.dest == "keep_key":
            assert "agreement-pairwise" not in action.choices, (
                "--keep-key is board-validated on the frozen unanimity key; leave it alone"
            )


def test_cd4_target_f_is_the_099_frame_and_solves_the_fraction_in_closed_form() -> None:
    """`f` is 099's ratio of median row norms; `CD4_FRACTION` is a differently-scaled knob."""
    # Two carrying rows and one that carries nothing, so the median lands on a carrying row --
    # the same situation as the real panel (276 of 300 complete).
    cd4 = M.EXP062.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray(
            [[3.0, 4.0], [3.0, 4.0], [0.0, 0.0]], dtype=np.float64
        ),
        gene_votes=np.asarray([[1, 1], [1, 1], [0, 0]], dtype=np.uint8),
        target_condition_count=np.asarray([3, 3, 1], dtype=np.uint8),
        condition_details=[],
    )
    # Whitened base with median row norm exactly 50.
    whitened = np.asarray(
        [[30.0, 40.0], [30.0, 40.0], [30.0, 40.0]], dtype=np.float64
    )

    fraction, frame = M.solve_cd4_fraction(whitened, cd4)
    # residual median row norm at unit fraction is 5, base is 50, so f = 0.1 per unit fraction.
    assert frame["residual_median_row_norm_at_unit_fraction"] == 5.0
    assert frame["whitened_base_median_row_norm"] == 50.0
    assert frame["f_per_unit_fraction"] == 0.1

    # Unset target: the constant is used untouched and `f` is only reported.
    assert M.CD4_TARGET_F is None or isinstance(M.CD4_TARGET_F, float)
    if M.CD4_TARGET_F is None:
        assert fraction == M.CD4_FRACTION
        assert frame["fraction_source"] == "VCC_CD4_FRACTION"
        # The knob and the frame are NOT the same number -- this is the 099 warning, asserted.
        assert frame["effective_f_at_cd4_fraction"] == M.CD4_FRACTION * 0.1
        assert frame["effective_f_at_cd4_fraction"] != M.CD4_FRACTION

    # The solve is closed-form and exact: applying it reproduces the requested `f`.
    for target in (0.05, 0.10, 0.15):
        solved = target / frame["f_per_unit_fraction"]
        added = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0) * solved
        realised = float(np.median(np.linalg.norm(added, axis=1))) / 50.0
        assert abs(realised - target) < 1e-12


def test_cd4_target_f_env_override_is_validated_and_fails_closed() -> None:
    """A target outside (0, 1) must stop the build at import, not silently clamp."""
    for bad in ("0", "-0.1", "1.0", "2", "nonsense"):
        env = {**os.environ, "VCC_CD4_TARGET_F": bad}
        proc = subprocess.run(
            [sys.executable, "-c", f"import runpy; runpy.run_path({str(MODULE_PATH)!r})"],
            env=env,
            capture_output=True,
            text=True,
        )
        assert proc.returncode != 0, f"VCC_CD4_TARGET_F={bad} was accepted"
        assert (
            "VCC_CD4_TARGET_F" in proc.stderr or "could not convert" in proc.stderr
            or "ValueError" in proc.stderr
        ), proc.stderr


def _mass_fix_fixture(seed: int = 11):
    """A control profile that is a genuine `log1p(TS * p)` object, as the real one is.

    This matters: the mass constraint is `sum(expm1(b)) == TS`, and it holds for the control
    EXACTLY only because `b_ctrl` is built as `log1p(TS * totals / totals.sum())`. A fixture
    with arbitrary `b` would make the constraint unreachable and the test vacuous.
    """
    rng = np.random.default_rng(seed)
    counts = rng.gamma(0.35, 900.0, size=400) + 1e-3
    b_ctrl = np.log1p(M.TS * counts / counts.sum())
    effect = rng.normal(0.0, 0.35, size=b_ctrl.size)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / expression.sum()
    gate = (cpm > M.MIN_CPM) & (b_ctrl > 0)
    assert gate.sum() >= 20 and not gate.all(), "fixture must gate some genes and keep some"
    assert abs(float(np.expm1(b_ctrl).sum()) - M.TS) < 1e-6
    return b_ctrl, effect, gate


def test_mass_fix_is_off_by_default_and_leaves_the_renderer_untouched() -> None:
    """Unset `VCC_MASS_FIX` must reproduce every prior build bit-for-bit."""
    assert M.MASS_FIX is None, "the module default must be off"
    b_ctrl, effect, _ = _mass_fix_fixture()
    baseline = M.profile_from_effect(effect, b_ctrl)
    assert np.array_equal(baseline, M.EXP061.profile_from_effect(effect, b_ctrl))

    saved = M.MASS_FIX
    try:
        M.MASS_FIX = "minnorm"
        fixed = M.profile_from_effect(effect, b_ctrl)
    finally:
        M.MASS_FIX = saved
    assert not np.array_equal(baseline, fixed), "the flag must actually reach the renderer"
    # and turning it back off restores the byte-exact original
    assert np.array_equal(baseline, M.profile_from_effect(effect, b_ctrl))


def test_mass_fix_makes_the_read_back_effect_exact() -> None:
    """The whole point: with the constraint met, the scorer reads back what we designed.

    Without the fix the renormalisation costs a measurable slice of the effect; with it the
    round trip `effect -> profile -> log1p(TS * p) - b_ctrl` is the identity.
    """
    b_ctrl, effect, gate = _mass_fix_fixture()
    gated = np.where(gate, effect, 0.0)

    plain = M.profile_from_effect(effect, b_ctrl)
    read_back = np.log1p(M.TS * plain) - b_ctrl
    drift = float(np.linalg.norm(read_back - gated) / np.linalg.norm(gated))
    assert drift > 1e-3, f"fixture does not exercise the renderer (drift {drift:.2e})"

    h = M.mass_consistent_effect(gated, b_ctrl, gate)
    # the renderer clips, so the constrained quantity is the CLIPPED mass
    assert abs(float(np.clip(np.expm1(b_ctrl + h), 0.0, None).sum()) - M.TS) <= 1e-6
    saved = M.MASS_FIX
    try:
        M.MASS_FIX = "minnorm"
        fixed = M.profile_from_effect(effect, b_ctrl)
    finally:
        M.MASS_FIX = saved
    assert abs(float(fixed.sum()) - 1.0) < 1e-12
    exact = np.log1p(M.TS * fixed) - b_ctrl
    # A gene the renderer clips to zero profile reads back as -b_ctrl whatever we do; those
    # are excluded, and the assertion is exactness on every gene that is actually renderable.
    live = (b_ctrl + h) > 0.0
    assert live.sum() > 0.9 * b_ctrl.size, "fixture clipped too much to be meaningful"
    assert float(np.abs(exact[live] - h[live]).max()) < 1e-9, "read-back must be exact"


def test_mass_fix_stays_inside_the_cpm_gate_and_is_idempotent() -> None:
    """A gene the CPM gate zeroed must stay at exactly zero, and re-applying is a no-op."""
    b_ctrl, effect, gate = _mass_fix_fixture()
    gated = np.where(gate, effect, 0.0)
    h = M.mass_consistent_effect(gated, b_ctrl, gate)
    assert np.array_equal(h[~gate], np.zeros(int((~gate).sum())))
    again = M.mass_consistent_effect(h, b_ctrl, gate)
    assert float(np.abs(again - h).max()) < 1e-12


def test_mass_fix_minnorm_is_the_smaller_correction_and_bad_modes_fail_closed() -> None:
    """`minnorm` is minimum-L2 by construction, so it must beat `const` on norm.

    `const` satisfies the identical scalar constraint and experiment 156 measured it as
    decisively WORSE on pds and transfer cosine -- it is kept only as that control, so this
    asserts the two are genuinely different objects rather than an alias.
    """
    b_ctrl, effect, gate = _mass_fix_fixture()
    gated = np.where(gate, effect, 0.0)
    a = M.mass_consistent_effect(gated, b_ctrl, gate, mode="minnorm")
    c = M.mass_consistent_effect(gated, b_ctrl, gate, mode="const")
    for h in (a, c):
        assert abs(float(np.clip(np.expm1(b_ctrl + h), 0.0, None).sum()) - M.TS) <= 1e-6
    assert np.linalg.norm(a - gated) < np.linalg.norm(c - gated)

    try:
        M.mass_consistent_effect(gated, b_ctrl, gate, mode="sideways")
    except ValueError:
        pass
    else:
        raise AssertionError("an unknown mass-fix mode must raise")


def load_module_with_env(**env: str):
    """Re-exec the builder with `env` applied, so its import-time constants are re-read.

    `VCC_DONOR_WEIGHTS` and `VCC_SIGN_REPAIR` are read once at module import, exactly as
    `VCC_ALPHA` and `VCC_MASS_FIX` are, so an in-process reload is the only way to test the
    ON state without a subprocess. The environment is restored afterwards.
    """
    saved = {key: os.environ.get(key) for key in env}
    os.environ.update(env)
    try:
        spec = importlib.util.spec_from_file_location(
            "test_target_make_submission_pm_cd4_variant", MODULE_PATH
        )
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    return module


class WeightRecordingModel:
    """Records the weight dict `build_v9_pre_whiten_effects` hands to `MultiDonorModel`."""

    names = list(M.PRODUCTION_DONORS)
    target_genes = np.asarray(["G0", "G1", "G2", "G3"])
    effect_scale = 6.817
    pert_index = {name: {"T0": 0, "T1": 1} for name in M.PRODUCTION_DONORS}
    # one distinguishable unit vector per donor, so a weighted mean is readable off the row
    rows = {
        "k562_gwps": np.asarray([1.0, 0.0, 0.0, 0.0]),
        "orion_hct116": np.asarray([0.0, 1.0, 0.0, 0.0]),
        "orion_hek293t": np.asarray([0.0, 0.0, 1.0, 0.0]),
    }

    def __init__(self) -> None:
        self.seen: list[dict[str, float]] = []

    def effect(self, target, alpha, beta, *, weights):
        self.seen.append(dict(weights))
        carriers = [n for n in self.names if target in self.pert_index[n]]
        total = sum(weights[n] for n in carriers) or 1.0
        out = np.zeros(4, dtype=np.float64)
        for n in carriers:
            out += (weights[n] / total) * self.rows[n]
        return alpha * out


def test_donor_weights_default_is_the_shipped_dict_and_env_overrides_it() -> None:
    """OFF must be the literal `{name: 1.0}`; ON must be the requested mix, in donor order."""
    assert M.DONOR_WEIGHTS is None, "the module default must be unset"
    model = WeightRecordingModel()
    assert M.donor_weights(model) == {name: 1.0 for name in M.PRODUCTION_DONORS}

    variant = load_module_with_env(VCC_DONOR_WEIGHTS="2,1,1")
    assert variant.DONOR_WEIGHTS == {
        "k562_gwps": 2.0,
        "orion_hct116": 1.0,
        "orion_hek293t": 1.0,
    }
    assert variant.donor_weights(model) == variant.DONOR_WEIGHTS
    # and the default is still the default in this process
    assert M.donor_weights(model) == {name: 1.0 for name in M.PRODUCTION_DONORS}

    # the order is positional, so a mis-assignment would show up as K562 not being the
    # up-weighted donor -- assert the mapping, not just the multiset of values
    lopsided = load_module_with_env(VCC_DONOR_WEIGHTS="1,1,3")
    assert lopsided.DONOR_WEIGHTS["orion_hek293t"] == 3.0
    assert lopsided.DONOR_WEIGHTS["k562_gwps"] == 1.0


def test_donor_weights_reach_the_effect_and_renormalise_over_carriers() -> None:
    """OFF is bit-identical to the shipped path; ON moves the mix by the intended amount."""
    targets = ("T0", "T1")
    base_model = WeightRecordingModel()
    base = M.build_v9_pre_whiten_effects(base_model, targets, alpha=1.0, beta=0.0)
    np.testing.assert_array_equal(base, np.tile([1 / 3, 1 / 3, 1 / 3, 0.0], (2, 1)))
    assert base_model.seen[0] == {name: 1.0 for name in M.PRODUCTION_DONORS}

    variant = load_module_with_env(VCC_DONOR_WEIGHTS="2,1,1")
    armed_model = WeightRecordingModel()
    armed = variant.build_v9_pre_whiten_effects(armed_model, targets, alpha=1.0, beta=0.0)
    np.testing.assert_allclose(armed, np.tile([0.5, 0.25, 0.25, 0.0], (2, 1)))
    assert armed_model.seen[0] == {"k562_gwps": 2.0, "orion_hct116": 1.0,
                                   "orion_hek293t": 1.0}
    # the weights are relative: (4,2,2) must land on exactly the same effect as (2,1,1)
    scaled = load_module_with_env(VCC_DONOR_WEIGHTS="4,2,2")
    np.testing.assert_array_equal(
        armed, scaled.build_v9_pre_whiten_effects(WeightRecordingModel(), targets,
                                                  alpha=1.0, beta=0.0)
    )
    assert not np.array_equal(base, armed), "the switch must actually move the mixture"


def test_donor_weights_env_is_validated_and_fails_closed() -> None:
    """A malformed or degenerate mix must stop the build at import, not silently default."""
    for bad in ("2,1", "2,1,1,1", "2,1,x", "-1,1,1", "0,0,0", ""):
        env = {**os.environ, "VCC_DONOR_WEIGHTS": bad}
        if bad == "":
            # the empty string is the documented OFF state, not an error
            proc = subprocess.run(
                [sys.executable, "-c",
                 f"import runpy; m = runpy.run_path({str(MODULE_PATH)!r}); "
                 "assert m['DONOR_WEIGHTS'] is None"],
                env=env, capture_output=True, text=True,
            )
            assert proc.returncode == 0, proc.stderr
            continue
        proc = subprocess.run(
            [sys.executable, "-c", f"import runpy; runpy.run_path({str(MODULE_PATH)!r})"],
            env=env, capture_output=True, text=True,
        )
        assert proc.returncode != 0, f"VCC_DONOR_WEIGHTS={bad} was accepted"
        assert (
            "VCC_DONOR_WEIGHTS" in proc.stderr
            or "could not convert" in proc.stderr
        ), proc.stderr


def _sign_repair_fixture(seed: int = 5):
    """A pre-whitening consensus and a `whitened` matrix that genuinely flips some signs.

    The fixture is only meaningful if the two disagree on some entries, agree on others, and
    have zeros on both sides -- otherwise `flip` is trivially all-true or all-false and the
    test would pass for a repair that did the wrong thing.
    """
    rng = np.random.default_rng(seed)
    pre = rng.normal(size=(6, 9))
    whitened = pre * rng.choice([-1.0, 1.0], size=pre.shape)
    pre[0, 0] = 0.0            # consensus silent -> must not be touched
    whitened[1, 1] = 0.0       # whitened exactly zero -> must not be touched
    columns = np.asarray([0, 1, 2, 3, 4, 5], dtype=np.int64)   # 6 on-axis, 3 off-axis
    flip = (np.sign(whitened[:, columns]) * np.sign(pre[:, columns])) < 0.0
    assert 0 < flip.sum() < flip.size, "fixture must both flip and not flip"
    return pre, whitened, columns, flip


def test_sign_repair_is_off_by_default_and_only_touches_flipped_on_axis_entries() -> None:
    assert M.SIGN_REPAIR is None, "the module default must be off"
    pre, whitened, columns, flip = _sign_repair_fixture()
    out, audit = M.sign_repair_to_consensus(whitened, pre, columns)

    # 1. magnitudes are preserved exactly -- this is a sign repair, not a shrinkage
    np.testing.assert_array_equal(np.abs(out), np.abs(whitened))
    # 2. every off-axis coordinate is untouched
    off = np.asarray([6, 7, 8])
    np.testing.assert_array_equal(out[:, off], whitened[:, off])
    # 3. on-axis, exactly the disagreeing entries moved, and they moved to the consensus sign
    moved = out[:, columns] != whitened[:, columns]
    np.testing.assert_array_equal(moved, flip)
    both = (np.sign(whitened[:, columns]) != 0) & (np.sign(pre[:, columns]) != 0)
    np.testing.assert_array_equal(
        np.sign(out[:, columns])[both], np.sign(pre[:, columns])[both]
    )
    # 4. the zero rows the fixture planted are still zero / still untouched
    assert out[0, 0] == whitened[0, 0]
    assert out[1, 1] == 0.0
    # 5. and the whole thing is NOT a no-op, which is what a broken mask would look like
    assert not np.array_equal(out, whitened)
    assert audit["flipped_entries"] == int(flip.sum())
    assert np.isclose(audit["flipped_fraction_on_axis"], flip.mean())
    assert audit["off_axis_columns_untouched"] == 3
    # applying it twice is a no-op: the signs already agree
    again, audit2 = M.sign_repair_to_consensus(out, pre, columns)
    np.testing.assert_array_equal(again, out)
    assert audit2["flipped_entries"] == 0


def test_sign_repair_switch_is_off_by_default_and_reaches_the_effect_matrix() -> None:
    """OFF: byte-identical effects and no `sign_repair` key. ON: the intended flips."""
    model = ParityModel()
    targets = ("T0", "T1")
    columns = np.asarray([0, 1, 2], dtype=np.int64)
    # a whitener that inverts the middle coordinate, so it must flip exactly that sign
    transform = LinearWhitener(np.diag([1.0, -1.0, 1.0]))

    base, base_audit = M.build_production_effects(
        model, targets, transform, columns, np.asarray([3], dtype=np.int64),
        pm=None, cd4=None,
    )
    assert "sign_repair" not in base_audit, "the OFF state must not add a meta key"

    variant = load_module_with_env(VCC_SIGN_REPAIR="cons")
    assert variant.SIGN_REPAIR == "cons"

    class _Whitener(LinearWhitener):
        def apply(self, values, gamma, lam, renorm=False):
            assert gamma == variant.WHITEN_GAMMA and lam == variant.WHITEN_LAM
            return np.asarray(values, dtype=np.float64) @ self.matrix

    armed, armed_audit = variant.build_production_effects(
        model, targets, _Whitener(np.diag([1.0, -1.0, 1.0])), columns,
        np.asarray([3], dtype=np.int64), pm=None, cd4=None,
    )
    assert armed_audit["sign_repair"]["enabled"] is True
    assert armed_audit["sign_repair"]["mode"] == "cons"
    assert armed_audit["sign_repair"]["flipped_entries"] == 2   # column 1, both targets

    pre = np.vstack([model.effect(t, variant.ALPHA, variant.BETA,
                                  weights={n: 1.0 for n in model.names})
                     for t in targets]).astype(np.float64)
    # magnitude identical to the un-repaired build, sign realigned to the consensus on-axis
    np.testing.assert_allclose(np.abs(armed), np.abs(base))
    live = (np.sign(pre[:, columns]) != 0) & (np.sign(base[:, columns]) != 0)
    np.testing.assert_array_equal(
        np.sign(armed[:, columns])[live], np.sign(pre[:, columns])[live]
    )
    assert not np.array_equal(armed, base), "the switch must move the effect matrix"
    # the OFF module in this process is unchanged
    again, again_audit = M.build_production_effects(
        model, targets, transform, columns, np.asarray([3], dtype=np.int64),
        pm=None, cd4=None,
    )
    np.testing.assert_array_equal(again, base)
    assert again_audit == base_audit


def test_sign_repair_env_is_validated_and_fails_closed() -> None:
    for bad in ("1", "true", "yes", "flip", "zero"):
        env = {**os.environ, "VCC_SIGN_REPAIR": bad}
        proc = subprocess.run(
            [sys.executable, "-c", f"import runpy; runpy.run_path({str(MODULE_PATH)!r})"],
            env=env, capture_output=True, text=True,
        )
        assert proc.returncode != 0, f"VCC_SIGN_REPAIR={bad} was accepted"
        assert "VCC_SIGN_REPAIR" in proc.stderr, proc.stderr


def _synthetic_bundle(root: Path, module, *, all_gated: bool, seed: int = 3,
                      contexts=("A", "B", "C"), perts=("G7", "G11")):
    """Three synthetic control contexts on the official gene dimension.

    `all_gated=False` reproduces the sparse bundle the shrink-threading test uses, where
    only ~40 genes clear the CPM gate, so a `--resample-keep-top` below that count leaves
    a genuine resampled tail. `all_gated=True` gives every gene enough mass to clear the
    gate, which is the only way `--resample-keep-top = n_genes` can select EVERY gene and
    so the only setting under which `resample_mix` can be bit-identical to `quasi`.
    """
    rng = np.random.default_rng(seed)
    n_genes = module.N_GENES
    genes = np.asarray([f"G{i}" for i in range(n_genes)])
    bundle = root / "bundle"
    bundle.mkdir()
    for context in contexts:
        if all_gated:
            # every gene expressed in every cell -> cpm ~ 54 > MIN_CPM = 10 for all
            base = np.ones((30, n_genes), dtype=np.float64)
            base += rng.integers(0, 3, base.shape)
            matrix = module.sp.csr_matrix(base)
        else:
            dense = module.sp.csr_matrix(
                (
                    np.ones(30 * 40),
                    (np.repeat(np.arange(30), 40), np.tile(np.arange(40), 30)),
                ),
                shape=(30, n_genes),
            )
            noise = module.sp.random(
                30, n_genes, density=0.02, format="csr",
                random_state=ord(context), data_rvs=lambda n: rng.integers(1, 40, n),
            )
            matrix = (dense + noise).tocsr()
        controls = module.ad.AnnData(X=matrix)
        controls.var_names = genes
        controls.write_h5ad(bundle / f"context_{context}.h5ad")
    # A real bundle also ships the three panel files, and the context adapter reads its
    # context ids out of the manifest rather than assuming A/B/C -- so the synthetic bundle
    # has to carry them too, or nothing that follows the manifest can be exercised here.
    (bundle / "manifest.json").write_text(json.dumps({
        "season": "test", "partition": "synthetic", "panel_id": "synthetic-1",
        "contexts": list(contexts), "pert_col": "target_gene", "context_col": "context",
        "control_label": "non-targeting", "n_genes": int(n_genes),
        "n_constructs": len(perts), "cells_per_pert": 4,
    }, indent=1) + "\n")
    (bundle / "gene_names.csv").write_text(
        "gene_name\n" + "\n".join(map(str, genes)) + "\n"
    )
    (bundle / "pert_counts.csv").write_text(
        "target_gene,n_cells\n" + "".join(f"{p},4\n" for p in perts)
    )
    spec = module.ChallengeSpec(
        genes=genes, perts=np.asarray(list(perts)),
        contexts=list(contexts), cells_per_pert=4, manifest={},
    )
    return bundle, spec


def _build_artifact(module, spec, bundle, path: Path, **kwargs) -> dict:
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        return module.write_submission_from_effects(
            spec, kwargs.pop("effects"), bundle, path, max_rss_gib=0.0, **kwargs
        )


def _panel_totals_by_block(path: Path, spec) -> dict:
    """Per-(context, perturbation) per-gene panel totals, read off the delivered file."""
    import h5py

    per_block = {}
    with h5py.File(path, "r") as handle:
        data = handle["X"]["data"][:]
        indices = handle["X"]["indices"][:]
        indptr = handle["X"]["indptr"][:]
    block = 0
    for context in spec.contexts:
        for pert in spec.perts:
            start = block * spec.cells_per_pert
            stop = start + spec.cells_per_pert
            lo, hi = int(indptr[start]), int(indptr[stop])
            totals = np.zeros(len(spec.genes), dtype=np.int64)
            np.add.at(totals, indices[lo:hi].astype(np.int64),
                      data[lo:hi].astype(np.int64))
            per_block[(str(context), str(pert))] = totals
            block += 1
    return per_block


def _load_pre_emitter_builder(markers=("--emitter", "resample_keep_top")):
    """The newest committed builder revision that predates the given flag.

    A same-module comparison cannot answer "did this patch change the DEFAULT emission",
    because whatever the patch does is on both sides of it. This resolves the real
    pre-change object out of git and compares the emitted ARTIFACT against it, which is
    the only form of the check that reads differently when the thing it tests is absent.
    It fails rather than skips if the revision cannot be found.
    """
    rel = "scripts/make_submission_pm_cd4.py"
    revs = subprocess.run(
        ["git", "-C", str(ROOT), "log", "--format=%H", "-n", "60", "--", rel],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    for rev in revs:
        source = subprocess.run(
            ["git", "-C", str(ROOT), "show", f"{rev}:{rel}"],
            capture_output=True, text=True, check=True,
        ).stdout
        if any(marker in source for marker in markers):
            continue
        target = HERE / f".pre_emitter_{rev[:12]}.py"
        target.write_text(source)
        try:
            spec = importlib.util.spec_from_file_location(
                "test_target_pre_emitter_builder", target
            )
            assert spec is not None and spec.loader is not None
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
        finally:
            target.unlink()
        return rev, module
    raise AssertionError(
        f"no committed revision of {rel} predates {markers}; the bit-identity check "
        "cannot be made and must not be silently skipped"
    )


def test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte() -> None:
    """Rule 10: the shipped configuration, through the patched builder, unchanged.

    `--emitter` defaults to `quasi`, and the claim is that an unflagged build is
    bit-identical to every prior build. This drives the REAL writer on a synthetic
    three-context bundle with both the patched module and the last pre-change revision
    out of git, and compares sha256 of the emitted .h5ad -- at the shipped emission
    policy (`--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`)
    as well as at the bare default.
    """
    rev, old = _load_pre_emitter_builder()
    assert old.N_GENES == M.N_GENES
    shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(17).normal(0.0, 0.5, (2, M.N_GENES))

        for tag, kwargs in (("bare", {}), ("shipped", shipped)):
            new = _build_artifact(M, spec, bundle, root / f"new_{tag}.h5ad",
                                  effects=effects, **kwargs)
            was = _build_artifact(old, old_spec, bundle, root / f"old_{tag}.h5ad",
                                  effects=effects, **kwargs)
            assert new["sha256"] == was["sha256"], (
                f"the patched builder changed the DEFAULT emission at {tag!r} against "
                f"revision {rev[:12]}: {new['sha256']} != {was['sha256']}"
            )
            assert new["nnz"] == was["nnz"] and new["nnz"] > 0
            # explicitly passing the default must land on the same bytes too
            explicit = _build_artifact(
                M, spec, bundle, root / f"explicit_{tag}.h5ad", effects=effects,
                emitter="quasi", **kwargs
            )
            assert explicit["sha256"] == new["sha256"]

        # ...and the check is not vacuous: the resampling law on the SAME inputs must
        # produce a different artifact, or the comparison above proves nothing.
        moved = _build_artifact(
            M, spec, bundle, root / "moved.h5ad", effects=effects,
            emitter="resample_mix", resample_keep_top=20, **shipped
        )
        assert moved["sha256"] != new["sha256"], "resample_mix did not reach the artifact"


def test_resample_mix_over_the_whole_gated_axis_is_the_quasi_law() -> None:
    """`param_mask` all-True must collapse `resample_mix` onto the shipped quasi law.

    `emit_counts` draws the parametric block FIRST, so with every gene parametric the
    resampled half is computed and then discarded, the splice is a no-op, and the
    exact-total repair has nothing to repair. If that identity does not hold, the mask is
    wired backwards or the two laws disagree on the block they are supposed to share.
    """
    profile = np.asarray([0.4, 0.3, 0.2, 0.1])
    libs = np.asarray([20_000.0, 21_000.0, 19_500.0, 20_400.0] * 100)
    control = np.random.default_rng(2).poisson(
        np.outer(libs, profile)
    ).astype(np.float64)
    spread = np.full(4, 0.25)

    quasi = M.emit_counts(profile, libs, 0.1, np.random.default_rng(9),
                          mode="quasi", spread=spread, overdispersion=4.0)
    mixed = M.emit_counts(profile, libs, 0.1, np.random.default_rng(9),
                          mode="resample_mix", spread=spread, overdispersion=4.0,
                          control_counts=control, control_profile=profile,
                          param_mask=np.ones(4, dtype=bool))
    for a, b in zip(quasi, mixed):
        np.testing.assert_array_equal(a, b)
    # and it is not vacuous: one gene off the parametric mask changes the emission
    partial = M.emit_counts(profile, libs, 0.1, np.random.default_rng(9),
                            mode="resample_mix", spread=spread, overdispersion=4.0,
                            control_counts=control, control_profile=profile,
                            param_mask=np.asarray([True, True, True, False]))
    assert not np.array_equal(quasi[0], partial[0])

    # End to end through the builder, on a bundle where every gene clears the CPM gate so
    # `--resample-keep-top = n_genes` really does select all of them.
    #
    # ONE emission block, deliberately. The writer drives all 900 blocks off a SINGLE
    # emission stream, and `resample_mix` draws `n * G` extra uniforms per block for its
    # stochastic rounding even when the splice discards them, so the streams diverge from
    # the second block onward. That is a stream offset, not a law change -- the per-gene
    # panel totals are a function of the picks (a separate stream) and the profile only,
    # which the next test checks on every block of a three-context build.
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, contexts=("A",),
                                         perts=("G7",))
        effects = np.random.default_rng(23).normal(0.0, 0.5, (1, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25,
                       overdispersion=4.0)
        base = _build_artifact(M, spec, bundle, root / "quasi.h5ad", effects=effects,
                               **shipped)
        allpar = _build_artifact(
            M, spec, bundle, root / "allpar.h5ad", effects=effects,
            emitter="resample_mix", resample_keep_top=M.N_GENES, **shipped
        )
        assert allpar["sha256"] == base["sha256"], (
            "resample_mix with every gene parametric is not the quasi law"
        )
        # not vacuous: leaving a tail off the mask must move the same single block
        tailed = _build_artifact(
            M, spec, bundle, root / "tailed.h5ad", effects=effects,
            emitter="resample_mix", resample_keep_top=100, **shipped
        )
        assert tailed["sha256"] != base["sha256"]
        # the mask really did cover the whole axis (not merely the gate's cap by accident)
        with M.h5py.File(bundle / "context_A.h5ad", "r") as handle:
            shape = tuple(handle["X"].attrs["shape"])
            matrix = M.sp.csr_matrix(
                (handle["X"]["data"][:], handle["X"]["indices"][:],
                 handle["X"]["indptr"][:]), shape=shape,
            )
        totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
        b_ctrl = np.log1p(M.TS * totals / totals.sum())
        mask = M.emission_param_mask(effects[0], b_ctrl, resample_keep_top=M.N_GENES)
        assert mask.all(), f"only {int(mask.sum())} of {mask.size} genes were parametric"


def test_resample_mix_delivers_the_quasi_panel_totals_on_every_block() -> None:
    """The identity `pds` and `mse`'s profile half ride on, checked on the ARTIFACTS.

    `emit_counts` splices the two laws and THEN repairs each gene's panel total to
    `floor(n * L * profile)`, and `L` is a function of the picked control libraries only,
    which the emission law does not touch. So the delivered per-gene panel total of every
    (context, perturbation) block must be identical between the two laws -- which is what
    makes the pseudobulk, and therefore `pds` and the profile half of `mse`, invariant by
    construction rather than by estimate.
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        effects = np.random.default_rng(31).normal(0.0, 0.5, (2, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25,
                       overdispersion=4.0)
        quasi_path = root / "q.h5ad"
        mix_path = root / "m.h5ad"
        quasi = _build_artifact(M, spec, bundle, quasi_path, effects=effects, **shipped)
        mixed = _build_artifact(M, spec, bundle, mix_path, effects=effects,
                                emitter="resample_mix", resample_keep_top=20, **shipped)
        assert quasi["sha256"] != mixed["sha256"], "the two laws emitted the same file"

        a = _panel_totals_by_block(quasi_path, spec)
        b = _panel_totals_by_block(mix_path, spec)
        assert set(a) == set(b) and len(a) == 6
        moved = 0
        for key, totals in a.items():
            np.testing.assert_array_equal(
                totals, b[key], err_msg=f"panel totals moved on {key}"
            )
            moved += int(totals.sum() > 0)
        assert moved == 6, "every block must carry counts, or the check is empty"


def test_rescaled_control_cells_are_rejected_before_they_reach_the_emitter() -> None:
    """Experiment 177 section 0.1, as a guard rather than as a comment.

    Rescaling the donated control cells to a nominal depth while `lib_sizes` still holds
    the real one makes the exact-total repair peel ~80%% of every gene's column one count
    at a time. That is variance-DESTROYING, it called 10,658 of 10,779 genes significant
    in 177's reproduction, and it inverts the sign of the whole arm -- while looking
    entirely reasonable from the outside.
    """
    rng = np.random.default_rng(4)
    libs = rng.integers(8_000, 30_000, 40).astype(np.float64)
    block = np.zeros((40, 60))
    for i, lib in enumerate(libs):
        idx = rng.integers(0, 60, int(lib))
        np.add.at(block[i], idx, 1.0)
    np.testing.assert_array_equal(block.sum(1), libs)
    M.assert_raw_control_counts(block, libs, where="ok")          # the passing case

    scaled = block * (100_000.0 / libs)[:, None]                  # 177 section 0.1
    cpm = 1e6 * block / libs[:, None]
    for bad, why in (
        (scaled, "rescaled to a nominal depth"),
        (cpm, "CPM-normalised"),
        (block + 0.5, "non-integer"),
        (-block, "negative"),
        (np.where(block > 0, np.inf, block), "non-finite"),
    ):
        try:
            M.assert_raw_control_counts(bad, libs, where="bad")
        except ValueError:
            continue
        raise AssertionError(f"a {why} control block was accepted")
    try:
        M.assert_raw_control_counts(block[:10], libs, where="bad")
    except ValueError:
        pass
    else:
        raise AssertionError("a control block with the wrong row count was accepted")

    # the guard has to be REACHED: a builder run with the resampling law exercises it, and
    # it is the only path that supplies control_counts at all.
    source = MODULE_PATH.read_text()
    assert "assert_raw_control_counts(" in source.split("def assert_raw_control_counts")[1]


def test_emitter_flags_default_to_the_shipped_law_and_fail_closed() -> None:
    parsed = M.build_parser().parse_args(["--output", "unused.h5ad"])
    assert parsed.emitter == "quasi"
    assert parsed.resample_keep_top is None

    # `--resample-keep-top` defaults to `--keep-top`, so the parametric half is exactly
    # the incumbent's keep block and only the tail's law changes.
    b_ctrl = np.log1p(np.asarray([0.0, 1e-6, 8.0, 80.0, 300.0]) * M.TS)
    effect = np.asarray([9.0, 8.0, -0.5, 0.25, 2.0])
    keep = M.emission_spread_vector(effect, b_ctrl, tail_spread=1.0, keep_top=2,
                                    keep_spread=0.25)
    mask = M.emission_param_mask(effect, b_ctrl, resample_keep_top=2)
    np.testing.assert_array_equal(mask, keep == 0.25)
    assert mask.tolist() == [False, False, True, False, True]
    # the gate caps it, exactly as `--keep-top` is capped
    assert M.emission_param_mask(effect, b_ctrl, resample_keep_top=99).sum() == 3
    assert not M.emission_param_mask(effect, b_ctrl, resample_keep_top=0).any()

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        effects = np.zeros((2, M.N_GENES))
        for i, bad in enumerate((
            {"emitter": "resample"},
            {"emitter": "quasi", "resample_keep_top": 10},
            {"emitter": "resample_mix", "resample_keep_top": -1},
        )):
            try:
                _build_artifact(M, spec, bundle, root / f"bad{i}.h5ad",
                                effects=effects, **bad)
            except ValueError:
                continue
            raise AssertionError(f"write_submission_from_effects accepted {bad}")


def test_bundle_file_list_follows_the_manifest_and_not_a_hardcoded_abc() -> None:
    """The bundle half of the frozen table has to be derivable for ANY context set.

    What this reads if the generalisation were absent: `bundle_file_names(["D","E","F"])`
    would either raise or hand back A/B/C, and `frozen_input_paths` would point the D/E/F
    build at `context_A.h5ad`.
    """
    assert M.bundle_file_names(["A", "B", "C"]) == [
        "manifest.json", "gene_names.csv", "pert_counts.csv",
        "context_A.h5ad", "context_B.h5ad", "context_C.h5ad",
    ]
    assert M.bundle_file_names(["D", "E", "F"])[3:] == [
        "context_D.h5ad", "context_E.h5ad", "context_F.h5ad"
    ]
    assert M.bundle_file_names(["Z"])[3:] == ["context_Z.h5ad"]

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle = root / "bundle"
        bundle.mkdir()
        (bundle / "manifest.json").write_text(json.dumps({"contexts": ["D", "E", "F"]}))
        assert M.read_bundle_contexts(bundle) == ("D", "E", "F")
        (bundle / "manifest.json").write_text(json.dumps({"contexts": []}))
        try:
            M.read_bundle_contexts(bundle)
        except ValueError:
            pass
        else:
            raise AssertionError("an empty context list was accepted")

    # and the path map really follows the argument, on every entry
    paths = M.frozen_input_paths(
        Path("/bundle"), ["D", "E", "F"], {"k562_gwps": Path("/k")},
        Path("/pm"), Path("/cd4"), Path("/w.npz"), Path("/w.json"),
    )
    assert {k for k in paths if k.startswith("bundle/")} == {
        "bundle/manifest.json", "bundle/gene_names.csv", "bundle/pert_counts.csv",
        "bundle/context_D.h5ad", "bundle/context_E.h5ad", "bundle/context_F.h5ad",
    }
    assert paths["bundle/context_E.h5ad"] == Path("/bundle/context_E.h5ad")


def test_frozen_bundle_table_fails_closed_on_every_shape_error() -> None:
    """A freeze that can be desynchronised silently is not a freeze."""
    good = json.loads(M.FROZEN_BUNDLE_PATH.read_text())
    assert good["contexts"] == ["A", "B", "C"]
    assert set(good["files"]) == {
        f"bundle/{name}" for name in M.bundle_file_names(good["contexts"])
    }
    # the A/B/C hashes must be preserved EXACTLY across the move out of source
    assert good["files"]["bundle/context_A.h5ad"] == (
        "f22e71968487d3da9402769dcd0c73a47851fad4cd7dd25510283c923031fe3e"
    )
    assert good["files"]["bundle/context_B.h5ad"] == (
        "41557555b3febf6db19fd930b4dd7db88941840077bf985d14dfff959c782380"
    )
    assert good["files"]["bundle/context_C.h5ad"] == (
        "e090e7c7cfdd0006f689214ba38ba0cd2669578e89c9bc325b65f0bd940c02e0"
    )
    assert good["files"]["bundle/manifest.json"] == (
        "2c2b5d425df817ca35e9edf08282d12494ce42913ba824375bf0f98a1cc676f5"
    )
    assert good["files"]["bundle/gene_names.csv"] == (
        "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"
    )
    assert good["files"]["bundle/pert_counts.csv"] == (
        "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
    )
    assert good["panel_support"] == {
        "perturbmatch_replacement_targets": M.EXPECTED_PM_REPLACEMENT_TARGETS,
        "cd4_complete_targets": M.EXPECTED_CD4_COMPLETE_TARGETS,
    }

    def broken(mutate):
        table = json.loads(json.dumps(good))
        mutate(table)
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            path = Path(directory) / "frozen.json"
            path.write_text(json.dumps(table))
            try:
                M.load_frozen_bundle(path)
            except ValueError:
                return
            raise AssertionError(f"a broken frozen table was accepted: {mutate}")

    broken(lambda t: t.__setitem__("contexts", []))
    broken(lambda t: t.__setitem__("contexts", ["A", "A", "C"]))
    broken(lambda t: t.__setitem__("contexts", ["A", "B"]))          # files now over-specify
    broken(lambda t: t.__setitem__("contexts", ["A", "B", "C", "D"]))  # ...and under-specify
    broken(lambda t: t["files"].pop("bundle/context_B.h5ad"))
    broken(lambda t: t["files"].__setitem__("bundle/context_D.h5ad", "0" * 64))
    broken(lambda t: t["files"].__setitem__("bundle/context_A.h5ad", "abc"))
    broken(lambda t: t["files"].__setitem__("bundle/context_A.h5ad", "z" * 64))
    broken(lambda t: t.__setitem__("n_constructs", 0))
    broken(lambda t: t.__setitem__("cells_per_pert", "400"))
    broken(lambda t: t.__setitem__("panel_support", {"cd4_complete_targets": 276}))
    broken(lambda t: t["panel_support"].__setitem__("cd4_complete_targets", "276"))

    # `panel_support: null` is a LEGAL state (a freshly re-pinned panel), not a shape error.
    table = json.loads(json.dumps(good))
    table["panel_support"] = None
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "frozen.json"
        path.write_text(json.dumps(table))
        assert M.load_frozen_bundle(path)["panel_support"] is None


def test_repin_frozen_inputs_moves_the_table_to_a_def_bundle_and_never_builds() -> None:
    """The documented re-pin step, rehearsed on a D/E/F bundle.

    Non-vacuity: before the re-pin the D/E/F bundle's own path set cannot even be verified
    against the A/B/C table (`verify_expected_hashes` reports the mismatch); after it, every
    hash matches the files on disk, `contexts` has moved, and `panel_support` has been
    dropped rather than carried over, because the panel changed.
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False, contexts=("D", "E", "F"))

        # (a) the A/B/C table cannot verify this bundle at all
        try:
            M.verify_expected_hashes(
                {f"bundle/{n}": bundle / n for n in M.bundle_file_names(spec.contexts)},
                {k: v for k, v in M.load_frozen_bundle()["files"].items()},
                label="frozen model input",
            )
        except ValueError as error:
            assert "path set mismatch" in str(error)
        else:
            raise AssertionError("the A/B/C table verified a D/E/F bundle")

        # (b) re-pin into a throwaway table, seeded from the real one
        path = root / "frozen.json"
        path.write_text(M.FROZEN_BUNDLE_PATH.read_text())
        before = M.FROZEN_BUNDLE_PATH.read_bytes()
        sink = io.StringIO()
        with contextlib.redirect_stdout(sink):
            table = M.repin_frozen_inputs(bundle, path)
        assert M.FROZEN_BUNDLE_PATH.read_bytes() == before, (
            "re-pinning a throwaway table rewrote the production one"
        )
        assert table["contexts"] == ["D", "E", "F"]
        assert table["n_constructs"] == len(spec.perts)
        assert table["panel_support"] is None, (
            "panel_support was carried onto a panel it was never measured on"
        )
        assert table["supersedes"] and len(table["supersedes"]) == 64
        assert set(table["files"]) == {
            f"bundle/{n}" for n in M.bundle_file_names(["D", "E", "F"])
        }
        for key, digest in table["files"].items():
            assert digest == M.file_sha256(bundle / key.split("/", 1)[1])

        # (c) the re-pinned table now verifies the bundle it was pinned from...
        observed = M.verify_expected_hashes(
            {f"bundle/{n}": bundle / n for n in M.bundle_file_names(["D", "E", "F"])},
            table["files"], label="frozen model input",
        )
        assert observed == table["files"]

        # (d) ...and one changed byte still aborts
        (bundle / "pert_counts.csv").write_text(
            (bundle / "pert_counts.csv").read_text() + "GX,4\n"
        )
        try:
            M.verify_expected_hashes(
                {f"bundle/{n}": bundle / n for n in M.bundle_file_names(["D", "E", "F"])},
                table["files"], label="frozen model input",
            )
        except ValueError as error:
            assert "hash drift" in str(error)
        else:
            raise AssertionError("a drifted bundle file was accepted")

        # (e) a second re-pin whose panel files are UNCHANGED keeps panel_support
        path2 = root / "frozen2.json"
        seeded = json.loads(json.dumps(table))
        seeded["panel_support"] = {
            "perturbmatch_replacement_targets": 7, "cd4_complete_targets": 9,
        }
        seeded["files"]["bundle/pert_counts.csv"] = M.file_sha256(
            bundle / "pert_counts.csv"
        )
        path2.write_text(json.dumps(seeded))
        with contextlib.redirect_stdout(io.StringIO()):
            again = M.repin_frozen_inputs(bundle, path2)
        assert again["panel_support"] == seeded["panel_support"]


def test_context_effects_default_is_none_and_only_moves_its_own_context() -> None:
    """`VCC_ALPHA=knee`'s only entry point into the writer, through the REAL writer.

    Three claims, each with a reading when the wiring is absent:
      * handing every context the same matrix must be bit-identical to handing none -- if
        `context_effects` were ignored this passes trivially, which is why (b) exists;
      * moving ONE context's matrix must move that context's blocks and no other's -- if
        the writer read the shared matrix, every block would be unchanged;
      * a missing or wrongly shaped context must abort.
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        effects = np.random.default_rng(41).normal(0.0, 0.5, (2, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)

        plain = _build_artifact(M, spec, bundle, root / "plain.h5ad",
                                effects=effects, **shipped)
        same = _build_artifact(
            M, spec, bundle, root / "same.h5ad", effects=effects,
            context_effects={c: effects for c in spec.contexts}, **shipped,
        )
        assert same["sha256"] == plain["sha256"], (
            "an all-identical context_effects map changed the artifact"
        )

        moved_effects = {c: effects for c in spec.contexts}
        moved_effects["B"] = effects * 1.35
        moved = _build_artifact(
            M, spec, bundle, root / "moved.h5ad", effects=effects,
            context_effects=moved_effects, **shipped,
        )
        assert moved["sha256"] != plain["sha256"], "context_effects never reached the writer"
        base_blocks = _panel_totals_by_block(root / "plain.h5ad", spec)
        moved_blocks = _panel_totals_by_block(root / "moved.h5ad", spec)
        for key in base_blocks:
            equal = np.array_equal(base_blocks[key], moved_blocks[key])
            if key[0] == "B":
                assert not equal, f"context B block {key} did not move"
            else:
                assert equal, f"context {key[0]} moved when only B's matrix changed"

        for bad in (
            {"A": effects, "B": effects},
            {"A": effects, "B": effects, "C": effects, "D": effects},
            {"A": effects, "B": effects, "C": effects[:1]},
        ):
            try:
                _build_artifact(M, spec, bundle, root / "bad.h5ad", effects=effects,
                                context_effects=bad, **shipped)
            except ValueError:
                pass
            else:
                raise AssertionError(f"a malformed context_effects map was accepted: {bad}")


def test_alpha_knee_env_is_parsed_and_fails_closed() -> None:
    """`VCC_ALPHA=knee` must be reachable, and the default must stay the shipped number."""
    default = load_module_with_env()
    assert default.ALPHA_MODE == "fixed" and default.ALPHA == 1.0796

    fixed = load_module_with_env(VCC_ALPHA="0.982")
    assert fixed.ALPHA_MODE == "fixed" and fixed.ALPHA == 0.982

    knee = load_module_with_env(VCC_ALPHA="knee")
    assert knee.ALPHA_MODE == "knee" and knee.ALPHA == 1.0796
    anchored = load_module_with_env(VCC_ALPHA="KNEE", VCC_ALPHA_REF="0.982")
    assert anchored.ALPHA_MODE == "knee" and anchored.ALPHA == 0.982

    for env in (
        {"VCC_ALPHA": "kneee"},
        {"VCC_ALPHA": ""},
        {"VCC_ALPHA": "-1"},
        {"VCC_ALPHA_REF": "0.9"},                      # ref without knee
        {"VCC_ALPHA": "0.982", "VCC_ALPHA_REF": "1.0"},
        {"VCC_ALPHA": "knee", "VCC_ALPHA_REF": "nope"},
    ):
        try:
            load_module_with_env(**env)
        except (ValueError, RuntimeError):
            pass
        else:
            raise AssertionError(f"VCC_ALPHA env {env} was accepted")


def test_knee_puts_the_unclamped_budget_ratio_on_one_and_is_not_a_no_op() -> None:
    """The knee rule, end to end on a synthetic context, with its own falsifier.

    `rho = min(1, B/claim)` reads 1.0 both ON the knee and ABOVE it, so it cannot be the
    check. `rho_unclamped = B/claim` can: it is < 1 below the knee and > 1 above it. The
    non-vacuity requirement is that the STARTING ratio is far from 1 -- otherwise "it lands
    on 1" is a statement about the starting point and not about the correction.
    """
    n_targets, n_genes, n_cells = 6, 400, 40
    rng = np.random.default_rng(19)
    b_ctrl = np.log1p(50_000.0 * rng.dirichlet(np.full(n_genes, 4.0)))
    effects = rng.normal(0.0, 0.35, (n_targets, n_genes))
    cols = {i: i for i in range(n_targets)}          # each target owns one gene column
    common = dict(
        min_cpm=5.0, effect_scale=1.0, shrink_scales=None, shrink_center=None,
        rank_scores=None, overdispersion=4.0, cells_per_pert=n_cells, subsample=3,
        spread_kwargs=dict(tail_spread=1.0, keep_top=200, keep_spread=0.25,
                           keep_top_frac=None, spread_map="block",
                           spread_lo=0.0, spread_hi=0.60),
    )

    def frame(cp):
        return M.ControlFrame(context="A", b_ctrl=b_ctrl, median_library=20_000.0,
                              cp=cp, cp_values=np.asarray([cp]), cp_source="synthetic")

    # `C_p` well above the emitted cells' own dispersion, so `min(C_hat, C_p)` does not bind
    # and `claim` is a property of the emission alone. The cap's own branch is checked below.
    cp = 1e6
    probe = M.knee_context_terms(effects, frame(cp), cols, **common)
    assert probe["cap_saturated_rows"] == 0, "the cap bound when it was meant not to"

    # Put the fixture a long way BELOW its knee by shrinking the amplitude: `B` is very
    # nearly quadratic in it, so this solves in one step and the residual is what the
    # correction below then has to remove.
    target_ratio = 0.80
    amplitude = float(np.sqrt(target_ratio * probe["claim"] / probe["B"]))
    effects = effects * amplitude
    before = M.knee_context_terms(effects, frame(cp), cols, **common)
    assert abs(before["rho_unclamped"] - target_ratio) < 0.05, before["rho_unclamped"]
    assert abs(before["rho_unclamped"] - 1.0) > 0.15, "the starting point is already on the knee"
    assert before["rho"] == before["rho_unclamped"]      # clamped only from above

    # The `min(C_hat, C_p)` branch, on the same fixture: a small `C_p` caps every row and
    # `claim` collapses onto the cap exactly.
    small = M.knee_context_terms(effects, frame(before["C_hat_min"] / 2.0), cols, **common)
    assert small["cap_saturated_rows"] == small["cap_rows_measured"]
    assert abs(small["claim"] - small["claim_cap"]) < 1e-12
    assert small["claim"] < before["claim"]

    scale = float(np.sqrt(before["claim"] / before["B"]))
    assert abs(scale - 1.0) > 0.05, "the correction is a no-op on this fixture"
    after = M.knee_context_terms(effects * scale, frame(cp), cols, **common)
    assert abs(after["rho_unclamped"] - 1.0) < 0.02, after["rho_unclamped"]

    # ...and over-correcting by the same factor again lands ABOVE the knee, where `rho`
    # clamps to 1 and `rho_unclamped` does not. That is the asymmetry the check relies on.
    over = M.knee_context_terms(effects * scale * scale, frame(cp), cols, **common)
    assert over["rho_unclamped"] > 1.05 and over["rho"] == 1.0

    # the solver turns the same terms into one alpha per context
    frames = {"A": frame(cp)}
    spec = M.ChallengeSpec(
        genes=np.asarray([f"G{i}" for i in range(n_genes)]),
        perts=np.asarray([f"G{i}" for i in range(n_targets)]),
        contexts=["A"], cells_per_pert=n_cells, manifest={},
    )
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        solved = M.solve_knee_alphas_with_columns(
            spec, effects, frames, cols, alpha_ref=0.982, **common
        )
    assert abs(solved["A"]["knee_scale"] - scale) < 1e-12
    assert abs(solved["A"]["alpha"] - 0.982 * scale) < 1e-12
    assert M.panel_gene_columns(spec) == cols


TESTS = (
    test_disabled_path_matches_v9_effect_semantics,
    test_perturbmatch_target_support_policy_and_source_share,
    test_panel_exclusion_for_operator_and_cd4_residual,
    test_complete_mean3_support_semantics,
    test_rng_streams_are_deterministically_separated,
    test_profile_matches_frozen_v9_gate,
    test_emission_spread_is_the_shipped_scalar_unless_a_tail_is_asked_for,
    test_emission_spread_tail_policy_keeps_the_top_k_gated_effects,
    test_atomic_h5_output_publishes_or_cleans_up,
    test_frozen_constants_match_reused_helpers,
    test_authorization_accepts_only_exploratory_scope,
    test_authorization_rejects_run_prefix_and_gate_spoofing,
    test_authorization_rejects_nonfinite_and_hash_relationship_drift,
    test_authorization_rejects_unsafe_or_duplicate_manifest_paths,
    test_frozen_input_and_helper_hash_drift_fails_before_use,
    test_partial_vcc_failure_never_publishes_final_package,
    test_cloud_contract_has_dual_locks_pinned_sources_and_success_only_vcc,
    test_keep_spread_lowers_the_kept_genes_instead_of_raising_the_tail,
    test_keep_spread_is_validated_and_fails_closed,
    test_whitener_registry_selects_a_pair_and_still_fails_closed,
    test_graded_spread_map_follows_the_rank_key_and_leaves_the_default_alone,
    test_graded_spread_map_fails_closed,
    test_donor_agreement_is_geometric_magnitude_gated_on_sign_unanimity,
    test_overdispersion_preserves_panel_totals_and_is_off_by_default,
    test_effect_scale_is_off_by_default_and_shrinks_toward_the_control,
    test_agreement_shrink_vector_is_the_measured_soft_alpha,
    test_agreement_shrink_is_off_by_default_and_shrinks_the_least_agreed_genes,
    test_shrink_toward_mean_lands_on_the_baseline_and_zero_lands_on_the_control,
    test_shrink_flags_are_rejected_without_the_agreement_key,
    test_shrink_threading_is_end_to_end_and_off_is_byte_identical,
    test_a_key_with_no_ordering_is_a_no_op_and_the_on_target_gene_is_protected,
    test_per_perturbation_weight_is_rank01_of_the_donor_cosine,
    test_the_two_shrink_halves_compose_and_each_is_reachable_alone,
    test_pairwise_agreement_matches_unanimity_at_two_donors_and_grades_more_at_three,
    test_cd4_target_f_is_the_099_frame_and_solves_the_fraction_in_closed_form,
    test_cd4_target_f_env_override_is_validated_and_fails_closed,
    test_mass_fix_is_off_by_default_and_leaves_the_renderer_untouched,
    test_mass_fix_makes_the_read_back_effect_exact,
    test_mass_fix_stays_inside_the_cpm_gate_and_is_idempotent,
    test_mass_fix_minnorm_is_the_smaller_correction_and_bad_modes_fail_closed,
    test_donor_weights_default_is_the_shipped_dict_and_env_overrides_it,
    test_donor_weights_reach_the_effect_and_renormalise_over_carriers,
    test_donor_weights_env_is_validated_and_fails_closed,
    test_sign_repair_is_off_by_default_and_only_touches_flipped_on_axis_entries,
    test_sign_repair_switch_is_off_by_default_and_reaches_the_effect_matrix,
    test_sign_repair_env_is_validated_and_fails_closed,
    test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte,
    test_resample_mix_over_the_whole_gated_axis_is_the_quasi_law,
    test_resample_mix_delivers_the_quasi_panel_totals_on_every_block,
    test_rescaled_control_cells_are_rejected_before_they_reach_the_emitter,
    test_emitter_flags_default_to_the_shipped_law_and_fail_closed,
    test_bundle_file_list_follows_the_manifest_and_not_a_hardcoded_abc,
    test_frozen_bundle_table_fails_closed_on_every_shape_error,
    test_repin_frozen_inputs_moves_the_table_to_a_def_bundle_and_never_builds,
    test_context_effects_default_is_none_and_only_moves_its_own_context,
    test_alpha_knee_env_is_parsed_and_fails_closed,
    test_knee_puts_the_unclamped_budget_ratio_on_one_and_is_not_a_no_op,
)

if __name__ == "__main__":
    for test in TESTS:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(TESTS)} tests passed; peak RSS={M.peak_rss_gib():.3f} GiB")


# --------------------------------------------------------------------------- 191: rc_hybrid


def test_rc_hybrid_is_off_by_default_and_the_default_path_is_byte_identical() -> None:
    """Rule 10: adding `--emitter rc_hybrid` must not move the shipped emission.

    Same form as the `--emitter` check above and for the same reason -- a same-module
    comparison cannot answer "did this patch change the DEFAULT", because the patch is on
    both sides of it. This resolves the newest committed builder that predates `rc_hybrid`
    out of git and compares sha256 of the emitted .h5ad, at the bare default AND at a15's
    shipped emission policy. It fails, never skips, if that revision cannot be found.
    """
    rev, old = _load_pre_emitter_builder(markers=("rc_hybrid",))
    assert "--emitter" in (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
    assert old.EMITTERS if hasattr(old, "EMITTERS") else True
    shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(17).normal(0.0, 0.5, (2, M.N_GENES))
        for tag, kwargs in (("bare", {}), ("shipped", shipped)):
            new = _build_artifact(M, spec, bundle, root / f"new_{tag}.h5ad",
                                  effects=effects, **kwargs)
            was = _build_artifact(old, old_spec, bundle, root / f"old_{tag}.h5ad",
                                  effects=effects, **kwargs)
            assert new["sha256"] == was["sha256"], (
                f"the patched builder changed the DEFAULT emission at {tag!r} against "
                f"revision {rev[:12]}: {new['sha256']} != {was['sha256']}"
            )
            assert new["nnz"] == was["nnz"] and new["nnz"] > 0
            explicit = _build_artifact(
                M, spec, bundle, root / f"explicit_{tag}.h5ad", effects=effects,
                emitter="quasi", **kwargs
            )
            assert explicit["sha256"] == new["sha256"]
            # ...and the resampling law still reaches the artifact, so the comparison is
            # not vacuous for the flag that was already there either.
            moved = _build_artifact(
                M, spec, bundle, root / f"rs_{tag}.h5ad", effects=effects,
                emitter="resample_mix", resample_keep_top=20, **shipped
            )
            assert moved["sha256"] != new["sha256"]


def test_rc_hybrid_delivers_the_quasi_pseudobulk_and_welds_every_cell() -> None:
    """The two identities the whole arm rests on, on a real three-context build.

    Per-gene panel totals identical to the quasi law's on every block (that is `pds` and
    `mse`'s profile half), every emitted cell at one library total (that is `nmae`), and a
    DIFFERENT file -- if the artifacts matched too, the totals identity would prove nothing.
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=9)
        effects = np.random.default_rng(23).normal(0.0, 0.5, (2, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
        base = _build_artifact(M, spec, bundle, root / "quasi.h5ad", effects=effects,
                               **shipped)
        rc = _build_artifact(M, spec, bundle, root / "rc.h5ad", effects=effects,
                             emitter="rc_hybrid", **shipped)
        assert rc["sha256"] != base["sha256"], "rc_hybrid did not reach the artifact"
        left = _panel_totals_by_block(root / "quasi.h5ad", spec)
        right = _panel_totals_by_block(root / "rc.h5ad", spec)
        assert set(left) == set(right) and len(left) == 6
        for block in left:
            assert np.array_equal(left[block], right[block]), block
        # every cell welded to one total, and the remainder cells one count above it
        import h5py

        with h5py.File(root / "rc.h5ad", "r") as handle:
            data = handle["X"]["data"][:]
            indptr = handle["X"]["indptr"][:]
        per_cell = np.add.reduceat(data.astype(np.int64), indptr[:-1])
        per_cell[indptr[1:] == indptr[:-1]] = 0
        for start in range(0, len(per_cell), spec.cells_per_pert):
            window = per_cell[start:start + spec.cells_per_pert]
            assert window.max() - window.min() <= 1, window
        assert rc["output" if "output" in rc else "nnz"] is not None
        summary = rc["rc_hybrid"]
        assert summary["blocks"] == 6
        for key in ("forced", "max_abs_resid", "moves", "unresolved"):
            assert summary[key]["max"] == 0.0, (key, summary[key])


def test_rc_hybrid_fails_closed_on_every_knob_that_would_move_its_key() -> None:
    """`rc_hybrid` keys the order on the incumbent's own gated effect. Anything that makes
    the rendered effect a different vector is an abort, not a silent re-key."""
    import pytest

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=4)
        effects = np.zeros((2, M.N_GENES))
        scales = np.full(effects.shape, 0.5)
        for emitter in M.RC_EMITTERS:
            for kwargs, needle in (
                (dict(effect_scale=0.5), "unshrunk"),
                (dict(shrink_scales=scales), "unshrunk"),
                (dict(resample_keep_top=10), "resample_mix"),
            ):
                with pytest.raises(ValueError, match=needle):
                    _build_artifact(M, spec, bundle, root / "x.h5ad", effects=effects,
                                    emitter=emitter, **kwargs)
        assert M.EMITTERS == ("quasi", "resample_mix", "rc_hybrid", "rc_hist_a",
                              "rc_callset")
        assert M.RC_EMITTERS == ("rc_hybrid", "rc_hist_a", "rc_callset")
        assert M.RC_CALL_SET_FIRST == {"rc_hybrid": True, "rc_hist_a": False,
                                       "rc_callset": True}
        assert M.build_parser().parse_args([]).emitter == "quasi"
        with pytest.raises(SystemExit):
            M.build_parser().parse_args(["--emitter", "nope"])


# ---------------------------------------------------------------------------- 201: rc_hist_a


class _PreCallSetFirstRcEmit:
    """The pre-201 `vcc2026.rc_emit`, with `call_set_first=True` accepted and dropped.

    That module has no `call_set_first` parameter -- the whole-axis re-key did not exist --
    so the shim exists only to let the CURRENT builder drive the OLD emitter on the
    `rc_hybrid` path. It refuses `call_set_first=False`, which is the case the old module
    genuinely cannot serve.
    """

    def __init__(self, module):
        self._module = module

    def __getattr__(self, name):
        return getattr(self._module, name)

    def rc_hybrid_block(self, *args, call_set_first=True, **kwargs):
        assert call_set_first is True, "the pre-201 module has no whole-axis re-key"
        return self._module.rc_hybrid_block(*args, **kwargs)


def _load_pre_201_rc_emit():
    """`src/vcc2026/rc_emit.py` as committed before 201, loaded out of git."""
    rel = "src/vcc2026/rc_emit.py"
    revs = subprocess.run(
        ["git", "-C", str(ROOT), "log", "--format=%H", "-n", "60", "--", rel],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    for rev in revs:
        source = subprocess.run(
            ["git", "-C", str(ROOT), "show", f"{rev}:{rel}"],
            capture_output=True, text=True, check=True,
        ).stdout
        if "call_set_first" in source:
            continue
        target = HERE / f".pre_201_rc_emit_{rev[:12]}.py"
        target.write_text(source)
        try:
            spec = importlib.util.spec_from_file_location(
                "test_target_pre_201_rc_emit", target
            )
            assert spec is not None and spec.loader is not None
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
        finally:
            target.unlink()
        return rev, module
    raise AssertionError(
        f"no committed revision of {rel} predates `call_set_first`; the rc_hybrid "
        "identity check cannot be made and must not be silently skipped"
    )


def test_rc_hist_a_is_off_by_default_and_the_default_path_is_byte_identical() -> None:
    """Rule 10: adding `--emitter rc_hist_a` must move NOTHING that already existed.

    Two claims, both against objects resolved out of git rather than against this same
    patch on both sides of a comparison:

      1. the DEFAULT emission is byte-identical to the newest committed builder that
         predates `rc_hist_a` -- at the bare default and at a15's shipped policy;
      2. the `rc_hybrid` path is byte-identical too, driven through the SAME pre-201
         builder AND the pre-201 `rc_emit` module, so neither half of the patch reached
         the arm 191 already built and the board already scored.

    Non-vacuity: `rc_hist_a` on the same inputs must produce a DIFFERENT artifact.
    """
    rev, old = _load_pre_emitter_builder(markers=("rc_hist_a",))
    emit_rev, old_rc_emit = _load_pre_201_rc_emit()
    assert old.EMITTERS == ("quasi", "resample_mix", "rc_hybrid")
    shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(17).normal(0.0, 0.5, (2, M.N_GENES))
        for tag, kwargs in (("bare", {}), ("shipped", shipped)):
            new = _build_artifact(M, spec, bundle, root / f"new_{tag}.h5ad",
                                  effects=effects, **kwargs)
            was = _build_artifact(old, old_spec, bundle, root / f"old_{tag}.h5ad",
                                  effects=effects, **kwargs)
            assert new["sha256"] == was["sha256"], (
                f"the patched builder changed the DEFAULT emission at {tag!r} against "
                f"revision {rev[:12]}: {new['sha256']} != {was['sha256']}"
            )
            assert new["nnz"] == was["nnz"] and new["nnz"] > 0

    # `rc_hybrid` unchanged, on the gated bundle the arm is actually exercised on. Both
    # halves of the patch are swapped out at once: the pre-201 builder AND, through the
    # shim, the pre-201 emitter module.
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=9)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(23).normal(0.0, 0.5, (2, M.N_GENES))
        shim = _PreCallSetFirstRcEmit(old_rc_emit)
        was_rc_emit = old.rc_emit
        old.rc_emit = shim
        try:
            was = _build_artifact(old, old_spec, bundle, root / "old_rc.h5ad",
                                  effects=effects, emitter="rc_hybrid", **shipped)
        finally:
            old.rc_emit = was_rc_emit
        new = _build_artifact(M, spec, bundle, root / "new_rc.h5ad", effects=effects,
                              emitter="rc_hybrid", **shipped)
        assert new["sha256"] == was["sha256"], (
            f"the 201 patch moved the `rc_hybrid` artifact against builder {rev[:12]} "
            f"and rc_emit {emit_rev[:12]}: {new['sha256']} != {was['sha256']}"
        )
        # ...and the new flag is not a no-op, or the identity above proves nothing.
        moved = _build_artifact(M, spec, bundle, root / "hist_a.h5ad", effects=effects,
                                emitter="rc_hist_a", **shipped)
        assert moved["sha256"] != new["sha256"], "rc_hist_a did not reach the artifact"
        assert moved["nnz"] > 0


def test_rc_hist_a_holds_the_totals_and_the_weld_and_re_keys_the_whole_axis() -> None:
    """`rc_hist_a` keeps every identity `rc_hybrid` has, and changes exactly one thing.

    Per-gene panel totals identical to the quasi law's on every block (`pds`, `mse`'s
    profile half), every cell welded to one library total (`nmae`), all four placement
    diagnostics zero -- and `call_set_first` recorded as 0.0 against `rc_hybrid`'s 1.0,
    which is the meta's own statement of which re-key ran.
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=9)
        effects = np.random.default_rng(23).normal(0.0, 0.5, (2, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
        base = _build_artifact(M, spec, bundle, root / "quasi.h5ad", effects=effects,
                               **shipped)
        hybrid = _build_artifact(M, spec, bundle, root / "hybrid.h5ad", effects=effects,
                                 emitter="rc_hybrid", **shipped)
        rc = _build_artifact(M, spec, bundle, root / "rc.h5ad", effects=effects,
                             emitter="rc_hist_a", **shipped)
        assert rc["sha256"] not in (base["sha256"], hybrid["sha256"])
        left = _panel_totals_by_block(root / "quasi.h5ad", spec)
        right = _panel_totals_by_block(root / "rc.h5ad", spec)
        assert set(left) == set(right) and len(left) == 6
        for block in left:
            assert np.array_equal(left[block], right[block]), block
        import h5py

        with h5py.File(root / "rc.h5ad", "r") as handle:
            data = handle["X"]["data"][:]
            indptr = handle["X"]["indptr"][:]
        per_cell = np.add.reduceat(data.astype(np.int64), indptr[:-1])
        per_cell[indptr[1:] == indptr[:-1]] = 0
        for start in range(0, len(per_cell), spec.cells_per_pert):
            window = per_cell[start:start + spec.cells_per_pert]
            assert window.max() - window.min() <= 1, window
        summary = rc["rc_hist_a"]
        assert "rc_hybrid" not in rc
        assert summary["blocks"] == 6
        for key in ("forced", "max_abs_resid", "moves", "unresolved"):
            assert summary[key]["max"] == 0.0, (key, summary[key])
        assert summary["call_set_first"]["max"] == 0.0
        assert hybrid["rc_hybrid"]["call_set_first"]["min"] == 1.0
        # the two arms carry the SAME per-block call count out of pass 1 -- the |z|
        # multiset and the BH threshold are the incumbent's under both laws
        assert (rc["rc_hist_a"]["n_pred_pass1"]
                == hybrid["rc_hybrid"]["n_pred_pass1"])


# --------------------------------------------------------------------------- 202: rc_callset


class _PreCallSetRcEmit:
    """The pre-202 `vcc2026.rc_emit`, with `call_set=None` accepted and dropped.

    That module has no `call_set` parameter -- the prescribed call set did not exist -- so
    the shim exists only to let the CURRENT builder drive the OLD emitter on the `rc_hybrid`
    and `rc_hist_a` paths. It refuses a non-None `call_set`, which is the case the old
    module genuinely cannot serve.
    """

    def __init__(self, module):
        self._module = module

    def __getattr__(self, name):
        return getattr(self._module, name)

    def rc_hybrid_block(self, *args, call_set=None, **kwargs):
        assert call_set is None, "the pre-202 module has no prescribed call set"
        return self._module.rc_hybrid_block(*args, **kwargs)


def _load_pre_202_rc_emit():
    """`src/vcc2026/rc_emit.py` as committed before 202, loaded out of git."""
    rel = "src/vcc2026/rc_emit.py"
    revs = subprocess.run(
        ["git", "-C", str(ROOT), "log", "--format=%H", "-n", "60", "--", rel],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    for rev in revs:
        source = subprocess.run(
            ["git", "-C", str(ROOT), "show", f"{rev}:{rel}"],
            capture_output=True, text=True, check=True,
        ).stdout
        if "CallSet" in source:
            continue
        if "call_set_first" not in source:
            # older still: predates 201 too, so it cannot drive `rc_hist_a` and is not the
            # comparator this test wants.
            break
        target = HERE / f".pre_202_rc_emit_{rev[:12]}.py"
        target.write_text(source)
        try:
            spec = importlib.util.spec_from_file_location(
                "test_target_pre_202_rc_emit", target
            )
            assert spec is not None and spec.loader is not None
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
        finally:
            target.unlink()
        return rev, module
    raise AssertionError(
        f"no committed revision of {rel} predates `CallSet` while still carrying "
        "`call_set_first`; the rc_hybrid / rc_hist_a identity check cannot be made and "
        "must not be silently skipped"
    )


def test_rc_callset_is_off_by_default_and_every_existing_path_is_byte_identical() -> None:
    """Rule 10: adding `--emitter rc_callset` must move NOTHING that already existed.

    Three claims, all against objects resolved out of git rather than against this same
    patch on both sides of a comparison:

      1. the DEFAULT emission is byte-identical to the newest committed builder that
         predates `rc_callset` -- at the bare default and at a15's shipped policy;
      2. the `rc_hybrid` path is byte-identical, driven through the SAME pre-202 builder
         AND the pre-202 `rc_emit` module, so neither half of the patch reached the arm
         191 built and the board scored at +0.0093325;
      3. the `rc_hist_a` path is byte-identical on the same terms.

    Non-vacuity: `rc_callset` on the same inputs must produce a DIFFERENT artifact.
    """
    rev, old = _load_pre_emitter_builder(markers=("rc_callset",))
    emit_rev, old_rc_emit = _load_pre_202_rc_emit()
    assert old.EMITTERS == ("quasi", "resample_mix", "rc_hybrid", "rc_hist_a")
    shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=False)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(17).normal(0.0, 0.5, (2, M.N_GENES))
        for tag, kwargs in (("bare", {}), ("shipped", shipped)):
            new = _build_artifact(M, spec, bundle, root / f"new_{tag}.h5ad",
                                  effects=effects, **kwargs)
            was = _build_artifact(old, old_spec, bundle, root / f"old_{tag}.h5ad",
                                  effects=effects, **kwargs)
            assert new["sha256"] == was["sha256"], (
                f"the patched builder changed the DEFAULT emission at {tag!r} against "
                f"revision {rev[:12]}: {new['sha256']} != {was['sha256']}"
            )
            assert new["nnz"] == was["nnz"] and new["nnz"] > 0

    # the two rank-controlled arms that already exist, on the gated bundle they are
    # actually exercised on. Both halves of the patch are swapped out at once.
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=9)
        old_spec = old.ChallengeSpec(
            genes=spec.genes, perts=spec.perts, contexts=list(spec.contexts),
            cells_per_pert=spec.cells_per_pert, manifest={},
        )
        effects = np.random.default_rng(23).normal(0.0, 0.5, (2, M.N_GENES))
        shim = _PreCallSetRcEmit(old_rc_emit)
        for emitter in ("rc_hybrid", "rc_hist_a"):
            was_rc_emit = old.rc_emit
            old.rc_emit = shim
            try:
                was = _build_artifact(old, old_spec, bundle, root / f"old_{emitter}.h5ad",
                                      effects=effects, emitter=emitter, **shipped)
            finally:
                old.rc_emit = was_rc_emit
            new = _build_artifact(M, spec, bundle, root / f"new_{emitter}.h5ad",
                                  effects=effects, emitter=emitter, **shipped)
            assert new["sha256"] == was["sha256"], (
                f"the 202 patch moved the `{emitter}` artifact against builder "
                f"{rev[:12]} and rc_emit {emit_rev[:12]}: "
                f"{new['sha256']} != {was['sha256']}"
            )
        # ...and the new flag is not a no-op, or the identities above prove nothing.
        head, cpm = M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM
        M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM = 8, 32
        try:
            moved = _build_artifact(M, spec, bundle, root / "callset.h5ad",
                                    effects=effects, emitter="rc_callset", **shipped)
        finally:
            M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM = head, cpm
        assert moved["sha256"] != new["sha256"], "rc_callset did not reach the artifact"
        assert moved["nnz"] > 0


def test_rc_callset_holds_the_totals_and_the_weld_and_prescribes_the_call_set() -> None:
    """`rc_callset` keeps every identity `rc_hybrid` has, and changes exactly one thing.

    Per-gene panel totals identical to the quasi law's on every block (`pds`, `mse`'s
    profile half), every cell welded to one library total (`nmae`), all four placement
    diagnostics zero -- and the block in front is the PRESCRIBED set, recorded per block
    at its prescribed size, with BH never delivering more than it (the one-sided guard:
    over-delivery means the complement was not suppressed).
    """
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        bundle, spec = _synthetic_bundle(root, M, all_gated=True, seed=9)
        effects = np.random.default_rng(23).normal(0.0, 0.5, (2, M.N_GENES))
        shipped = dict(tail_spread=1.0, keep_top=8000, keep_spread=0.25, overdispersion=4.0)
        base = _build_artifact(M, spec, bundle, root / "quasi.h5ad", effects=effects,
                               **shipped)
        hybrid = _build_artifact(M, spec, bundle, root / "hybrid.h5ad", effects=effects,
                                 emitter="rc_hybrid", **shipped)
        head, cpm = M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM
        M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM = 8, 32
        try:
            rc = _build_artifact(M, spec, bundle, root / "rc.h5ad", effects=effects,
                                 emitter="rc_callset", **shipped)
        finally:
            M.RC_CALLSET_HEAD, M.RC_CALLSET_CPM = head, cpm
        assert rc["sha256"] not in (base["sha256"], hybrid["sha256"])
        left = _panel_totals_by_block(root / "quasi.h5ad", spec)
        right = _panel_totals_by_block(root / "rc.h5ad", spec)
        assert set(left) == set(right) and len(left) == 6
        for block in left:
            assert np.array_equal(left[block], right[block]), block
        import h5py

        with h5py.File(root / "rc.h5ad", "r") as handle:
            data = handle["X"]["data"][:]
            indptr = handle["X"]["indptr"][:]
        per_cell = np.add.reduceat(data.astype(np.int64), indptr[:-1])
        per_cell[indptr[1:] == indptr[:-1]] = 0
        for start in range(0, len(per_cell), spec.cells_per_pert):
            window = per_cell[start:start + spec.cells_per_pert]
            assert window.max() - window.min() <= 1, window
        summary = rc["rc_callset"]
        assert "rc_hybrid" not in rc and "rc_hist_a" not in rc
        assert summary["blocks"] == 6
        for key in ("forced", "max_abs_resid", "moves", "unresolved"):
            assert summary[key]["max"] == 0.0, (key, summary[key])
        # the prescribed set, at its prescribed size, on every block, and the one mechanism
        # claim that is a property of THIS CODE rather than of the bundle: every gene
        # outside the prescribed set was targeted below the cap.
        assert summary["call_set_prescribed"]["min"] == 40.0
        assert summary["call_set_prescribed"]["max"] == 40.0
        assert summary["call_set_z_scale"]["max"] < 1.0, "the complement was never rescaled"
        assert summary["call_set_z_target_max"]["max"] <= M.RC_CALLSET_Z_CAP + 1e-9
        assert summary["call_set_first"]["min"] == 1.0
        # What BH then DELIVERS is a property of the bundle: a gene's reachable |z| has a
        # floor set by its own panel total against the control distribution, and this
        # bundle's floor is enormous (30 control cells, 4 emitted cells per perturbation),
        # so BH here rejects hundreds of genes no allocation can suppress. Delivery is
        # asserted where the block has real power instead --
        # `test_rc_callset_prescribes_the_union_and_suppresses_the_complement` above, and
        # on h1-test in `experiments/202-callset-rider` (median 2,635 delivered against a
        # prescribed 2,600, nothing missed on 85/85). Here it is only recorded -- on this
        # bundle BH can call anything from 0 to hundreds and both are the bundle talking.
        assert set(summary["call_set_delivered"]) == {"median", "max", "min"}
        assert summary["call_set_jaccard"]["max"] <= 1.0
        # pass 1 is untouched: the same block, the same |z| multiset, the same BH count
        assert (summary["n_pred_pass1"] == hybrid["rc_hybrid"]["n_pred_pass1"])


def test_rc_callset_prescribes_the_union_and_suppresses_the_complement() -> None:
    """The construction itself, at the emitter's own level, with real DE power.

    A synthetic control block with a sub-gate tail and enough cells that the Wilcoxon can
    resolve the BH line, driven through `rc_emit.rc_hybrid_block` directly. Five claims:

      * the prescribed set is exactly `top n_head by |H| UNION top n_cpm by control CPM`
        over the scorer's DE universe, at size `n_head + n_cpm`;
      * BH on the emitted block delivers that set -- up to the genes whose reachable |z|
        FLOOR sits above the BH line, which no allocation can suppress. This block is
        deliberately coarse (1,000 control cells against the real 18,400), so it carries
        more of those than the real surface does; `experiments/202-callset-rider` measures
        the real number on h1-test;
      * nothing prescribed is MISSED, which is the direction that would mean the
        mechanism is broken rather than the block is coarse;
      * the per-gene totals are the `rc_hybrid` run's, count for count;
      * the complement's order is `rc_hybrid`'s -- the common-factor rescaling is the
        monotone map it claims to be -- and the order inside the prescribed block is
        `|H|`.
    """
    from scipy.stats import spearmanr

    sys.path.insert(0, str(ROOT / "src"))
    from vcc2026 import rc_emit
    from vcc2026.donor_transfer import emit_counts

    rng = np.random.default_rng(202)
    n_expressed, n_low, n_controls, n_cells = 3000, 3000, 1000, 400
    n_genes = n_expressed + n_low
    mu = np.concatenate([np.exp(rng.normal(1.8, 1.3, n_expressed)),
                         rng.uniform(0.005, 0.03, n_low)])
    controls = rng.poisson(mu[None, :], size=(n_controls, n_genes)).astype(np.float64)
    control = rc_emit.control_block_from_counts(controls, min_cpm=5.0)
    keep = control.keep
    assert 2000 < keep.sum() < n_genes, "the block needs a sub-gate tail to weld against"

    profile = mu * np.exp(rng.normal(0.0, 0.10, n_genes))
    profile = profile / profile.sum()
    libraries = np.full(n_cells, 20_000, dtype=np.int64)
    key = np.where(keep, rng.normal(0.0, 1.0, n_genes), 0.0)
    spread = np.full(n_genes, 0.25)

    n_head, n_cpm = 120, 656
    call = rc_emit.CallSet(n_head, n_cpm, 1.90)
    common = dict(spread=spread, overdispersion=4.0, control=control, key=key,
                  jitter_seed=1000, emit_counts=emit_counts)
    hybrid, _ = rc_emit.rc_hybrid_block(
        profile, libraries, 0.1, np.random.default_rng(5), **common)
    block, diag = rc_emit.rc_hybrid_block(
        profile, libraries, 0.1, np.random.default_rng(5), call_set=call, **common)

    # the totals identity, against the arm this one is marginal to
    assert np.array_equal(block.sum(0), hybrid.sum(0))
    assert not np.array_equal(block, hybrid)
    assert int(block.sum(1).max() - block.sum(1).min()) <= 1

    # the prescribed set, recomputed independently of the emitter
    key_abs = np.abs(key[keep])
    head = np.zeros(int(keep.sum()), bool)
    head[np.argsort(-key_abs, kind="stable")[:n_head]] = True
    want = head.copy()
    want[np.argsort(-np.where(head, -np.inf, control.ref_lin),
                    kind="stable")[:n_cpm]] = True
    assert want.sum() == n_head + n_cpm

    _z, p, q, _lfc = rc_emit.de_from_block(block, control)
    delivered = q < 0.05
    assert diag["call_set_prescribed"] == n_head + n_cpm
    assert diag["call_set_delivered"] == int(delivered.sum())
    assert int((want & ~delivered).sum()) == 0, "a prescribed gene was not called"
    jaccard = float((delivered & want).sum() / (delivered | want).sum())
    assert abs(diag["call_set_jaccard"] - jaccard) < 1e-12
    assert jaccard > 0.80, (
        f"BH delivered {int(delivered.sum())} genes against a prescribed "
        f"{n_head + n_cpm}: jaccard {jaccard:.4f}"
    )
    assert diag["call_set_z_scale"] < 1.0, "the complement was never rescaled"

    # the order, against `rc_hybrid` on the same block rather than an absolute bar
    lp = -np.log10(np.maximum(p, 1e-300))
    _zh, ph, qh, _lh = rc_emit.de_from_block(hybrid, control)
    lph = -np.log10(np.maximum(ph, 1e-300))
    inc = qh < 0.05
    assert int(inc.sum()) > int(delivered.sum()), "the call set did not shrink"
    rho_out = spearmanr(lp[~delivered], key_abs[~delivered]).correlation
    rho_out_h = spearmanr(lph[~inc], key_abs[~inc]).correlation
    assert abs(rho_out - rho_out_h) < 0.02, (rho_out_h, rho_out)
    assert spearmanr(lp[want], key_abs[want]).correlation > 0.90
