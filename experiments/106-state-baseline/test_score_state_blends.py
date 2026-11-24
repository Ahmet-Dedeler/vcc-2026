#!/usr/bin/env python3
"""Synthetic tests for STATE/incumbent effect blending."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from cell_eval2.metrics.discrimination import discrimination_score


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "score_state_blends", HERE / "score_state_blends.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "score_state_blends.py")
SCORE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SCORE
SPEC.loader.exec_module(SCORE)


def fixtures() -> tuple[np.ndarray, ...]:
    incumbent = np.asarray([[0.4, -0.2, 0.3], [0.1, 0.2, -0.4]])
    state = np.asarray([[0.1, 0.5, 0.0], [-0.3, 0.1, 0.0]])
    bctrl = np.log1p(np.asarray([10_000.0, 20_000.0, 20_000.0]))
    ctrl_cpm = np.asarray([200_000.0, 400_000.0, 400_000.0])
    supported = np.asarray([True, True, False])
    panel = np.asarray([False, False, False])
    return incumbent, state, bctrl, ctrl_cpm, supported, panel


def test_effect_blend_preserves_unsupported_coordinates() -> None:
    incumbent, state, _, _, supported, _ = fixtures()
    mixed = SCORE.effect_blend(incumbent, state, supported, 0.5)
    np.testing.assert_array_equal(mixed[:, ~supported], incumbent[:, ~supported])


def test_direction_blend_restores_incumbent_norm() -> None:
    incumbent, state, _, _, supported, panel = fixtures()
    mixed = SCORE.direction_blend(incumbent, state, supported, ~panel, 0.35)
    for row in range(len(mixed)):
        np.testing.assert_allclose(
            np.linalg.norm(mixed[row, supported]),
            np.linalg.norm(incumbent[row, supported]),
        )


def test_composition_blend_preserves_unsupported_mass() -> None:
    incumbent, state, bctrl, _, supported, _ = fixtures()
    mixed = SCORE.composition_blend(incumbent, state, bctrl, supported, 0.5)
    q_inc = SCORE.effect_to_profile(incumbent, bctrl)
    q_mixed = SCORE.effect_to_profile(mixed, bctrl)
    np.testing.assert_allclose(q_mixed[:, ~supported], q_inc[:, ~supported])
    np.testing.assert_allclose(q_mixed[:, supported].sum(1), q_inc[:, supported].sum(1))


def test_zero_weight_is_realized_incumbent() -> None:
    incumbent, state, bctrl, ctrl_cpm, supported, panel = fixtures()
    raw = SCORE.effect_blend(incumbent, state, supported, 0.0)
    observed = SCORE.realize_effect(raw, incumbent, bctrl, ctrl_cpm, supported, panel)
    expected = SCORE.realize_effect(
        incumbent, incumbent, bctrl, ctrl_cpm, supported, panel
    )
    np.testing.assert_allclose(observed, expected)


def test_candidate_grid_preserves_unsupported_and_panel() -> None:
    incumbent, state, bctrl, ctrl_cpm, supported, panel = fixtures()
    panel[1] = True
    candidates = SCORE.candidate_effects(
        incumbent, state, bctrl, ctrl_cpm, supported, panel
    )
    identity = candidates["incumbent"]
    for candidate in candidates.values():
        np.testing.assert_allclose(candidate[:, ~supported], identity[:, ~supported])
        np.testing.assert_allclose(candidate[:, panel], identity[:, panel])
        implied = np.expm1(bctrl[None, :] + candidate).sum(axis=1)
        np.testing.assert_allclose(implied, SCORE.TS, rtol=1e-10, atol=1e-6)


def test_pds_rewards_correct_identity() -> None:
    real = np.eye(3)
    score = SCORE.pds_per_target(real, real, np.ones(3, dtype=bool))
    np.testing.assert_array_equal(score, np.ones(3))
    swapped = real[[1, 2, 0]]
    assert SCORE.pds_per_target(swapped, real, np.ones(3, dtype=bool)).mean() < 1


def test_pds_zero_vectors_use_official_midrank_convention() -> None:
    zeros = np.zeros((3, 2))
    observed = SCORE.pds_per_target(zeros, zeros, np.ones(2, dtype=bool))
    np.testing.assert_array_equal(observed, 0.5)


def test_pds_matches_cell_eval2_on_zero_and_nonzero_rows() -> None:
    predicted = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, -1.0]])
    real = np.asarray([[0.0, 0.0], [0.8, 0.2], [0.1, -0.9]])
    labels = np.asarray(["non-targeting", "A", "B", "C"])
    zero = np.zeros((1, 2))
    official = discrimination_score(
        pred_bulk=(labels, np.vstack((zero, predicted))),
        real_bulk=(labels, np.vstack((zero, real))),
        distance="cosine",
        rank_denominator="n-1",
        tie_policy="midrank",
        exclude_target_gene=False,
        control_source="pred",
    )
    expected = np.asarray([official[label] for label in labels[1:]])
    observed = SCORE.pds_per_target(predicted, real, np.ones(2, dtype=bool))
    np.testing.assert_allclose(observed, expected)


def test_common_target_axis_preserves_first_surface_order() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        first = root / "first"
        second = root / "second"
        first.mkdir()
        second.mkdir()
        for directory, targets in (
            (first, ("B", "A", "C")),
            (second, ("C", "B", "D")),
        ):
            np.savez(directory / "frozen.npz", targets=np.asarray(targets, dtype="U"))
            (directory / "manifest.json").write_text(
                json.dumps(
                    {
                        "stage": "effect-freeze",
                        "surface": directory.name,
                        "target_order_sha256": SCORE.sha256_strings(targets),
                    }
                )
            )
        output = root / "common.txt"
        expected = SCORE.sha256_strings(("B", "C"))
        SCORE.common_targets(
            (first, second),
            output,
            expected_axes=None,
            expected_count=2,
            expected_sha256=expected,
        )
        assert output.read_text().splitlines() == ["B", "C"]


def test_common_target_axis_rejects_count_or_hash_drift() -> None:
    expected = SCORE.sha256_strings(("A", "B"))
    SCORE.validate_common_target_axis(
        ("A", "B"), expected_count=2, expected_sha256=expected
    )
    for values in (("A",), ("B", "A")):
        try:
            SCORE.validate_common_target_axis(
                values, expected_count=2, expected_sha256=expected
            )
        except ValueError:
            pass
        else:
            raise AssertionError("drifted common target axis was accepted")


def test_plugin_nmse_proxy_zero_for_identity() -> None:
    real = np.asarray([[1.0, -1.0], [2.0, 3.0]])
    observed = SCORE.plugin_nmse_proxy_per_target(real, real, np.ones(2, dtype=bool))
    np.testing.assert_array_equal(observed, 0)


def test_truth_projection_uses_last_duplicate_and_measured_mask() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        (root / "genes.txt").write_text("B\nA\nB\n")
        (root / "perts.txt").write_text("T1\nT2\n")
        np.save(
            root / "delta.npy",
            np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        )
        values, measured, audit = SCORE.load_truth_effects(
            root, ("T2",), ("A", "B", "C")
        )
        np.testing.assert_array_equal(values, [[5.0, 6.0, 0.0]])
        np.testing.assert_array_equal(measured, [True, True, False])
        assert audit["duplicate_rows"] == 1
        assert audit["duplicate_policy"] == "last occurrence"


def test_bootstrap_rejects_nonfinite_values() -> None:
    try:
        SCORE.bootstrap(np.asarray([0.0, np.nan]), "nonfinite")
    except ValueError:
        pass
    else:
        raise AssertionError("non-finite bootstrap values were accepted")


def test_directional_authorization_excludes_hct116() -> None:
    provenance = {
        "checkpoint_manifest_sha256": "checkpoint-manifest",
        "checkpoint_sha256": "checkpoint",
        "panel_sha256": "panel",
        "scorer_source_sha256": "scorer",
        "cell_eval2_version": "test",
        "pairwise_full_source_sha256": "pairwise",
    }
    reports = [
        {
            "surface": "hct116",
            "target_order_sha256": "same",
            "target_allowlist_sha256": "allowlist",
            "evaluation_run_id": "run-1",
            "provenance": provenance,
            "candidates": {
                "x": {
                    "delta_vs_incumbent": {
                        metric: {"estimate": value}
                        for metric, value in {
                            "pds": 2.0,
                            "cosine": 2.0,
                            "plugin_nmse_proxy": 2.0,
                            "credit": 100.0,
                        }.items()
                    }
                }
            },
        },
        {
            "surface": "hek293t",
            "target_order_sha256": "same",
            "target_allowlist_sha256": "allowlist",
            "evaluation_run_id": "run-1",
            "provenance": provenance,
            "candidates": {
                "x": {
                    "delta_vs_incumbent": {
                        metric: {"estimate": value}
                        for metric, value in {
                            "pds": 4.0,
                            "cosine": 4.0,
                            "plugin_nmse_proxy": 4.0,
                            "credit": 7.0,
                        }.items()
                    }
                }
            },
        },
    ]
    with tempfile.TemporaryDirectory() as temp:
        dirs = []
        for report in reports:
            directory = Path(temp) / report["surface"]
            directory.mkdir()
            (directory / "blend-metrics.json").write_text(json.dumps(report))
            dirs.append(directory)
        output = Path(temp) / "summary.json"
        SCORE.summarize(dirs, output)
        result = json.loads(output.read_text())
        candidate = result["candidates"]["x"]
        assert candidate["diagnostic_macro_delta"]["credit"] == 53.5
        assert "credit" not in candidate["authorization_delta"]
        assert candidate["directional_authorization_delta"] == {
            "metric": "credit",
            "surface": "hek293t",
            "estimate": 7.0,
            "pending": "H1 confirmation",
        }
        assert result["hct116_directional_credit"] == "diagnostic-only"


def test_summary_rejects_mixed_evaluation_runs() -> None:
    provenance = {
        "checkpoint_manifest_sha256": "checkpoint-manifest",
        "checkpoint_sha256": "checkpoint",
        "panel_sha256": "panel",
        "scorer_source_sha256": "scorer",
        "cell_eval2_version": "test",
        "pairwise_full_source_sha256": "pairwise",
    }
    reports = []
    for surface, run_id in (("hct116", "run-1"), ("hek293t", "run-2")):
        reports.append(
            {
                "surface": surface,
                "target_order_sha256": "same",
                "target_allowlist_sha256": "allowlist",
                "evaluation_run_id": run_id,
                "provenance": provenance,
                "candidates": {},
            }
        )
    with tempfile.TemporaryDirectory() as temp:
        dirs = []
        for report in reports:
            directory = Path(temp) / report["surface"]
            directory.mkdir()
            (directory / "blend-metrics.json").write_text(json.dumps(report))
            dirs.append(directory)
        try:
            SCORE.summarize(dirs, Path(temp) / "summary.json")
        except ValueError:
            pass
        else:
            raise AssertionError("summary accepted reports from different runs")


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
