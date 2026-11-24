#!/usr/bin/env python3
"""Synthetic invariants for experiment 103; no real or KOLF data are opened."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exp103_run", HERE / "run.py")
assert spec and spec.loader
run = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = run
spec.loader.exec_module(run)


def test_forbidden_coordinates_cannot_change_axis_or_projection() -> None:
    allowed = np.asarray([True, False, True, False])
    clean = np.asarray(
        [
            [3.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 4.0, 0.0],
            [2.0, 0.0, 2.0, 0.0],
        ]
    )
    contaminated = clean.copy()
    contaminated[:, 1] = [1e8, -2e8, 7e8]
    contaminated[:, 3] = [-9e9, 4e9, 3e9]
    axis_clean = run.fit_axis_from_rows(clean, allowed, label="clean")
    axis_dirty = run.fit_axis_from_rows(contaminated, allowed, label="dirty")
    np.testing.assert_allclose(axis_clean, axis_dirty)
    assert np.all(axis_dirty[~allowed] == 0)
    vector = np.asarray([5.0, 8e12, -2.0, -9e12])
    masked = vector.copy()
    masked[~allowed] = 0
    residual_clean, loading_clean = run.residualize(masked, axis_clean)
    residual_dirty, loading_dirty = run.residualize(masked, axis_dirty)
    np.testing.assert_allclose(residual_clean, residual_dirty)
    assert loading_clean == loading_dirty


def test_residual_is_orthogonal() -> None:
    axis = run.normalize(np.asarray([2.0, -1.0, 3.0]), "axis")
    residual, loading = run.residualize(
        np.asarray([5.0, 7.0, -2.0]), axis
    )
    assert np.isfinite(loading)
    assert abs(float(residual @ axis)) < 1e-12


def test_off_reconstructs_baseline_exactly() -> None:
    baseline = np.arange(20, dtype=np.float64).reshape(4, 5)
    correction = np.full_like(baseline, 0.25)
    off = run.apply_projection_correction(
        baseline, correction, enabled=False
    )
    assert np.array_equal(off, baseline)
    run.verify_exact("synthetic OFF", off, baseline)
    candidate = run.apply_projection_correction(
        baseline, correction, enabled=True
    )
    np.testing.assert_allclose(candidate, baseline - correction)

    class FakeWhitener:
        pass

    class FakeOperator:
        whitener = FakeWhitener()
        columns = np.asarray([0, 1, 2], dtype=np.int64)

    class FakeProduction:
        @staticmethod
        def apply_v9_whitening(
            effects: np.ndarray,
            whitener: FakeWhitener,
            columns: np.ndarray,
        ) -> tuple[np.ndarray, dict[str, float]]:
            assert isinstance(whitener, FakeWhitener)
            result = np.asarray(effects, dtype=np.float64).copy()
            result[:, columns] *= 2.0
            return result, {"operator_preserve_scale": 2.0}

        @staticmethod
        def apply_frozen_cd4_residual(
            whitened: np.ndarray,
            cd4: np.ndarray,
            columns: np.ndarray,
            panel_positions: np.ndarray,
        ) -> tuple[np.ndarray, dict[str, object]]:
            assert np.array_equal(columns, FakeOperator.columns)
            result = np.asarray(whitened) + np.asarray(cd4)
            result[:, panel_positions] = whitened[:, panel_positions]
            return result, {"fully_supported_targets": len(result)}

    panel_positions = np.asarray([4], dtype=np.int64)
    cd4 = np.full_like(baseline, 3.0)
    expected = baseline.copy()
    expected[:, FakeOperator.columns] *= 2.0
    expected += cd4
    expected[:, panel_positions] = baseline[:, panel_positions]
    rebuilt, _ = run.finish_prewhiten(
        FakeProduction(),
        off,
        FakeOperator(),
        panel_positions,
        cd4,
    )
    run.verify_exact("synthetic full OFF path", rebuilt, expected)


def test_sourcewise_null_preserves_marginals_and_coverage() -> None:
    loadings = np.asarray(
        [[1.0, 2.0, 3.0, 0.0], [8.0, 0.0, 5.0, 4.0]]
    )
    carrier = np.asarray(
        [[True, True, True, False], [True, False, True, True]]
    )
    permuted = run.permute_source_loadings(loadings, carrier, draw=7)
    for source in range(2):
        rows = carrier[source]
        np.testing.assert_allclose(
            np.sort(permuted[source, rows]), np.sort(loadings[source, rows])
        )
        assert np.all(permuted[source, ~rows] == 0)
    np.testing.assert_array_equal(permuted != 0, carrier)


def test_pm_survival_changes_only_retained_real_k562_share() -> None:
    class FakeModel:
        names = ["k562_gwps", "orion_hct116"]
        target_genes = np.asarray(["g1", "g2"])
        pert_index = {
            "k562_gwps": {"P": 0},
            "orion_hct116": {"P": 0},
        }
        _global = 2.0
        _renorm = {"k562_gwps": 3.0, "orion_hct116": 5.0}

    axes = {
        "k562_gwps": run.normalize(np.asarray([1.0, 1.0]), "synthetic K562"),
        "orion_hct116": np.asarray([0.0, 1.0]),
    }
    loadings = np.asarray([[4.0], [7.0]])
    carrier = np.ones((2, 1), dtype=bool)
    correction = run.projection_correction(
        FakeModel(),
        ["P"],
        axes,
        loadings,
        carrier,
        np.asarray([True]),
        np.asarray([True, False]),
        alpha=1.0,
        pm_blend=0.5,
    )
    # Equal source weights: PM-covered g1 keeps half of K562's projection;
    # PM-missing g2 keeps its full K562 projection. HCT116 is unchanged.
    k562 = 2.0 * 0.5 * 3.0 * 4.0 * axes["k562_gwps"]
    hct = 2.0 * 0.5 * 5.0 * 7.0 * axes["orion_hct116"]
    expected = k562.copy()
    expected[0] *= 0.5
    expected += hct
    np.testing.assert_allclose(correction, expected[None, :])


def test_freeze_before_truth_boundary() -> None:
    events: list[str] = []

    def freeze_builder() -> dict[str, object]:
        events.append("freeze")
        return {"sealed_before_truth": True}

    def truth_loader(_: dict[str, object]) -> str:
        events.append("truth")
        return "loaded"

    frozen, truth = run.execute_after_freeze(freeze_builder, truth_loader)
    assert frozen["sealed_before_truth"] is True
    assert truth == "loaded"
    assert events == ["freeze", "truth"]

    called = False

    def bad_truth(_: dict[str, object]) -> None:
        nonlocal called
        called = True

    try:
        run.execute_after_freeze(
            lambda: {"sealed_before_truth": False}, bad_truth
        )
    except ValueError:
        pass
    else:
        raise AssertionError("unsealed freeze reached truth")
    assert called is False


def test_no_post_truth_fitting() -> None:
    state = run.ProtocolState()
    state.require_pretruth_fit("synthetic fit")
    state.seal()
    state.open_truth()
    try:
        run.fit_axis_from_rows(
            np.eye(3),
            np.ones(3, dtype=bool),
            label="forbidden",
            state=state,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("post-truth donor-axis fit was accepted")


def test_pds_and_cluster_bootstrap() -> None:
    truth = np.eye(6)
    np.testing.assert_allclose(
        run.pds_per_target(truth, truth), np.ones(6)
    )
    deltas = np.asarray([0.01, 0.01, 0.02, 0.02, 0.03, 0.03])
    labels = np.asarray([0, 0, 1, 1, 2, 2])
    result = run.cluster_bootstrap(
        deltas, labels, repetitions=500, seed=3
    )
    assert result["valid"] is True
    assert result["lower_95"] > 0


def test_null_threshold_uses_conservative_order_statistic() -> None:
    truth = np.eye(6)
    baseline = truth.copy()
    candidate = truth.copy()
    candidate[0] = np.roll(candidate[0], 1)
    nulls = [baseline.copy() for _ in range(19)]
    extreme = baseline.copy()
    extreme[0] = np.roll(extreme[0], 1)
    nulls.append(extreme)
    held = run.HeldTruth(
        values=truth,
        panel_rows=np.arange(6),
        official_columns=np.arange(6),
        targets=tuple(f"p{i}" for i in range(6)),
        score_genes=tuple(f"g{i}" for i in range(6)),
        hashes={},
    )
    result = run.evaluate_frozen(
        baseline,
        candidate,
        nulls,
        np.arange(6),
        held,
    )
    assert result["null_percentile_method"] == "higher"
    assert result["null_p95"] == max(result["null_mean_deltas"])
    assert result["null_one_sided_p"] >= 1 / 21


def test_memory_ceiling() -> None:
    assert run.validate_memory_limit(24.0) == 24.0
    try:
        run.validate_memory_limit(24.1)
    except ValueError:
        pass
    else:
        raise AssertionError("memory limit above 24 GiB was accepted")


def main() -> None:
    test_forbidden_coordinates_cannot_change_axis_or_projection()
    test_residual_is_orthogonal()
    test_off_reconstructs_baseline_exactly()
    test_sourcewise_null_preserves_marginals_and_coverage()
    test_pm_survival_changes_only_retained_real_k562_share()
    test_freeze_before_truth_boundary()
    test_no_post_truth_fitting()
    test_pds_and_cluster_bootstrap()
    test_null_threshold_uses_conservative_order_statistic()
    test_memory_ceiling()
    print("experiment-103 synthetic invariant tests: PASS")


if __name__ == "__main__":
    main()
