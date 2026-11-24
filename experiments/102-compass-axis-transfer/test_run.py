#!/usr/bin/env python3
"""Synthetic invariant tests for experiment 102.

No donor array, harness, held truth, scorer, network service, or challenge
artifact is opened.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exp102_run", HERE / "run.py")
assert spec and spec.loader
run = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = run
spec.loader.exec_module(run)


def test_coordinate_mask_and_shared_rows() -> None:
    genes = np.asarray(["safe1", "panel", "held", "safe2"])
    allowed = run.build_allowed_mask(genes, {"panel"}, {"held"})
    np.testing.assert_array_equal(allowed, [True, False, False, True])
    shared = run.common_fit_targets(
        {
            "a": ("P1", "P2", "P3", "P4", "ONLY_A"),
            "b": ("P4", "P3", "P2", "P1", "ONLY_B"),
        },
        {"P2"},
        {"P3"},
    )
    assert shared == ("P1", "P4")


def test_mask_precedes_geometry_and_projection() -> None:
    allowed = np.asarray([True, False, True, False])
    safe = np.asarray(
        [
            [3.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 4.0, 0.0],
            [2.0, 0.0, 2.0, 0.0],
        ]
    )
    contaminated = safe.copy()
    contaminated[:, 1] = [1e6, -2e6, 7e6]
    contaminated[:, 3] = [-9e7, 4e7, 3e7]
    control = np.asarray([2.0, 8e8, 3.0, -5e8])
    clean = run.fit_source_geometry_from_rows("clean", safe, control, allowed)
    dirty = run.fit_source_geometry_from_rows("dirty", contaminated, control, allowed)
    np.testing.assert_allclose(clean.mean_scaled, dirty.mean_scaled)
    np.testing.assert_allclose(clean.axis, dirty.axis)
    np.testing.assert_allclose(clean.scale, dirty.scale)
    np.testing.assert_allclose(clean.beta_mean, dirty.beta_mean)
    assert np.all(dirty.mean_scaled[~allowed] == 0)
    assert np.all(dirty.axis[~allowed] == 0)
    assert np.all(dirty.control[~allowed] == 0)


def test_streaming_geometry_matches_pure_equations() -> None:
    class FakeDonor:
        name = "fake"
        perts = ("P1", "P2", "P3")
        pert_index = {"P1": 0, "P2": 1, "P3": 2}
        delta = np.asarray(
            [
                [3.0, 1000.0, 1.0],
                [1.0, -2000.0, 4.0],
                [2.0, 7000.0, 2.0],
            ]
        )

        @staticmethod
        def projection(
            target_genes: np.ndarray, allowed_mask: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray]:
            rows = np.flatnonzero(allowed_mask)
            return rows, rows

        @staticmethod
        def project_control(
            target_genes: np.ndarray, allowed_mask: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray]:
            control = np.asarray([2.0, 8e8, 3.0])
            control[~allowed_mask] = 0
            return control, allowed_mask.copy()

    allowed = np.asarray([True, False, True])
    genes = np.asarray(["g1", "forbidden", "g2"])
    pure = run.fit_source_geometry_from_rows(
        "fake",
        FakeDonor.delta,
        np.asarray([2.0, 8e8, 3.0]),
        allowed,
    )
    streaming = run.fit_source_geometry_streaming(
        FakeDonor(),
        FakeDonor.perts,
        genes,
        allowed,
        max_rss_gib=24.0,
        block_rows=2,
    )
    np.testing.assert_allclose(streaming.mean_scaled, pure.mean_scaled)
    np.testing.assert_allclose(streaming.axis, pure.axis)
    np.testing.assert_allclose(streaming.beta_mean, pure.beta_mean)
    np.testing.assert_allclose(streaming.scale, pure.scale)
    np.testing.assert_allclose(streaming.control, pure.control)


def geometry(
    name: str,
    mean: np.ndarray,
    axis: np.ndarray,
    beta_mean: float,
    control: np.ndarray,
) -> run.SourceGeometry:
    return run.SourceGeometry(
        name=name,
        mean_scaled=np.asarray(mean, dtype=float),
        axis=np.asarray(axis, dtype=float),
        beta_mean=float(beta_mean),
        scale=1.0,
        control=np.asarray(control, dtype=float),
        gene_mask=np.ones(len(mean), dtype=bool),
        n_fit=3,
    )


def test_exact_arm_equations_and_orthogonality() -> None:
    allowed = np.asarray([True, True, True, False])
    geometries = [
        geometry("s1", [2, 0, 0, 0], [1, 0, 0, 0], 2, [3, 1, 0, 99]),
        geometry("s2", [0, 3, 0, 0], [0, 1, 0, 0], 3, [1, 4, 0, -99]),
    ]
    delta = np.asarray(
        [
            [[5, 1, 2, 1e9], [4, 2, -1, -1e9]],
            [[2, 7, 4, 5e8], [1, 8, 3, -5e8]],
        ],
        dtype=float,
    )
    carrier = np.asarray([[True, True], [True, False]])
    target_control = np.asarray([2.0, 5.0, 1.0, 8e8])
    arms, parts = run.compose_arms(delta, carrier, geometries, target_control, allowed)

    masked = delta.copy()
    masked[:, :, ~allowed] = 0
    means = np.stack([g.mean_scaled for g in geometries])
    axes = np.stack([g.axis for g in geometries])
    beta = np.einsum("ctg,cg->ct", masked, axes)
    residual = masked - beta[:, :, None] * axes[:, None, :]
    expected_a = run.carrier_mean(
        (masked - means[:, None, :]) * carrier[:, :, None], carrier
    )
    expected_b = run.carrier_mean(masked * carrier[:, :, None], carrier)
    expected_c = run.carrier_mean(residual * carrier[:, :, None], carrier)
    primary_axis = run.normalize(axes.mean(axis=0), "test primary")
    centered_beta = beta - np.asarray([2.0, 3.0])[:, None]
    expected_d = (
        expected_c + run.carrier_mean(centered_beta, carrier)[:, None] * primary_axis
    )
    expected_e = expected_c + run.carrier_mean(beta, carrier)[:, None] * primary_axis

    np.testing.assert_allclose(arms["A"], expected_a)
    np.testing.assert_allclose(arms["B"], expected_b)
    np.testing.assert_allclose(arms["C"], expected_c)
    np.testing.assert_allclose(arms["D"], expected_d)
    np.testing.assert_allclose(arms["E"], expected_e)
    assert not np.allclose(arms["D"], arms["E"])
    assert not np.allclose(arms["D"], arms["F"])
    for source in range(2):
        for target in np.flatnonzero(carrier[source]):
            assert abs(parts["source_residual"][source, target] @ axes[source]) < 1e-12
    for value in arms.values():
        assert np.all(value[:, ~allowed] == 0)
    np.testing.assert_array_equal(parts["coverage"], [2, 1])


def test_beta_permutation_preserves_marginals_and_coverage() -> None:
    centered = np.asarray([[1.0, 2.0, 3.0, 0.0], [8.0, 0.0, 5.0, 4.0]])
    carrier = np.asarray([[True, True, True, False], [True, False, True, True]])
    permuted = run.permuted_centered_beta(
        centered, carrier, surface="synthetic", draw=4
    )
    for source in range(2):
        rows = carrier[source]
        np.testing.assert_allclose(
            np.sort(permuted[source, rows]), np.sort(centered[source, rows])
        )
        assert np.all(permuted[source, ~rows] == 0)
    np.testing.assert_array_equal(permuted != 0, carrier)


def test_uncovered_target_stays_zero_in_every_arm() -> None:
    allowed = np.asarray([True, True])
    geometries = [
        geometry("s1", [1, 0], [1, 0], 1, [1, 1]),
        geometry("s2", [0, 1], [0, 1], 1, [1, 1]),
    ]
    delta = np.asarray(
        [
            [[2, 1], [0, 0]],
            [[1, 2], [0, 0]],
        ],
        dtype=float,
    )
    carrier = np.asarray([[True, False], [True, False]])
    arms, parts = run.compose_arms(
        delta,
        carrier,
        geometries,
        np.asarray([1.0, 1.0]),
        allowed,
    )
    assert parts["coverage"][1] == 0
    for value in arms.values():
        np.testing.assert_allclose(value[1], 0.0)


def test_complete_gallery_pds() -> None:
    real = np.eye(5)
    predicted = np.eye(5)
    scores = run.pds_per_target(predicted, real, np.ones(5, dtype=bool))
    np.testing.assert_allclose(scores, 1.0)


def test_primary_gate() -> None:
    deltas = {
        "h1": np.full(20, 0.004),
        "hct116": np.full(20, 0.005),
        "hek293t": np.full(20, 0.006),
    }
    nulls = {name: np.linspace(-0.001, 0.002, 20) for name in run.SURFACES}
    coverage = {name: True for name in run.SURFACES}
    result = run.evaluate_primary_gate(deltas, nulls, coverage, n_boot=500)
    assert result["decision"] == "GO"
    assert all(result["checks"].values())

    failed = {name: values.copy() for name, values in deltas.items()}
    failed["hct116"][:] = -0.0001
    result = run.evaluate_primary_gate(failed, nulls, coverage, n_boot=500)
    assert result["decision"] == "NO-GO"
    assert not result["checks"]["D_beats_A_every_surface"]


def test_memory_ceiling() -> None:
    assert run.validate_memory_limit(24.0) == 24.0
    try:
        run.validate_memory_limit(24.01)
    except ValueError:
        pass
    else:
        raise AssertionError("memory limit above 24 GiB was accepted")


def main() -> None:
    test_coordinate_mask_and_shared_rows()
    test_mask_precedes_geometry_and_projection()
    test_streaming_geometry_matches_pure_equations()
    test_exact_arm_equations_and_orthogonality()
    test_beta_permutation_preserves_marginals_and_coverage()
    test_uncovered_target_stays_zero_in_every_arm()
    test_complete_gallery_pds()
    test_primary_gate()
    test_memory_ceiling()
    print("experiment-102 synthetic invariant tests: PASS")


if __name__ == "__main__":
    main()
