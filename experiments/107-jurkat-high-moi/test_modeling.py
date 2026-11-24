#!/usr/bin/env python3
"""Synthetic-only tests for the Experiment 107 estimator algebra."""

from __future__ import annotations

import dataclasses
import importlib.util
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from scipy import sparse


HERE = Path(__file__).resolve().parent


def load_modeling():
    spec = importlib.util.spec_from_file_location(
        "exp107_modeling",
        HERE / "modeling.py",
    )
    if spec is None or spec.loader is None:
        raise ImportError(HERE / "modeling.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


MODELING = load_modeling()
ASSERTIONS = unittest.TestCase()


def matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        result = np.matmul(left, right)
    assert np.all(np.isfinite(result))
    return result


def synthetic_low_rank_data(seed: int = 107):
    rng = np.random.default_rng(seed)
    channel_count = 6
    rows_per_channel = 90
    rows = channel_count * rows_per_channel
    target_count = 8
    nuisance_count = 2
    response_count = 12
    latent_rank = 3

    channels = np.repeat(np.arange(channel_count), rows_per_channel)
    probabilities = np.linspace(0.12, 0.55, target_count)
    targets = (rng.random((rows, target_count)) < probabilities[None, :]).astype(
        np.float64
    )
    nuisance = rng.normal(size=(rows, nuisance_count))
    left = rng.normal(size=(target_count, latent_rank))
    right = rng.normal(size=(latent_rank, response_count))
    coefficients = matmul(left, right)
    nuisance_coefficients = rng.normal(size=(nuisance_count, response_count))
    channel_offsets = rng.normal(size=(channel_count, response_count))
    response = (
        matmul(targets, coefficients)
        + matmul(nuisance, nuisance_coefficients)
        + channel_offsets[channels]
    )
    return response, targets, nuisance, channels, coefficients


def synthetic_production_data(seed: int = 107):
    rng = np.random.default_rng(seed)
    rows_per_channel = 16
    channels = np.repeat(np.arange(1, 17), rows_per_channel)
    target_count = 5
    targets_dense = (rng.random((channels.size, target_count)) < 0.35).astype(
        np.float64
    )
    for channel in range(1, 17):
        channel_rows = np.flatnonzero(channels == channel)
        for target in range(target_count):
            targets_dense[channel_rows[:3], target] = 1.0
    targets = sparse.csr_matrix(targets_dense)
    guide_presence = sparse.csr_matrix(np.repeat(targets_dense, 3, axis=1))
    guide_targets = np.repeat(np.arange(target_count), 3)
    nuisance = rng.normal(size=(channels.size, 3))
    coefficients = rng.normal(size=(target_count, 9))
    nuisance_coefficients = rng.normal(size=(3, 9))
    channel_offsets = rng.normal(size=(16, 9))
    response = (
        matmul(targets_dense, coefficients)
        + matmul(nuisance, nuisance_coefficients)
        + channel_offsets[channels - 1]
        + 0.01 * rng.normal(size=(channels.size, 9))
    )
    evidence = MODELING.GuideEvidence(
        guide_presence=guide_presence,
        guide_target_indices=guide_targets,
    )
    return response, targets, nuisance, channels, evidence


def evidence_for_targets(targets: sparse.csr_matrix) -> object:
    target_count = targets.shape[1]
    guide_presence = sparse.csr_matrix(np.repeat(targets.toarray(), 3, axis=1))
    return MODELING.GuideEvidence(
        guide_presence=guide_presence,
        guide_target_indices=np.repeat(np.arange(target_count), 3),
    )


def assert_same_fit(first, second) -> None:
    np.testing.assert_array_equal(
        first.diagnostics.requested_support,
        second.diagnostics.requested_support,
    )
    np.testing.assert_array_equal(
        first.diagnostics.transferable,
        second.diagnostics.transferable,
    )
    np.testing.assert_array_equal(
        first.diagnostics.fitted,
        second.diagnostics.fitted,
    )
    np.testing.assert_allclose(
        first.diagnostics.target_scales,
        second.diagnostics.target_scales,
        rtol=0.0,
        atol=0.0,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        first.diagnostics.projection_diagonal,
        second.diagnostics.projection_diagonal,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        first.diagnostics.singular_values,
        second.diagnostics.singular_values,
        rtol=0.0,
        atol=0.0,
    )
    assert first.diagnostics.effective_rank == second.diagnostics.effective_rank
    assert first.diagnostics.condition_number == second.diagnostics.condition_number
    np.testing.assert_allclose(
        first.nullspace_basis,
        second.nullspace_basis,
        rtol=0.0,
        atol=0.0,
    )
    for name in (
        "nuisance_response_slopes",
        "nuisance_target_slopes",
        "response_basis",
        "coefficients_binary_units",
        "projected_coefficients",
        "exported_effects",
    ):
        np.testing.assert_allclose(
            getattr(first, name),
            getattr(second, name),
            rtol=0.0,
            atol=0.0,
        )
    if first.identifiable_projection is None:
        assert second.identifiable_projection is None
    else:
        np.testing.assert_allclose(
            first.identifiable_projection,
            second.identifiable_projection,
            rtol=0.0,
            atol=0.0,
        )


def test_complete_channel_fold_has_no_held_parameter_or_support_leakage() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    held_channels = MODELING.FROZEN_HELD_CHANNEL_FOLDS[0]
    first = MODELING.fit_channel_fold(
        response,
        targets,
        nuisance,
        channels,
        evidence,
        held_channels=held_channels,
        ridge_lambda=3.0,
    )

    held = np.isin(channels, held_channels)
    redistributed_guides = evidence.guide_presence.toarray()
    target_values = targets.toarray()
    for row in np.flatnonzero(held):
        for target in np.flatnonzero(target_values[row]):
            redistributed_guides[row, 3 * target : 3 * target + 3] = 0.0
            redistributed_guides[row, 3 * target + row % 3] = 1.0
    guide_changed = MODELING.GuideEvidence(
        guide_presence=sparse.csr_matrix(redistributed_guides),
        guide_target_indices=evidence.guide_target_indices,
    )
    guide_only = MODELING.fit_channel_fold(
        response,
        targets,
        nuisance,
        channels,
        guide_changed,
        held_channels=held_channels,
        ridge_lambda=3.0,
    )
    assert_same_fit(first.estimator, guide_only.estimator)
    np.testing.assert_array_equal(
        first.target_support.supported,
        guide_only.target_support.supported,
    )

    changed_response = response.copy()
    changed_targets_dense = targets.toarray()
    changed_nuisance = nuisance.copy()
    changed_response[held] = 1e8 * np.arange(response.shape[1])[None, :]
    changed_targets_dense[held] = 1.0 - changed_targets_dense[held]
    changed_nuisance[held] = -1e7 + 3e6 * changed_nuisance[held]
    changed_targets = sparse.csr_matrix(changed_targets_dense)
    changed_evidence = evidence_for_targets(changed_targets)
    second = MODELING.fit_channel_fold(
        changed_response,
        changed_targets,
        changed_nuisance,
        channels,
        changed_evidence,
        held_channels=held_channels,
        ridge_lambda=3.0,
    )

    assert first.estimator.ridge_lambda == 3.0
    assert second.estimator.ridge_lambda == 3.0
    assert_same_fit(first.estimator, second.estimator)
    assert not np.array_equal(
        first.held.predicted_residual_response,
        second.held.predicted_residual_response,
    )


def test_nullspace_gauge_shift_is_removed_by_projection() -> None:
    design = np.array(
        [
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 2.0],
            [-1.0, -1.0, -2.0],
            [2.0, -1.0, 1.0],
        ]
    )
    coefficients = np.arange(15, dtype=np.float64).reshape(3, 5) / 7.0
    gauge = np.array([1.0, 1.0, -1.0])[:, None]
    response_shift = np.array([[7.0, -3.0, 2.0, 11.0, -5.0]])

    projected, diagnostics = MODELING.project_identifiable_coefficients(
        design,
        coefficients,
    )
    shifted, shifted_diagnostics = MODELING.project_identifiable_coefficients(
        design,
        coefficients + matmul(gauge, response_shift),
    )

    np.testing.assert_allclose(matmul(design, gauge), 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(projected, shifted, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        diagnostics.projection,
        shifted_diagnostics.projection,
        rtol=0.0,
        atol=0.0,
    )
    assert diagnostics.effective_rank == 2
    assert diagnostics.retained[-1] == np.bool_(False)


def test_held_channels_use_own_centering_and_training_nuisance_slopes() -> None:
    response, targets, nuisance, channels, _ = synthetic_low_rank_data()
    training = channels < 4
    model = MODELING._fit_estimator_nonproduction(
        response[training],
        targets[training],
        nuisance[training],
        channels[training],
        ridge_lambda=1.0,
        rank=3,
        support=np.ones(targets.shape[1], dtype=bool),
        projection_diagonal_threshold=0.0,
    )
    held_response = response[~training]
    held_targets = targets[~training]
    held_nuisance = nuisance[~training]
    held_channels = channels[~training]
    baseline = model.evaluate_held(
        held_response,
        held_targets,
        held_nuisance,
        held_channels,
    )

    response_shift = np.zeros_like(held_response)
    nuisance_shift = np.zeros_like(held_nuisance)
    for channel in np.unique(held_channels):
        selected = held_channels == channel
        response_shift[selected] = (channel + 1) * np.arange(
            1,
            held_response.shape[1] + 1,
        )
        nuisance_shift[selected] = (channel + 2) * np.array([13.0, -17.0])
    shifted = model.evaluate_held(
        held_response + response_shift,
        held_targets,
        held_nuisance + nuisance_shift,
        held_channels,
    )

    for channel in np.unique(held_channels):
        selected = held_channels == channel
        np.testing.assert_allclose(
            baseline.centered_response[selected].mean(axis=0),
            0.0,
            rtol=0.0,
            atol=2e-15,
        )
        np.testing.assert_allclose(
            baseline.centered_target_design[selected].mean(axis=0),
            0.0,
            rtol=0.0,
            atol=2e-15,
        )
        np.testing.assert_allclose(
            baseline.centered_nuisance[selected].mean(axis=0),
            0.0,
            rtol=0.0,
            atol=2e-15,
        )
    np.testing.assert_allclose(
        baseline.residual_response,
        shifted.residual_response,
        rtol=0.0,
        atol=3e-13,
    )
    np.testing.assert_allclose(
        baseline.predicted_residual_response,
        shifted.predicted_residual_response,
        rtol=0.0,
        atol=3e-13,
    )


def test_low_rank_fit_recovers_coefficients_in_original_binary_units() -> None:
    response, targets, nuisance, channels, truth = synthetic_low_rank_data()
    model = MODELING._fit_estimator_nonproduction(
        response,
        targets,
        nuisance,
        channels,
        ridge_lambda=0.0,
        rank=3,
        support=np.ones(targets.shape[1], dtype=bool),
        projection_diagonal_threshold=0.0,
    )

    assert model.effective_response_rank == 3
    assert model.diagnostics.effective_rank == targets.shape[1]
    assert np.ptp(model.diagnostics.target_scales) > 0.05
    np.testing.assert_allclose(
        model.coefficients_binary_units,
        truth,
        rtol=0.0,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        model.projected_coefficients,
        truth,
        rtol=0.0,
        atol=3e-12,
    )


def test_supplied_lambda_matches_the_standardized_ridge_objective() -> None:
    response, targets, nuisance, channels, _ = synthetic_low_rank_data()
    ridge_lambda = 7.0
    model = MODELING._fit_estimator_nonproduction(
        response,
        targets,
        nuisance,
        channels,
        ridge_lambda=ridge_lambda,
        rank=3,
        support=np.ones(targets.shape[1], dtype=bool),
        projection_diagonal_threshold=0.0,
    )

    centered_response = MODELING.within_channel_center(response, channels)
    centered_targets = MODELING.within_channel_center(targets, channels)
    centered_nuisance = MODELING.within_channel_center(nuisance, channels)
    response_slope = np.linalg.lstsq(
        centered_nuisance,
        centered_response,
        rcond=None,
    )[0]
    target_slope = np.linalg.lstsq(
        centered_nuisance,
        centered_targets,
        rcond=None,
    )[0]
    residual_response = centered_response - matmul(
        centered_nuisance,
        response_slope,
    )
    residual_targets = centered_targets - matmul(
        centered_nuisance,
        target_slope,
    )
    scales = np.std(residual_targets, axis=0, ddof=0)
    standardized = residual_targets / scales
    gram = matmul(standardized.T, standardized)
    gram.flat[:: gram.shape[0] + 1] += ridge_lambda
    expected = (
        np.linalg.solve(
            gram,
            matmul(standardized.T, residual_response),
        )
        / scales[:, None]
    )

    assert model.ridge_lambda == ridge_lambda
    np.testing.assert_allclose(
        model.coefficients_binary_units,
        expected,
        rtol=0.0,
        atol=4e-12,
    )


def test_sparse_nuisance_slopes_match_float64_lstsq_when_rank_deficient() -> None:
    rng = np.random.default_rng(1707)
    channels = np.repeat(np.arange(1, 5), 20)
    base = rng.normal(size=(channels.size, 2))
    nuisance = np.c_[base, base[:, 0] - 2.0 * base[:, 1]]
    centered_nuisance = MODELING.within_channel_center(nuisance, channels)
    dense_targets = (
        rng.random((channels.size, 17)) < np.linspace(0.1, 0.5, 17)
    ).astype(np.float64)
    sparse_targets = sparse.csr_matrix(dense_targets)

    _, _, pseudoinverse, observed = MODELING._sparse_target_statistics(
        sparse_targets,
        centered_nuisance,
        channels,
    )
    centered_targets = MODELING.within_channel_center(dense_targets, channels)
    expected = np.linalg.lstsq(
        centered_nuisance,
        centered_targets,
        rcond=None,
    )[0]

    np.testing.assert_allclose(
        pseudoinverse @ centered_targets,
        expected,
        rtol=2e-12,
        atol=2e-12,
    )
    np.testing.assert_allclose(observed, expected, rtol=2e-12, atol=2e-12)


def test_large_sparse_nullspace_is_deterministic_and_validated() -> None:
    rng = np.random.default_rng(2701)
    rows = 384
    columns = 257
    channels = np.repeat(np.arange(1, 17), rows // 16)
    dense = (rng.random((rows, columns)) < 0.12).astype(np.float64)
    dense[:, -1] = dense[:, 0]
    operator = MODELING._SparseResidualizedDesign(
        target_design=sparse.csr_matrix(dense),
        channels=channels,
        centered_nuisance=np.zeros((rows, 0), dtype=np.float64),
        nuisance_pseudoinverse=np.zeros((0, rows), dtype=np.float64),
        target_scales=np.ones(columns, dtype=np.float64),
    )

    first = MODELING._sparse_identifiable_nullspace(operator)
    second = MODELING._sparse_identifiable_nullspace(operator)

    np.testing.assert_allclose(first[0], second[0], rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(first[1], second[1], rtol=0.0, atol=1e-12)
    assert first[2:] == second[2:]
    assert first[0].shape[0] >= 1
    residual = operator.matmat(first[0].T)
    assert np.max(np.linalg.norm(residual, axis=0) / first[1][0]) <= (
        MODELING.IDENTIFIABLE_RELATIVE_TOLERANCE + 1e-10
    )
    np.testing.assert_allclose(
        first[0] @ first[0].T,
        np.eye(first[0].shape[0]),
        rtol=0.0,
        atol=1e-10,
    )


def test_production_projection_enforces_all_postconditions() -> None:
    rng = np.random.default_rng(3701)
    rows = 128
    channels = np.repeat(np.arange(1, 17), rows // 16)
    dense = (rng.random((rows, 4)) < 0.35).astype(np.float64)
    dense[:, 3] = dense[:, 0]
    targets = sparse.csr_matrix(dense)
    nuisance = rng.normal(size=(rows, 3))
    response = rng.normal(size=(rows, 7))
    geometry = MODELING._prepare_sparse_fit_geometry(
        response,
        targets,
        nuisance,
        channels,
        support=np.ones(targets.shape[1], dtype=bool),
        rank=4,
    )
    estimator, _ = MODELING._fit_prepared_sparse_geometry(
        geometry,
        ridge_lambda=3.0,
    )

    self_projection = estimator.project_coefficients(estimator.projected_coefficients)
    np.testing.assert_allclose(
        self_projection,
        estimator.projected_coefficients,
        rtol=0.0,
        atol=MODELING.PROJECTION_POSTCONDITION_TOLERANCE,
    )
    if estimator.nullspace_basis.size:
        np.testing.assert_allclose(
            estimator.nullspace_basis @ estimator.projected_coefficients,
            0.0,
            rtol=0.0,
            atol=MODELING.PROJECTION_POSTCONDITION_TOLERANCE,
        )
        gauge = estimator.nullspace_basis.T @ np.ones(
            (estimator.nullspace_basis.shape[0], response.shape[1])
        )
        np.testing.assert_allclose(
            estimator.project_coefficients(estimator.coefficients_binary_units + gauge),
            estimator.projected_coefficients,
            rtol=0.0,
            atol=MODELING.PROJECTION_POSTCONDITION_TOLERANCE,
        )
    assert (
        estimator.diagnostics.nullspace_annihilation_error
        <= MODELING.PROJECTION_POSTCONDITION_TOLERANCE
    )
    assert (
        estimator.diagnostics.projection_idempotence_error
        <= MODELING.PROJECTION_POSTCONDITION_TOLERANCE
    )
    assert (
        estimator.diagnostics.gauge_shift_invariance_error
        <= MODELING.PROJECTION_POSTCONDITION_TOLERANCE
    )


def test_support_diagnostics_and_exact_unsupported_fallback() -> None:
    response, targets, nuisance, channels, _ = synthetic_low_rank_data()
    # Target 2 is exactly nuisance-aliased; target 5 fails an external support
    # gate. Neither may be exported.
    nuisance_with_alias = np.c_[nuisance, targets[:, 2]]
    support = np.ones(targets.shape[1], dtype=bool)
    support[5] = False
    model = MODELING._fit_estimator_nonproduction(
        response,
        targets,
        nuisance_with_alias,
        channels,
        ridge_lambda=1.0,
        rank=3,
        support=support,
        projection_diagonal_threshold=0.0,
    )

    assert model.diagnostics.nuisance_aliased[2]
    assert not model.diagnostics.fitted[2]
    assert not model.diagnostics.fitted[5]
    assert not model.diagnostics.transferable[2]
    assert not model.diagnostics.transferable[5]
    np.testing.assert_array_equal(model.exported_effects[[2, 5]], 0.0)

    incumbent = np.arange(
        targets.shape[1] * response.shape[1],
        dtype=np.float64,
    ).reshape(targets.shape[1], response.shape[1])
    incumbent[2, 0] = np.float64(-0.0)
    source = incumbent + 1000.0
    source[2, 0] = np.float64(0.0)
    coordinate_support = np.ones(response.shape[1], dtype=bool)
    coordinate_support[-2:] = False
    candidate = MODELING.exact_unsupported_fallback(
        incumbent,
        source,
        model.diagnostics.transferable[:, None] & coordinate_support[None, :],
    )
    used = model.diagnostics.transferable[:, None] & coordinate_support[None, :]
    assert candidate[~used].tobytes() == incumbent[~used].tobytes()
    np.testing.assert_array_equal(candidate[used], source[used])


def test_relative_singular_cutoff_is_strict() -> None:
    design = np.diag([1.0, 1e-6, 1.000001e-6])
    diagnostics = MODELING.identifiable_row_space(design)
    np.testing.assert_array_equal(
        diagnostics.retained,
        [True, True, False],
    )
    assert diagnostics.singular_values[1] > 1e-6
    assert diagnostics.singular_values[2] == 1e-6


def test_frozen_support_thresholds_are_exact_and_independent() -> None:
    channels = np.repeat(np.arange(1, 5), 10)
    targets = np.zeros((channels.size, 5), dtype=np.float64)
    carrier_rows = [
        np.concatenate(
            [np.flatnonzero(channels == channel)[:5] for channel in range(1, 5)]
        ),
        np.concatenate(
            [np.flatnonzero(channels == channel)[:5] for channel in range(1, 5)]
        ),
        np.concatenate(
            [np.flatnonzero(channels == channel)[:5] for channel in range(1, 5)]
        ),
        np.concatenate(
            [
                *[np.flatnonzero(channels == channel)[:5] for channel in range(1, 4)],
                np.flatnonzero(channels == 4)[:4],
            ]
        ),
        np.concatenate(
            [
                np.flatnonzero(channels == 1)[:7],
                np.flatnonzero(channels == 2)[:7],
                np.flatnonzero(channels == 3)[:6],
            ]
        ),
    ]
    guide_distributions = (
        (10, 9, 1),
        (10, 10, 0),
        (18, 1, 1),
        (9, 9, 1),
        (10, 9, 1),
    )
    guides = np.zeros((channels.size, 15), dtype=np.float64)
    for target, (rows, distribution) in enumerate(
        zip(carrier_rows, guide_distributions, strict=True)
    ):
        targets[rows, target] = 1.0
        offset = 0
        for guide, count in enumerate(distribution):
            selected = rows[offset : offset + count]
            guides[selected, 3 * target + guide] = 1.0
            offset += count

    target_design = sparse.csr_matrix(targets)
    evidence = MODELING.GuideEvidence(
        guide_presence=sparse.csr_matrix(guides),
        guide_target_indices=np.repeat(np.arange(5), 3),
    )
    support = MODELING.derive_target_support(
        evidence,
        target_design,
        channels,
        training_channels=(1, 2, 3, 4),
    )

    np.testing.assert_array_equal(
        support.distinct_called_guides,
        [3, 2, 3, 3, 3],
    )
    np.testing.assert_array_equal(
        support.guides_with_three_untreated_cells,
        [2, 2, 1, 2, 2],
    )
    np.testing.assert_array_equal(
        support.untreated_carrier_cells,
        [20, 20, 20, 19, 20],
    )
    np.testing.assert_array_equal(
        support.training_channel_count,
        [4, 4, 4, 4, 3],
    )
    np.testing.assert_array_equal(
        support.supported,
        [True, False, False, False, False],
    )


def test_production_selector_rejects_every_frozen_contract_drift() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "lambda_grid"):
        MODELING.select_production_estimator(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            lambda_grid=(1.0, 3.0, 10.0, 30.0),
        )
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "rank"):
        MODELING.select_production_estimator(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            rank=63,
        )
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "folds"):
        MODELING.select_production_estimator(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            folds=tuple(reversed(MODELING.FROZEN_HELD_CHANNEL_FOLDS)),
        )
    drifted_channels = channels.copy()
    drifted_channels[drifted_channels == 16] = 17
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "exactly channels"):
        MODELING.fit_estimator(
            response,
            targets,
            nuisance,
            drifted_channels,
            evidence,
            ridge_lambda=3.0,
        )
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "guide_evidence"):
        MODELING.fit_estimator(
            response,
            targets,
            nuisance,
            channels,
            None,
            ridge_lambda=3.0,
        )
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "CSR"):
        MODELING.fit_estimator(
            response,
            targets.toarray(),
            nuisance,
            channels,
            evidence,
            ridge_lambda=3.0,
        )
    dense_evidence = MODELING.GuideEvidence(
        guide_presence=evidence.guide_presence.toarray(),
        guide_target_indices=evidence.guide_target_indices,
    )
    with ASSERTIONS.assertRaisesRegex(MODELING.ModelingError, "CSR"):
        MODELING.fit_estimator(
            response,
            targets,
            nuisance,
            channels,
            dense_evidence,
            ridge_lambda=3.0,
        )
    with ASSERTIONS.assertRaisesRegex(TypeError, "support"):
        MODELING.fit_channel_fold(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            held_channels=MODELING.FROZEN_HELD_CHANNEL_FOLDS[0],
            ridge_lambda=3.0,
            support=np.ones(targets.shape[1], dtype=bool),
        )


def test_largest_lambda_tie_break_is_inclusive_and_frozen() -> None:
    selected = MODELING._largest_frozen_lambda_within_relative_gain(
        np.array([1.0, 0.9995, 0.9990, 0.9989, 0.5])
    )
    assert selected == 10.0
    assert (
        MODELING._largest_frozen_lambda_within_relative_gain(
            np.array([-1.0, -1.0005, -1.001, -1.002, -2.0])
        )
        == 10.0
    )


def test_variance_weighted_gain_is_mean_over_the_fixed_held_cells() -> None:
    response = np.array([[1.0, 2.0], [-1.0, -2.0]])
    prediction = np.array([[0.5, 1.0], [-0.5, -1.0]])
    weights = MODELING._held_response_inverse_variance_weights(response)
    observed = MODELING._mean_held_cell_variance_weighted_sse_gain(
        response,
        prediction,
        weights,
    )
    expected = np.mean(
        np.sum(
            (response**2 - (response - prediction) ** 2) * weights,
            axis=1,
        )
        / np.sum(weights)
    )
    assert observed == expected


def test_sparse_production_inputs_are_never_converted_to_dense() -> None:
    class NoDenseCSR(sparse.csr_matrix):
        def __array__(self, dtype=None, copy=None):
            raise AssertionError("full sparse design was densified")

        def toarray(self, order=None, out=None):
            raise AssertionError("full sparse design was densified")

        def todense(self, order=None, out=None):
            raise AssertionError("full sparse design was densified")

    response, targets, nuisance, channels, evidence = synthetic_production_data()
    production_target_count = 18_465
    padded_targets = sparse.hstack(
        [
            targets,
            sparse.csr_matrix(
                (targets.shape[0], production_target_count - targets.shape[1])
            ),
        ],
        format="csr",
    )
    padded_guides = sparse.hstack(
        [
            evidence.guide_presence,
            sparse.csr_matrix(
                (
                    targets.shape[0],
                    3 * (production_target_count - targets.shape[1]),
                )
            ),
        ],
        format="csr",
    )
    guarded_targets = NoDenseCSR(padded_targets)
    guarded_guides = NoDenseCSR(padded_guides)
    guarded_evidence = MODELING.GuideEvidence(
        guide_presence=guarded_guides,
        guide_target_indices=np.repeat(np.arange(production_target_count), 3),
    )
    result = MODELING.fit_channel_fold(
        response,
        guarded_targets,
        nuisance,
        channels,
        guarded_evidence,
        held_channels=MODELING.FROZEN_HELD_CHANNEL_FOLDS[0],
        ridge_lambda=3.0,
    )
    assert result.estimator.exported_effects.shape == (
        production_target_count,
        response.shape[1],
    )


def test_valid_selector_uses_all_frozen_folds_and_fixed_populations() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    result = MODELING.select_production_estimator(
        response,
        targets,
        nuisance,
        channels,
        evidence,
    )
    assert tuple(result.lambda_grid) == MODELING.FROZEN_LAMBDA_GRID
    assert result.held_channel_folds == MODELING.FROZEN_HELD_CHANNEL_FOLDS
    np.testing.assert_array_equal(result.held_cell_counts, 64)
    assert result.fold_mean_gains.shape == (4, 5)
    assert result.selected_lambda in MODELING.FROZEN_LAMBDA_GRID
    assert result.estimator.requested_response_rank == 64
    assert result.fold_selected_effects.shape == (
        4,
        targets.shape[1],
        response.shape[1],
    )
    assert len(result.fold_geometry_diagnostics) == 4
    for fold_diagnostics in result.fold_solver_diagnostics:
        assert not fold_diagnostics[0].warm_started
        assert all(item.warm_started for item in fold_diagnostics[1:])
        assert all(np.all(item.converged) for item in fold_diagnostics)


def test_batched_sparse_ridge_matches_dense_solution_for_every_lambda() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    held_channels = MODELING.FROZEN_HELD_CHANNEL_FOLDS[0]
    training_rows = ~np.isin(channels, held_channels)
    training_channels = tuple(
        channel for channel in MODELING.FROZEN_CHANNELS if channel not in held_channels
    )
    support = MODELING.derive_target_support(
        evidence,
        targets,
        channels,
        training_channels=training_channels,
    )
    geometry = MODELING._prepare_sparse_fit_geometry(
        response[training_rows],
        targets[training_rows].tocsr(),
        nuisance[training_rows],
        channels[training_rows],
        support=support.supported,
        rank=MODELING.FROZEN_RESPONSE_RANK,
    )
    assert geometry.active_operator is not None
    active_count = geometry.active_indices.size
    dense_design = geometry.active_operator.matmat(
        np.eye(active_count, dtype=np.float64)
    )
    gram = dense_design.T @ dense_design
    warm_start = None
    for penalty in MODELING.FROZEN_LAMBDA_GRID:
        estimator, warm_start = MODELING._fit_prepared_sparse_geometry(
            geometry,
            ridge_lambda=penalty,
            initial_scaled_coefficients=warm_start,
        )
        expected = np.linalg.solve(
            gram + penalty * np.eye(active_count),
            geometry.ridge_rhs,
        )
        np.testing.assert_allclose(
            warm_start,
            expected,
            rtol=2e-9,
            atol=2e-10,
        )
        diagnostics = estimator.ridge_solve_diagnostics
        assert diagnostics is not None
        assert np.all(diagnostics.converged)
        assert np.max(diagnostics.relative_residual_norms) <= 1e-8


def test_batched_sparse_ridge_fails_closed_without_convergence() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    support = MODELING.derive_target_support(
        evidence,
        targets,
        channels,
        training_channels=MODELING.FROZEN_CHANNELS,
    )
    geometry = MODELING._prepare_sparse_fit_geometry(
        response,
        targets,
        nuisance,
        channels,
        support=support.supported,
        rank=MODELING.FROZEN_RESPONSE_RANK,
    )
    assert geometry.active_operator is not None
    with ASSERTIONS.assertRaisesRegex(
        MODELING.ModelingError,
        "failed strict convergence",
    ):
        MODELING._batched_spd_ridge_solve(
            geometry.active_operator,
            geometry.ridge_rhs,
            1.0,
            max_iterations=1,
        )


def test_selector_reuses_prevalidated_design_geometry() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    audit = MODELING.audit_production_design(
        targets,
        nuisance,
        channels,
        evidence,
    )
    with (
        mock.patch.object(
            MODELING,
            "_prepare_sparse_fit_geometry",
            wraps=MODELING._prepare_sparse_fit_geometry,
        ) as prepare,
        mock.patch.object(
            MODELING,
            "_sparse_identifiable_nullspace",
            wraps=MODELING._sparse_identifiable_nullspace,
        ) as nullspace,
        mock.patch.object(
            MODELING,
            "_canonical_response_basis_from_covariance",
            wraps=MODELING._canonical_response_basis_from_covariance,
        ) as response_basis,
        mock.patch.object(
            MODELING,
            "_sparse_target_statistics",
            wraps=MODELING._sparse_target_statistics,
        ) as target_statistics,
        mock.patch.object(
            MODELING,
            "_batched_spd_ridge_solve",
            wraps=MODELING._batched_spd_ridge_solve,
        ) as ridge_solve,
    ):
        MODELING.select_production_estimator(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            design_audit=audit,
        )

    assert prepare.call_count == 5
    assert nullspace.call_count == 0
    assert response_basis.call_count == 5
    assert target_statistics.call_count == 5
    assert ridge_solve.call_count == 21


def test_invalid_design_audit_fails_before_response_basis_fit() -> None:
    response, targets, nuisance, channels, evidence = synthetic_production_data()
    audit = MODELING.audit_production_design(
        targets,
        nuisance,
        channels,
        evidence,
    )
    changed_diagonal = audit.final_geometry.projection_diagonal.copy()
    changed_diagonal[0] -= 0.25
    invalid = dataclasses.replace(
        audit,
        final_geometry=dataclasses.replace(
            audit.final_geometry,
            projection_diagonal=changed_diagonal,
        ),
    )
    with mock.patch.object(
        MODELING,
        "_canonical_response_basis_from_covariance",
        wraps=MODELING._canonical_response_basis_from_covariance,
    ) as response_basis:
        with ASSERTIONS.assertRaisesRegex(
            MODELING.ModelingError,
            "projection diagonal",
        ):
            MODELING.select_production_estimator(
                response,
                targets,
                nuisance,
                channels,
                evidence,
                design_audit=invalid,
            )
    response_basis.assert_not_called()


class TestModeling(unittest.TestCase):
    test_complete_channel_fold_has_no_held_parameter_or_support_leakage = staticmethod(
        test_complete_channel_fold_has_no_held_parameter_or_support_leakage
    )
    test_nullspace_gauge_shift_is_removed_by_projection = staticmethod(
        test_nullspace_gauge_shift_is_removed_by_projection
    )
    test_held_channels_use_own_centering_and_training_nuisance_slopes = staticmethod(
        test_held_channels_use_own_centering_and_training_nuisance_slopes
    )
    test_low_rank_fit_recovers_coefficients_in_original_binary_units = staticmethod(
        test_low_rank_fit_recovers_coefficients_in_original_binary_units
    )
    test_supplied_lambda_matches_the_standardized_ridge_objective = staticmethod(
        test_supplied_lambda_matches_the_standardized_ridge_objective
    )
    test_sparse_nuisance_slopes_match_float64_lstsq_when_rank_deficient = staticmethod(
        test_sparse_nuisance_slopes_match_float64_lstsq_when_rank_deficient
    )
    test_large_sparse_nullspace_is_deterministic_and_validated = staticmethod(
        test_large_sparse_nullspace_is_deterministic_and_validated
    )
    test_production_projection_enforces_all_postconditions = staticmethod(
        test_production_projection_enforces_all_postconditions
    )
    test_support_diagnostics_and_exact_unsupported_fallback = staticmethod(
        test_support_diagnostics_and_exact_unsupported_fallback
    )
    test_relative_singular_cutoff_is_strict = staticmethod(
        test_relative_singular_cutoff_is_strict
    )
    test_frozen_support_thresholds_are_exact_and_independent = staticmethod(
        test_frozen_support_thresholds_are_exact_and_independent
    )
    test_production_selector_rejects_every_frozen_contract_drift = staticmethod(
        test_production_selector_rejects_every_frozen_contract_drift
    )
    test_largest_lambda_tie_break_is_inclusive_and_frozen = staticmethod(
        test_largest_lambda_tie_break_is_inclusive_and_frozen
    )
    test_variance_weighted_gain_is_mean_over_the_fixed_held_cells = staticmethod(
        test_variance_weighted_gain_is_mean_over_the_fixed_held_cells
    )
    test_sparse_production_inputs_are_never_converted_to_dense = staticmethod(
        test_sparse_production_inputs_are_never_converted_to_dense
    )
    test_valid_selector_uses_all_frozen_folds_and_fixed_populations = staticmethod(
        test_valid_selector_uses_all_frozen_folds_and_fixed_populations
    )
    test_batched_sparse_ridge_matches_dense_solution_for_every_lambda = staticmethod(
        test_batched_sparse_ridge_matches_dense_solution_for_every_lambda
    )
    test_batched_sparse_ridge_fails_closed_without_convergence = staticmethod(
        test_batched_sparse_ridge_fails_closed_without_convergence
    )
    test_selector_reuses_prevalidated_design_geometry = staticmethod(
        test_selector_reuses_prevalidated_design_geometry
    )
    test_invalid_design_audit_fails_before_response_basis_fit = staticmethod(
        test_invalid_design_audit_fails_before_response_basis_fit
    )


if __name__ == "__main__":
    unittest.main()
