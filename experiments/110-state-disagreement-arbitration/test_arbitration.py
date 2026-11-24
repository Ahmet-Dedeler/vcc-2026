#!/usr/bin/env python3
"""Synthetic, source-free tests for experiment-110 numerical utilities."""

from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import arbitration as ARB  # noqa: E402


def synthetic_context(name: str, scale: float = 1.0) -> ARB.DevelopmentContext:
    targets = 4
    genes = 7
    signs = np.asarray(
        [
            [1, -1, 1, -1, 1, -1, 0],
            [-1, 1, -1, 1, -1, 1, 0],
            [1, 1, -1, -1, 1, -1, 0],
            [-1, -1, 1, 1, -1, 1, 0],
        ],
        dtype=np.float64,
    )
    magnitude = np.asarray([0.45, 0.25, 0.35, 0.2, 0.3, 0.15, 0.0])
    incumbent = scale * signs * magnitude[None, :]
    state = incumbent.copy()
    state[:, :5] *= -np.asarray([1.1, 0.8, 1.25, 0.9, 1.05])[None, :]
    state[:, 5] *= -1.0  # Unsupported disagreement; must never enter training.
    bank_multiplier = np.asarray([0.94, 0.98, 1.02, 1.06])
    banks = bank_multiplier[:, None, None] * state[None, :, :]
    cross_donor = np.broadcast_to(
        np.linspace(-0.8, 0.8, genes), (targets, genes)
    ).copy()
    cross_donor += np.arange(targets)[:, None] * 0.03
    control_log_cpm = np.linspace(2.0, 8.0, genes)

    state_is_correct = np.asarray(
        [
            [True, True, False, False, True, False, False],
            [False, True, True, False, True, False, False],
            [True, False, True, False, False, True, False],
            [False, True, False, True, True, False, False],
        ]
    )
    reference = np.where(state_is_correct, state, incumbent)
    significant = np.ones_like(reference, dtype=bool)
    supported = np.asarray([True, True, True, True, True, False, False])
    forbidden = np.asarray([False, False, False, False, False, False, True])
    own_target = np.zeros_like(reference, dtype=bool)
    own_target[np.arange(targets), np.arange(targets)] = True
    return ARB.make_development_context(
        name=name,
        incumbent_lfc=incumbent,
        state_bank_lfc=banks,
        a14_cross_donor_agreement=cross_donor,
        control_log_cpm=control_log_cpm,
        reference_lfc=reference,
        reference_significant=significant,
        state_supported_gene_mask=supported,
        forbidden_gene_mask=forbidden,
        own_target_coordinate_mask=own_target,
        target_blocks_by_row=np.asarray([f"T{index}" for index in range(targets)]),
    )


def contexts() -> tuple[ARB.DevelopmentContext, ...]:
    return (
        synthetic_context("rpe1", 0.95),
        synthetic_context("jurkat", 1.00),
        synthetic_context("hepg2", 1.05),
    )


def test_fixed_grid_matches_frozen_protocol() -> None:
    protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
    grid = protocol["model"]
    assert tuple(grid["ridge_lambda"]) == ARB.RIDGE_LAMBDAS
    assert tuple(grid["probability_threshold"]) == ARB.PROBABILITY_THRESHOLDS
    assert tuple(grid["state_blend_fraction"]) == ARB.STATE_BLEND_FRACTIONS
    assert (
        tuple(grid["maximum_cosine_distance_from_a14"]) == ARB.MAXIMUM_COSINE_DISTANCES
    )
    assert len(ARB.FIXED_GRID) == 4 * 3 * 4 * 3 == 144
    assert len({tuple(config.as_array()) for config in ARB.FIXED_GRID}) == 144


def test_disagreement_mask_matches_sign_definition_and_support() -> None:
    incumbent = np.asarray([[1.0, -1.0, 0.0, -0.0]])
    state = np.asarray([[-2.0, -3.0, 4.0, 0.0]])
    mutable = np.asarray([[True, True, False, True]])
    observed = ARB.disagreement_mask(incumbent, state, mutable_mask=mutable)
    np.testing.assert_array_equal(observed, [[True, False, False, False]])


def test_permitted_feature_matrix_has_no_identity_channel() -> None:
    context = synthetic_context("rpe1")
    parameters = set(inspect.signature(ARB.build_permitted_feature_matrix).parameters)
    assert "genes" not in parameters
    assert "targets" not in parameters
    assert "context" not in parameters
    assert context.features.shape[1] == len(ARB.FEATURE_NAMES) == 9
    assert tuple(ARB.FEATURE_NAMES) == (
        "a14_abs_lfc",
        "state_abs_lfc",
        "signed_log_magnitude_ratio",
        "absolute_log_magnitude_ratio",
        "a14_cross_donor_agreement",
        "control_log_cpm",
        "target_level_a14_state_cosine",
        "state_control_bank_sign_agreement",
        "state_control_bank_variance",
    )
    assert np.isfinite(context.features).all()
    assert np.all(
        np.sign(context.incumbent[context.pair_rows, context.pair_cols])
        != np.sign(context.state[context.pair_rows, context.pair_cols])
    )


def test_target_equal_weights_assign_one_unit_per_target() -> None:
    blocks = np.asarray(["A", "A", "A", "B", "C", "C"])
    weights = ARB.target_equal_sample_weights(blocks)
    assert weights.sum() == 3.0
    for block in ("A", "B", "C"):
        np.testing.assert_allclose(weights[blocks == block].sum(), 1.0)


def test_pooled_weight_blocks_remain_target_equal_across_contexts() -> None:
    pooled_features, _, pooled_blocks = ARB._pool_contexts(contexts()[:2])
    weights = ARB.target_equal_sample_weights(pooled_blocks)
    assert len(weights) == len(pooled_features)
    for block in np.unique(pooled_blocks):
        np.testing.assert_allclose(weights[pooled_blocks == block].sum(), 1.0)


def test_standardization_is_fit_on_training_rows_only() -> None:
    train = np.arange(45, dtype=np.float64).reshape(5, 9)
    weights = np.ones(5)
    fitted = ARB.fit_standardizer(train, weights)
    original_mean = fitted.mean.copy()
    held = np.full((3, 9), 1e9)
    transformed = fitted.transform(held)
    np.testing.assert_array_equal(fitted.mean, original_mean)
    assert np.all(transformed > 1e6)


def test_weighted_ridge_logistic_uses_target_equal_weights() -> None:
    rng = np.random.default_rng(8)
    features = rng.normal(size=(60, len(ARB.FEATURE_NAMES)))
    labels = features[:, 0] + 0.5 * features[:, 1] > 0
    blocks = np.repeat(np.arange(6), 10)
    model = ARB.fit_weighted_ridge_logistic(features, labels, blocks, 0.1)
    probabilities = model.predict_probability(features)
    assert probabilities.shape == (len(features),)
    assert np.all((probabilities > 0) & (probabilities < 1))
    assert np.mean((probabilities >= 0.5) == labels) > 0.9


def test_fixed_grid_is_scored_in_every_loco_rotation() -> None:
    result = ARB.score_loco_grid(contexts())
    assert result.contexts == ARB.DEVELOPMENT_CONTEXTS
    assert result.grid == ARB.FIXED_GRID
    assert result.rotation_scores.shape == (3, 144)
    assert result.rotation_winner_indices.shape == (3,)
    assert np.isfinite(result.rotation_scores).all()


def test_selection_requires_one_identical_rotation_winner() -> None:
    source_contexts = contexts()
    scores = np.zeros((3, len(ARB.FIXED_GRID)), dtype=np.float64)
    inconsistent = ARB.GridScores(
        ARB.DEVELOPMENT_CONTEXTS,
        ARB.FIXED_GRID,
        scores,
        np.asarray([0, 1, 0]),
        np.ones(3, dtype=bool),
    )
    with mock.patch.object(ARB, "score_loco_grid", return_value=inconsistent):
        try:
            ARB.select_loco_configuration(source_contexts)
        except ARB.SelectionError:
            pass
        else:
            raise AssertionError("inconsistent LOCO winners were accepted")

    consistent = ARB.GridScores(
        ARB.DEVELOPMENT_CONTEXTS,
        ARB.FIXED_GRID,
        scores,
        np.zeros(3, dtype=np.int64),
        np.ones(3, dtype=bool),
    )
    with mock.patch.object(ARB, "score_loco_grid", return_value=consistent):
        selected = ARB.select_loco_configuration(source_contexts)
    assert selected.selected_index == 0
    assert selected.selected_config == ARB.FIXED_GRID[0]


def test_sparse_replacement_changes_only_selected_coordinates() -> None:
    incumbent = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    state = -incumbent
    mask = np.asarray([[True, False, False], [False, True, False]])
    observed = ARB.sparse_state_replacement(incumbent, state, mask, 0.5)
    expected = incumbent.copy()
    expected[mask] = 0.0
    np.testing.assert_array_equal(observed, expected)
    np.testing.assert_array_equal(observed[~mask], incumbent[~mask])


def test_row_cosine_cap_preserves_sparse_delta() -> None:
    incumbent = np.asarray([[1.0, 0.0, 0.4, 0.0]])
    candidate = np.asarray([[-1.0, 0.0, 0.4, 0.0]])
    geometry = np.asarray([True, True, True, False])
    capped, scales = ARB.apply_row_cosine_distance_cap(
        incumbent, candidate, geometry, 0.003
    )
    distance = ARB.cosine_distance_rows(incumbent, capped, geometry)[0]
    assert distance <= 0.003 + 5e-13
    assert 0 < scales[0] < 1
    np.testing.assert_array_equal(capped[:, 1:], incumbent[:, 1:])
    reconstructed = ARB.apply_row_cap_scales(incumbent, candidate, scales)
    np.testing.assert_array_equal(capped, reconstructed)


def test_coordinate_contract_is_exact_and_fails_closed() -> None:
    incumbent = np.asarray([[1.0, 2.0, 0.0], [3.0, 4.0, 0.0]])
    candidate = incumbent + 1.0
    mutable = np.asarray([[True, False, False], [True, False, False]])
    forbidden = np.asarray([False, False, True])
    observed = ARB.enforce_coordinate_contract(candidate, incumbent, mutable, forbidden)
    np.testing.assert_array_equal(observed[:, 1], incumbent[:, 1])
    np.testing.assert_array_equal(observed[:, 2], 0.0)
    bad_incumbent = incumbent.copy()
    bad_incumbent[0, 2] = 0.1
    try:
        ARB.enforce_coordinate_contract(candidate, bad_incumbent, mutable, forbidden)
    except ValueError:
        pass
    else:
        raise AssertionError("nonzero forbidden incumbent coordinate was accepted")


def test_target_block_bootstrap_is_deterministic_and_block_level() -> None:
    values = np.asarray([1.0, 3.0, 10.0, 20.0, 22.0])
    blocks = np.asarray(["A", "A", "B", "C", "C"])
    first = ARB.target_block_bootstrap(values, blocks, repetitions=250, label="test")
    second = ARB.target_block_bootstrap(values, blocks, repetitions=250, label="test")
    np.testing.assert_array_equal(first["draws"], second["draws"])
    assert first["blocks"] == 3
    np.testing.assert_allclose(first["estimate"], np.mean([2.0, 10.0, 21.0]))


def test_exactly_100_within_target_permutations_are_deterministic() -> None:
    labels = np.asarray([True, False, True, False, False, True, True, False])
    blocks = np.asarray(["A", "A", "A", "A", "B", "B", "B", "B"])
    first = list(ARB.iter_within_target_label_permutations(labels, blocks))
    second = list(ARB.iter_within_target_label_permutations(labels, blocks))
    assert len(first) == len(second) == 100
    for left, right in zip(first, second, strict=True):
        np.testing.assert_array_equal(left, right)
        for block in ("A", "B"):
            assert left[blocks == block].sum() == labels[blocks == block].sum()


def test_full_selection_null_repeats_complete_selection_100_times() -> None:
    source_contexts = contexts()
    score_calls = 0
    fit_calls = 0

    def fake_score(
        permuted_contexts: tuple[ARB.DevelopmentContext, ...],
        *,
        grid: tuple[ARB.CandidateConfig, ...] = ARB.FIXED_GRID,
    ) -> ARB.GridScores:
        nonlocal score_calls
        score_calls += 1
        signature = sum(int(context.labels.sum()) for context in permuted_contexts)
        winner = signature % len(grid)
        scores = np.zeros((3, len(grid)), dtype=np.float64)
        scores[:, winner] = 1.0
        return ARB.GridScores(
            ARB.DEVELOPMENT_CONTEXTS,
            tuple(grid),
            scores,
            np.full(3, winner, dtype=np.int64),
            np.ones(3, dtype=bool),
        )

    def fake_fit(
        features: np.ndarray,
        labels: np.ndarray,
        target_blocks: np.ndarray,
        ridge_lambda: float,
    ) -> ARB.RidgeLogisticModel:
        nonlocal fit_calls
        fit_calls += 1
        width = len(ARB.FEATURE_NAMES)
        return ARB.RidgeLogisticModel(
            ridge_lambda,
            ARB.Standardizer(np.zeros(width), np.ones(width)),
            np.zeros(width),
            0.0,
        )

    with (
        mock.patch.object(ARB, "score_loco_grid", side_effect=fake_score),
        mock.patch.object(ARB, "fit_weighted_ridge_logistic", side_effect=fake_fit),
    ):
        observed = ARB.run_full_selection_null(source_contexts)
    assert score_calls == 100
    assert fit_calls == 100
    assert observed.selected_configs.shape == (100, 4)
    assert observed.rotation_winner_indices.shape == (100, 3)
    assert observed.label_sha256.shape == (100, 3)
    assert observed.scaler_mean.shape == (100, len(ARB.FEATURE_NAMES))


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
