#!/usr/bin/env python3
"""Strong synthetic tests for experiment 109; no real expression is opened."""

from __future__ import annotations

import dataclasses
import importlib.util
import inspect
import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise ImportError(HERE / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MODEL = load("exp109_model", "model.py")
sys.modules["model"] = MODEL
FIT = load("exp109_fit_freeze", "fit_freeze.py")
sys.modules["fit_freeze"] = FIT
SOURCE_ONLY = load("exp109_source_only", "source_only.py")
EVALUATE = load("exp109_evaluate", "evaluate.py")


def synthetic_data(seed: int = 109):
    rng = np.random.default_rng(seed)
    n_targets = 30
    guides_per_target = 3
    n_features = 8
    n_genes = 18
    latent_rank = 3
    unique_targets = tuple(f"T{i:02d}" for i in range(n_targets))
    targets = tuple(
        target for target in unique_targets for _ in range(guides_per_target)
    )
    target_effect = rng.normal(size=(n_targets, n_genes))
    features = rng.normal(size=(len(targets), n_features))
    # The final feature is exactly target-constant and therefore absent from
    # every within-target contrast despite varying across targets.
    features[:, -1] = np.repeat(rng.normal(size=n_targets), guides_per_target)
    left = rng.normal(size=(n_features - 1, latent_rank))
    right, _ = np.linalg.qr(rng.normal(size=(n_genes, latent_rank)))
    nuisance = features[:, :-1] @ left @ right.T
    effects = np.repeat(target_effect, guides_per_target, axis=0)
    effects += 0.8 * nuisance
    effects += rng.normal(scale=0.01, size=effects.shape)
    counts = rng.integers(40, 500, size=len(targets))
    identifiable = np.ones(len(targets), dtype=bool)
    gene_axis = tuple(f"G{i:02d}" for i in range(n_genes))
    panel_targets = gene_axis[:2]
    development_targets = gene_axis[2:4]
    response_genes = gene_axis[4:]
    response_axis = MODEL.authenticate_response_axis(
        gene_axis,
        panel_targets,
        development_targets,
        gene_axis_contract=MODEL.AxisContract.from_values(gene_axis),
        panel_targets_contract=MODEL.AxisContract.from_values(panel_targets),
        development_targets_contract=MODEL.AxisContract.from_values(
            development_targets
        ),
        response_axis_contract=MODEL.AxisContract.from_values(response_genes),
    )
    return (
        effects,
        features,
        targets,
        counts,
        identifiable,
        response_axis,
        unique_targets,
    )


def fitted_model(
    *,
    eta: float = 0.75,
    seed: int = 109,
) -> tuple[object, tuple[object, ...]]:
    data = synthetic_data(seed)
    effects, features, targets, _, identifiable, response_axis, unique = data
    candidate = MODEL.Candidate(
        response_rank=3,
        ridge=0.1,
        eta=eta,
        weight_power=0.0,
        huber_c=2.0,
    )
    fitted = MODEL.fit_stable_subspace_model(
        effects,
        features,
        targets,
        identifiable,
        unique[:18],
        response_axis,
        candidate,
        seed=seed,
    )
    return fitted, data


def test_preregistered_grid_and_contract_are_frozen() -> None:
    protocol = FIT.validate_preregistration()
    assert (
        protocol["status"]
        == "source_only_runner_ready_after_prefit_extraction_repair"
    )
    assert protocol["amendments"][-1]["source_cache_published_before_repair"] is False
    assert protocol["amendments"][-1]["model_fit_before_repair"] is False
    assert protocol["amendments"][-1]["held_outcomes_opened"] is False
    assert (
        protocol["reuse"]["experiment_083_extractor_pooling_contract"][
            "maximum_absolute_error"
        ]
        == 2e-5
    )
    assert len(MODEL.candidate_grid()) == 49
    assert MODEL.candidate_grid()[0].is_identity
    assert {candidate.weight_power for candidate in MODEL.candidate_grid()} == {0.0}
    assert "hct116_directional_credit_delta_minimum" not in protocol["promotion_gates"]
    assert "h1_directional_credit_delta_minimum" in protocol["promotion_gates"]
    assert (
        protocol["production_adapter"]["zero_mean_stage"]
        == "raw K562 singleton response correction before PerturbMatch target-by-gene retention"
    )
    assert "feature_nullspace_prediction_change_abs_max" in protocol["promotion_gates"]
    assert protocol["prohibited"][-1] == "push"


def test_gauge_shift_leaves_deployed_corrections_unchanged() -> None:
    fitted, data = fitted_model()
    _, features, _, _, _, _, _ = data
    target_features = features[::3]
    support = np.asarray([index % 4 != 0 for index in range(len(target_features))])
    nuisance = fitted.nuisance(target_features)
    shift = 1e6 + np.linspace(-3.0, 2.0, nuisance.shape[1])
    target_ids = tuple(f"T{i:02d}" for i in range(len(target_features)))
    supported_ids = tuple(
        target for index, target in enumerate(target_ids) if index % 4 != 0
    )
    contract = MODEL.AxisContract.from_values(supported_ids)
    first = MODEL.centered_singleton_corrections(
        nuisance,
        target_ids,
        support,
        eta=fitted.candidate.eta,
        centering_axis_contract=contract,
    )
    shifted = MODEL.centered_singleton_corrections(
        nuisance + shift[None, :],
        target_ids,
        support,
        eta=fitted.candidate.eta,
        centering_axis_contract=contract,
    )
    tolerance = MODEL.gauge_invariance_tolerance(nuisance, nuisance + shift)
    np.testing.assert_allclose(first, shifted, rtol=0.0, atol=tolerance)


def test_full_scale_centering_meets_frozen_numerical_gates() -> None:
    rng = np.random.default_rng(MODEL.SPLIT_SEED)
    nuisance = rng.normal(size=(6689, 188))
    support = np.ones(len(nuisance), dtype=bool)
    target_ids = tuple(f"S{i:04d}" for i in range(len(nuisance)))
    contract = MODEL.AxisContract.from_values(target_ids)
    shift = 1e6 + rng.normal(size=nuisance.shape[1])
    correction = MODEL.centered_singleton_corrections(
        nuisance,
        target_ids,
        support,
        eta=0.75,
        centering_axis_contract=contract,
    )
    shifted_input = nuisance + shift
    shifted = MODEL.centered_singleton_corrections(
        shifted_input,
        target_ids,
        support,
        eta=0.75,
        centering_axis_contract=contract,
    )
    assert float(np.max(np.abs(correction.sum(axis=0)))) <= 1e-12
    tolerance = MODEL.gauge_invariance_tolerance(nuisance, shifted_input)
    assert float(np.max(np.abs(correction - shifted))) <= tolerance


def test_corrections_are_zero_mean_and_exactly_zero_off_support() -> None:
    fitted, data = fitted_model()
    _, features, _, _, _, _, _ = data
    target_features = features[::3]
    support = np.asarray([True] * 17 + [False] * 13)
    target_ids = tuple(f"T{i:02d}" for i in range(len(target_features)))
    contract = MODEL.AxisContract.from_values(target_ids[:17])
    corrections = MODEL.deploy_singleton_corrections(
        fitted,
        target_features,
        target_ids,
        support,
        centering_axis_contract=contract,
    )
    np.testing.assert_array_equal(
        corrections[~support], np.zeros_like(corrections[~support])
    )
    np.testing.assert_allclose(
        corrections[support].sum(axis=0), 0.0, rtol=0.0, atol=1e-14
    )


def test_model_uses_only_stable_feature_contrast_subspace() -> None:
    fitted, data = fitted_model()
    _, features, _, _, _, _, _ = data
    assert fitted.stable_rank < int(fitted.feature_keep.sum())
    kept = np.flatnonzero(fitted.feature_keep)
    standardized = fitted.standardized_features(features[:5])
    projector = fitted.stable_feature_basis @ fitted.stable_feature_basis.T
    residual_projector = np.eye(len(kept)) - projector
    direction = residual_projector[:, -1]
    direction /= np.linalg.norm(direction)
    changed = features[:5].copy()
    changed[:, kept] = fitted.feature_mean + fitted.feature_scale * (
        standardized + 100.0 * direction
    )
    np.testing.assert_allclose(
        fitted.nuisance(features[:5]),
        fitted.nuisance(changed),
        rtol=0.0,
        atol=2e-12,
    )
    projected = fitted.projected_standardized_features(features[:5])
    residual = standardized - projected
    np.testing.assert_allclose(
        residual @ fitted.stable_feature_basis,
        0.0,
        rtol=0.0,
        atol=2e-12,
    )


def test_extrapolation_diagnostics_detect_score_space_outliers() -> None:
    fitted, data = fitted_model()
    _, features, _, _, _, _, _ = data
    reference = features[:60]
    query = features[60:].copy()
    query[0] += 1000.0
    diagnostics = MODEL.stable_score_extrapolation_diagnostics(
        fitted,
        reference,
        query,
    )
    assert diagnostics["query_rows"] == len(query)
    assert diagnostics["query_distance_max"] > diagnostics["reference_distance_p99"]
    assert diagnostics["query_fraction_above_reference_p99"] > 0.0


def test_fit_is_isolated_from_tune_and_test_rows() -> None:
    fitted, data = fitted_model(seed=110)
    effects, features, targets, _, identifiable, response_axis, unique = data
    mutated_effects = effects.copy()
    mutated_features = features.copy()
    nonfit = np.asarray([target not in set(unique[:18]) for target in targets])
    mutated_effects[nonfit] = 1e9
    mutated_features[nonfit] = -1e8
    candidate = fitted.candidate
    refit = MODEL.fit_stable_subspace_model(
        mutated_effects,
        mutated_features,
        targets,
        identifiable,
        unique[:18],
        response_axis,
        candidate,
        seed=110,
    )
    for field in (
        "feature_mean",
        "feature_scale",
        "feature_keep",
        "stable_feature_basis",
        "stable_singular_values",
        "response_mask",
        "response_basis",
        "coefficients",
    ):
        np.testing.assert_array_equal(getattr(fitted, field), getattr(refit, field))


def test_eta_zero_is_exact_identity_everywhere() -> None:
    fitted, data = fitted_model(eta=0.0)
    effects, features, targets, counts, identifiable, _, unique = data
    support = np.ones(len(features[::3]), dtype=bool)
    target_ids = tuple(f"T{i:02d}" for i in range(len(features[::3])))
    contract = MODEL.AxisContract.from_values(target_ids)
    corrections = MODEL.deploy_singleton_corrections(
        fitted,
        features[::3],
        target_ids,
        support,
        centering_axis_contract=contract,
    )
    np.testing.assert_array_equal(corrections, np.zeros_like(corrections))
    singleton_mean = MODEL.singleton_nuisance_mean(
        fitted,
        features[::3],
        target_ids,
        support,
        centering_axis_contract=contract,
    )
    truth, baseline, candidate, labels, _ = MODEL.held_guide_predictions(
        fitted,
        effects,
        features,
        targets,
        counts,
        identifiable,
        unique[18:24],
        singleton_mean,
    )
    assert len(labels) > 0 and len(truth) > 0
    np.testing.assert_array_equal(candidate, baseline)


def test_split_is_target_level_disjoint_and_panel_safe() -> None:
    panel = ("P0", "P1")
    targets = [f"T{i:03d}" for i in range(30)] + list(panel)
    splits = FIT.development_split(targets, panel)
    MODEL.assert_split_contract(splits, panel)
    assert not (set().union(*map(set, splits.values())) & set(panel))
    assert sum(map(len, splits.values())) == 30


def test_source_null_assignments_are_outcome_free_split_local_derangements() -> None:
    signature = inspect.signature(FIT.build_split_local_cross_target_strata)
    assert "guide_effects" not in signature.parameters
    assert "cell_counts" not in signature.parameters
    targets_by_split = {
        "fit": ("A", "B", "C", "D"),
        "tune": ("E", "F", "G", "H"),
        "test": ("I", "J", "K", "L"),
    }
    guide_targets = tuple(
        target
        for split_name in FIT.REQUIRED_SPLITS
        for target in targets_by_split[split_name]
        for _ in range(2)
    )
    features = np.arange(len(guide_targets) * 3, dtype=np.float64).reshape(
        len(guide_targets),
        3,
    )
    identifiable = np.ones(len(guide_targets), dtype=bool)
    strata, audit = FIT.build_split_local_cross_target_strata(
        features,
        guide_targets,
        identifiable,
        targets_by_split,
    )
    assert audit["outcome_values_used"] is False
    assert audit["expression_values_used"] is False
    assert audit["split_local"] is True
    assert {
        details["coarsening_level"] for details in audit["splits"].values()
    } == {"guide_count"}
    labels = np.asarray(guide_targets)
    split_for_target = {
        target: split_name
        for split_name, targets in targets_by_split.items()
        for target in targets
    }
    for replicate in range(5):
        seed = MODEL.EXP083.stable_seed(
            MODEL.SPLIT_SEED,
            "exp109-cross_target-null",
            replicate,
        )
        permutation = FIT.cross_target_permutation_indices(
            guide_targets,
            strata,
            seed=seed,
            row_mask=identifiable,
        )
        for row, donor_row in enumerate(permutation):
            if guide_targets[row] == "I":
                assert labels[donor_row] != "I"
            if guide_targets[row] in split_for_target:
                assert split_for_target[labels[donor_row]] == split_for_target[
                    guide_targets[row]
                ]
        within = FIT.within_target_permutation_indices(
            guide_targets,
            seed=seed,
            row_mask=identifiable,
            targets=targets_by_split["test"],
        )
        test_rows = np.asarray(
            [target in set(targets_by_split["test"]) for target in guide_targets]
        )
        assert np.all(within[test_rows] != np.flatnonzero(test_rows))


def test_panel_rows_require_canonical_order_or_explicit_remap() -> None:
    canonical = ("A", "B", "C")
    reversed_targets = tuple(reversed(canonical))
    values = np.asarray([[3.0], [2.0], [1.0]])
    try:
        EVALUATE.canonicalize_panel_rows(
            values,
            reversed_targets,
            canonical,
        )
    except ValueError as exc:
        assert "explicit remapping" in str(exc)
    else:
        raise AssertionError("reversed panel rows were accepted implicitly")
    remapped = EVALUATE.canonicalize_panel_rows(
        values,
        reversed_targets,
        canonical,
        allow_remap=True,
    )
    np.testing.assert_array_equal(remapped[:, 0], [1.0, 2.0, 3.0])


def test_excluded_coordinates_and_pm_retention_are_exact() -> None:
    fitted, _ = fitted_model()
    response = np.ones((2, fitted.response_width), dtype=np.float64)
    full = MODEL.expand_response_corrections(fitted, response)
    np.testing.assert_array_equal(
        full[:, ~fitted.response_mask],
        np.zeros((2, int(np.sum(~fitted.response_mask)))),
    )
    pm_support = np.zeros(full.shape[1], dtype=bool)
    pm_support[::2] = True
    pm_target_support = np.asarray([True, False])
    retained = MODEL.apply_k562_slot_retention(
        full,
        pm_target_support,
        pm_support,
    )
    active = fitted.response_mask
    np.testing.assert_array_equal(retained[:, ~active], 0.0)
    np.testing.assert_array_equal(retained[0, active & pm_support], 0.5)
    np.testing.assert_array_equal(retained[0, active & ~pm_support], 1.0)
    np.testing.assert_array_equal(retained[1, active], 1.0)


def test_excluded_response_coordinates_cannot_affect_fit() -> None:
    fitted, data = fitted_model(seed=111)
    effects, features, targets, _, identifiable, response_axis, unique = data
    mutated = effects.copy()
    mutated[:, ~response_axis.response_mask] = 1e12
    refit = MODEL.fit_stable_subspace_model(
        mutated,
        features,
        targets,
        identifiable,
        unique[:18],
        response_axis,
        fitted.candidate,
        seed=111,
    )
    for field in (
        "feature_mean",
        "feature_scale",
        "feature_keep",
        "stable_feature_basis",
        "stable_singular_values",
        "response_mask",
        "response_basis",
        "coefficients",
    ):
        np.testing.assert_array_equal(getattr(fitted, field), getattr(refit, field))
    forged = dataclasses.replace(
        response_axis,
        response_mask=np.ones_like(response_axis.response_mask),
    )
    try:
        MODEL.fit_stable_subspace_model(
            effects,
            features,
            targets,
            identifiable,
            unique[:18],
            forged,
            fitted.candidate,
            seed=111,
        )
    except ValueError as exc:
        assert "caller-supplied response mask" in str(exc)
    else:
        raise AssertionError("forged response mask was accepted")


def test_centering_axis_rejects_subset_duplicates_and_reordering() -> None:
    nuisance = np.arange(20, dtype=np.float64).reshape(5, 4)
    targets = ("A", "B", "C", "D", "E")
    support = np.ones(5, dtype=bool)
    contract = MODEL.AxisContract.from_values(targets)
    MODEL.centered_singleton_corrections(
        nuisance,
        targets,
        support,
        eta=0.5,
        centering_axis_contract=contract,
    )
    for changed_targets, changed_support in (
        (targets[:-1], support[:-1]),
        (("A", "B", "B", "D", "E"), support),
        (tuple(reversed(targets)), support),
    ):
        try:
            MODEL.centered_singleton_corrections(
                nuisance[: len(changed_targets)],
                changed_targets,
                changed_support,
                eta=0.5,
                centering_axis_contract=contract,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("invalid centering population was accepted")


def test_raw_reconstruction_gain_does_not_substitute_for_residual_gate() -> None:
    effects = np.asarray(
        [
            [0.9574436971, 2.5742170264, -1.4890638997, -1.5116905383],
            [-0.5037588707, -0.7891256446, -0.2406150545, 0.1388340711],
            [-0.4283496749, -0.2584153502, 0.1015522669, 0.5786678080],
        ]
    )
    nuisance = np.asarray(
        [
            [0.2596235228, -0.4692765083, -0.6156055525, -1.3556699608],
            [0.2169340870, 0.8300470096, 1.4350183979, -1.0245757361],
            [-0.3220521930, 0.5844500736, -0.4938091640, 0.2931887607],
        ]
    )
    singleton_mean = np.asarray(
        [-0.2433260452, -0.3149217091, 1.1890383371, 0.0833709347]
    )
    fitted = MODEL.StableSubspaceModel(
        candidate=MODEL.Candidate(4, 0.0, 1.0, 0.0),
        feature_mean=np.zeros(4),
        feature_scale=np.ones(4),
        feature_keep=np.ones(4, dtype=bool),
        stable_feature_basis=np.eye(4),
        stable_singular_values=np.ones(4),
        response_mask=np.ones(4, dtype=bool),
        gene_axis_sha256="0" * 64,
        panel_targets_sha256="1" * 64,
        development_targets_sha256="2" * 64,
        response_axis_sha256="3" * 64,
        fit_targets_sha256="4" * 64,
        response_basis=np.eye(4),
        coefficients=np.eye(4),
    )
    args = (
        fitted,
        effects,
        nuisance,
        ("T", "T", "T"),
        np.ones(3),
        np.ones(3, dtype=bool),
        ("T",),
    )
    raw = MODEL.held_guide_reconstruction_diagnostic(*args)
    deployed = MODEL.held_guide_predictions(*args, singleton_mean)
    raw_metrics = MODEL.target_macro_metrics(*raw)
    deployed_metrics = MODEL.target_macro_metrics(*deployed[:4])
    assert raw_metrics["cosine_delta"] > 0.3
    assert not FIT.candidate_passes_tune(deployed_metrics)


def test_target_constant_nuisance_receives_zero_source_credit() -> None:
    rng = np.random.default_rng(109)
    effects = rng.normal(size=(3, 4))
    nuisance = np.repeat(
        np.asarray([[2.0, -1.0, 0.5, 3.0]]),
        len(effects),
        axis=0,
    )
    fitted = MODEL.StableSubspaceModel(
        candidate=MODEL.Candidate(4, 0.0, 1.0, 0.0),
        feature_mean=np.zeros(4),
        feature_scale=np.ones(4),
        feature_keep=np.ones(4, dtype=bool),
        stable_feature_basis=np.eye(4),
        stable_singular_values=np.ones(4),
        response_mask=np.ones(4, dtype=bool),
        gene_axis_sha256="0" * 64,
        panel_targets_sha256="1" * 64,
        development_targets_sha256="2" * 64,
        response_axis_sha256="3" * 64,
        fit_targets_sha256="4" * 64,
        response_basis=np.eye(4),
        coefficients=np.eye(4),
    )
    result = MODEL.held_guide_predictions(
        fitted,
        effects,
        nuisance,
        ("T", "T", "T"),
        np.ones(3),
        np.ones(3, dtype=bool),
        ("T",),
        nuisance[0],
    )
    metrics = MODEL.target_macro_metrics(*result[:4])
    assert metrics["cosine_delta"] == 0.0
    assert metrics["nmse_delta"] == 0.0


def test_complete_selector_bootstrap_and_null_gates() -> None:
    records = []
    winner = MODEL.Candidate(8, 0.1, 0.25, 0.0)
    for candidate in MODEL.candidate_grid():
        metrics = {
            "targets": 10,
            "guide_holdouts": 30,
            "cosine_delta": 0.02 if candidate == winner else -0.01,
            "nmse_delta": -0.01 if candidate == winner else 0.01,
        }
        records.append((candidate, metrics))
    selected, metrics = FIT.select_candidate(records)
    assert selected == winner
    assert metrics["cosine_delta"] == 0.02
    bootstrap = FIT.target_cluster_bootstrap(
        np.asarray([0.01, 0.03, 0.02, 0.04]),
        ("A", "A", "B", "B"),
        repetitions=1000,
        seed=109,
    )
    decision = FIT.source_test_decision(
        {
            "targets": 102,
            "guide_holdouts": 306,
            "cosine_delta": 0.025,
            "nmse_delta": -0.01,
        },
        bootstrap_samples=bootstrap,
        within_target_null=np.asarray([-0.02, -0.01, 0.0]),
        cross_target_null=np.asarray([-0.03, -0.02, -0.01]),
        direction_preserving_shrink_delta=0.0,
    )
    assert decision["passes"] is True
    FIT.validate_source_test(decision)
    forged = json.loads(json.dumps(decision))
    forged["bootstrap_lower_95"] = -1.0
    try:
        FIT.validate_source_test(forged)
    except ValueError as exc:
        assert "bootstrap gate" in str(exc)
    else:
        raise AssertionError("self-inconsistent source gate was accepted")


def test_large_validator_allclose_is_numpy_version_independent() -> None:
    shape = (300, 2_000)
    actual = np.zeros(shape, dtype=np.float64)
    expected = np.zeros(shape, dtype=np.float64)
    FIT.assert_finite_allclose(
        actual,
        expected,
        label="large synthetic correction",
        rtol=0.0,
        atol=1e-12,
    )
    actual[-1, -1] = 1e-6
    try:
        FIT.assert_finite_allclose(
            actual,
            expected,
            label="large synthetic correction",
            rtol=0.0,
            atol=1e-12,
        )
    except ValueError as exc:
        assert "differs at flat index" in str(exc)
    else:
        raise AssertionError("large-array corruption passed the validator")


def test_frozen_model_round_trip() -> None:
    fitted, _ = fitted_model()
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "model.npz"
        FIT.save_frozen_model(path, fitted)
        restored = FIT.load_frozen_model(path)
    assert restored.candidate == fitted.candidate
    for field in (
        "feature_mean",
        "feature_scale",
        "feature_keep",
        "stable_feature_basis",
        "stable_singular_values",
        "response_mask",
        "gene_axis_sha256",
        "panel_targets_sha256",
        "development_targets_sha256",
        "response_axis_sha256",
        "fit_targets_sha256",
        "response_basis",
        "coefficients",
    ):
        np.testing.assert_array_equal(getattr(restored, field), getattr(fitted, field))


def _full_freeze_fixture(
    directory: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    source_axes = FIT.load_canonical_source_axes()
    genes, gene_contract = FIT.canonical_axis("gene_axis", source_axes)
    panel, panel_contract = FIT.canonical_axis("panel_targets", source_axes)
    development, development_contract = FIT.canonical_axis(
        "development_targets", source_axes
    )
    response, response_contract = FIT.canonical_axis("response_axis", source_axes)
    response_axis = MODEL.authenticate_response_axis(
        genes,
        panel,
        development,
        gene_axis_contract=gene_contract,
        panel_targets_contract=panel_contract,
        development_targets_contract=development_contract,
        response_axis_contract=response_contract,
    )
    centering, centering_contract = FIT.canonical_axis(
        "singleton_centering_targets", source_axes
    )
    FIT.write_axis_payload(
        directory / "axes.json",
        response_axis,
        centering,
        centering_contract,
    )
    splits = FIT.development_split(development, panel)
    winner = MODEL.Candidate(8, 0.1, 0.25, 0.0)
    fit_tune = tuple(sorted((*splits["fit"], *splits["tune"])))
    feature_keep = np.arange(229) < 225
    common = {
        "candidate": winner,
        "feature_mean": np.zeros(225),
        "feature_scale": np.ones(225),
        "feature_keep": feature_keep,
        "stable_feature_basis": np.eye(225, 188),
        "stable_singular_values": np.ones(188),
        "response_mask": response_axis.response_mask,
        "gene_axis_sha256": response_axis.gene_axis_contract.sha256,
        "panel_targets_sha256": response_axis.panel_targets_contract.sha256,
        "development_targets_sha256": (
            response_axis.development_targets_contract.sha256
        ),
        "response_axis_sha256": response_axis.response_axis_contract.sha256,
        "response_basis": np.eye(8, len(response)),
        "coefficients": np.zeros((188, 8)),
    }
    source_test_model = MODEL.StableSubspaceModel(
        **common,
        fit_targets_sha256=FIT.sha256_lines(fit_tune),
    )
    final_model = MODEL.StableSubspaceModel(
        **common,
        fit_targets_sha256=development_contract.sha256,
    )
    FIT.save_frozen_model(
        directory / "source-test-model.npz",
        source_test_model,
        role="source_test_model",
    )
    FIT.save_frozen_model(directory / "selected-model.npz", final_model)
    selection_records = []
    for candidate in MODEL.candidate_grid():
        selection_records.append(
            {
                "candidate": FIT.candidate_payload(candidate),
                "metrics": {
                    "targets": 101,
                    "guide_holdouts": 303,
                    "cosine_delta": 0.02 if candidate == winner else -0.01,
                    "nmse_delta": -0.01 if candidate == winner else 0.01,
                },
            }
        )
    (directory / "tune-selection.json").write_text(
        json.dumps(
            {
                "candidate_records": selection_records,
                "selected_candidate": FIT.candidate_payload(winner),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    singleton_features = np.zeros((len(centering), 229), dtype=np.float64)
    development_ids = np.asarray([target for target in development for _ in range(2)])
    development_features = np.zeros((len(development_ids), 229), dtype=np.float64)
    singleton_pos = {target: index for index, target in enumerate(centering)}
    panel_index = np.asarray(
        [singleton_pos.get(target, -1) for target in panel],
        dtype=np.int64,
    )
    projection_ratio = 0.995
    tail_value = np.sqrt(
        188.0 * (projection_ratio**-2 - 1.0) / (225.0 - 188.0)
    )
    panel_template = np.zeros(229, dtype=np.float64)
    panel_template[:188] = 1.0
    panel_template[188:225] = tail_value
    singleton_features[panel_index[panel_index >= 0]] = panel_template
    panel_singleton_features = singleton_features[panel_index[panel_index >= 0]]
    ordered_development = (
        *splits["test"],
        *splits["fit"],
        *splits["tune"],
    )
    development_guide_targets = tuple(
        target for target in ordered_development for _ in range(2)
    )
    test_labels = tuple(target for target in splits["test"] for _ in range(2))
    guide_targets = (
        *development_guide_targets,
        *(["FIXTURE"] * (7701 - len(development_guide_targets))),
    )
    guide_ids = tuple(f"guide-{index:05d}" for index in range(len(guide_targets)))
    guide_cache = directory / "source-cache" / "guide"
    feature_cache = directory / "source-cache" / "feature"
    guide_cache.mkdir(parents=True)
    feature_cache.mkdir(parents=True)
    np.save(guide_cache / "guide_delta.npy", np.zeros((1, 1)), allow_pickle=False)
    np.save(
        guide_cache / "target_pooled_delta.npy",
        np.zeros((1, 1)),
        allow_pickle=False,
    )
    np.save(
        guide_cache / "guide_cells.npy",
        np.ones(len(guide_targets), dtype=np.int64),
        allow_pickle=False,
    )
    np.save(
        guide_cache / "target_cells.npy",
        np.ones(1, dtype=np.int64),
        allow_pickle=False,
    )
    (guide_cache / "genes.txt").write_text("\n".join(genes) + "\n")
    (guide_cache / "guide_ids.txt").write_text("\n".join(guide_ids) + "\n")
    (guide_cache / "guide_targets.txt").write_text("\n".join(guide_targets) + "\n")
    (guide_cache / "target_names.txt").write_text(
        "\n".join(ordered_development) + "\n"
    )
    (guide_cache / "metadata.json").write_text("{}\n")
    feature_names = tuple(f"F{index:03d}" for index in range(229))
    null_features = np.zeros((len(guide_targets), len(feature_names)), dtype=np.float64)
    np.save(
        feature_cache / "guide_features.npy",
        null_features,
        allow_pickle=False,
    )
    np.save(
        feature_cache / "guide_identifiable.npy",
        np.ones(len(guide_targets), dtype=bool),
        allow_pickle=False,
    )
    (feature_cache / "feature_names.txt").write_text(
        "\n".join(feature_names) + "\n"
    )
    (feature_cache / "guide_ids.txt").write_text("\n".join(guide_ids) + "\n")
    (feature_cache / "seed-feature-provenance.json").write_text("{}\n")
    guide_members = FIT.source_cache_member_hashes(
        guide_cache,
        feature_cache,
        exact_members=True,
    )
    guide_axis_sha256 = FIT.sha256_lines(guide_ids)
    source_input = {
        "schema_version": FIT.SOURCE_INPUT_SCHEMA_VERSION,
        "experiment": FIT.EXPERIMENT,
        "guide_axis": {"count": 7701, "sha256": guide_axis_sha256},
        "gene_axis": {
            "count": gene_contract.count,
            "sha256": gene_contract.sha256,
        },
        "feature_axis": {
            "count": len(feature_names),
            "sha256": FIT.sha256_lines(feature_names),
        },
        "cache_locations": {
            "guide_cache": str(guide_cache.resolve()),
            "feature_cache": str(feature_cache.resolve()),
        },
        "guide_cache": {
            "metadata_sha256": guide_members["guide"]["metadata.json"],
            "guide_axis_sha256": guide_axis_sha256,
            "gene_axis_sha256": gene_contract.sha256,
            "outputs": {
                name: guide_members["guide"][name]
                for name in FIT.GUIDE_CACHE_OUTPUTS
            },
        },
        "feature_cache": {
            "provenance_sha256": guide_members["feature"][
                "seed-feature-provenance.json"
            ],
            "guide_axis_sha256": guide_axis_sha256,
            "feature_axis_sha256": FIT.sha256_lines(feature_names),
            "outputs": {
                name: guide_members["feature"][name]
                for name in FIT.FEATURE_CACHE_OUTPUTS
            },
        },
        "eligibility": {
            "predicate": FIT.FEATURE_ELIGIBILITY_PREDICATE,
            "identifiable_guides_count": 7701,
            "development_targets_count": len(development),
            "singleton_targets_count": len(centering),
            "panel_singletons_count": int(np.sum(panel_index >= 0)),
            "development_target_axis": {
                "count": development_contract.count,
                "sha256": development_contract.sha256,
            },
            "singleton_target_axis": {
                "count": centering_contract.count,
                "sha256": centering_contract.sha256,
            },
            "panel_target_axis": {
                "count": panel_contract.count,
                "sha256": panel_contract.sha256,
            },
        },
        "derived_arrays": {},
        "raw_sources": {
            "source_h5ad_sha256": source_axes["source"][
                "h5ad_sha256_from_pinned_source_audit"
            ],
            "source_h5ad_bytes": int(source_axes["source"]["h5ad_bytes"]),
            "seed_table_sha256": source_axes["source"]["seed_metadata"][
                "compressed_sha256"
            ],
            "panel_csv_sha256": source_axes["source"]["panel_csv_sha256"],
            "perturbmatch_perts_sha256": source_axes["source"][
                "perturbmatch_perts_sha256"
            ],
            "perturbmatch_genes_sha256": source_axes["source"][
                "perturbmatch_genes_sha256"
            ],
            "experiment_083_artifact_manifest_sha256": FIT.sha256_file(
                FIT.ROOT
                / "experiments/083-guide-latent-decontamination/artifact_manifest.json"
            ),
        },
        "source_code_sha256": {
            name: FIT.sha256_file(path)
            for name, path in FIT.SOURCE_FILES.items()
            if name.startswith("experiment_083_")
        },
    }
    derived = {
        "singleton_features": singleton_features,
        "panel_singleton_features": panel_singleton_features,
        "development_features": development_features,
        "development_target_ids": development_ids,
        "panel_singleton_source_index": panel_index,
    }
    source_input["derived_arrays"] = {
        name: {
            "shape": list(values.shape),
            "dtype": str(values.dtype),
            "sha256": FIT.sha256_array(values),
        }
        for name, values in derived.items()
    }
    FIT.validate_source_input_manifest(source_input)
    (directory / "source-input-manifest.json").write_text(
        json.dumps(source_input, indent=2, sort_keys=True) + "\n"
    )
    snapshot_receipt = {
        "schema_version": FIT.SOURCE_CACHE_SNAPSHOT_SCHEMA_VERSION,
        "experiment": FIT.EXPERIMENT,
        "guide_cache": "source-cache/guide",
        "feature_cache": "source-cache/feature",
        "members": guide_members,
        "source_manifest_content_sha256": FIT.sha256_json(
            FIT.source_manifest_content(source_input)
        ),
        "original_member_hashes_before": guide_members,
        "original_member_hashes_after": guide_members,
    }
    (directory / "source-cache-snapshot.json").write_text(
        json.dumps(snapshot_receipt, indent=2, sort_keys=True) + "\n"
    )

    truth = np.zeros((len(test_labels), len(response)), dtype=np.float64)
    truth[:, 0] = 1.0
    predicted = truth.copy()
    shrink_predicted = np.zeros_like(truth)
    with (directory / "source-test-holdouts.npz").open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(FIT.SOURCE_TEST_HOLDOUT_SCHEMA_VERSION),
            raw_residual=truth,
            recovered_nuisance=predicted,
            shrink_recovered_nuisance=shrink_predicted,
            target_labels=np.asarray(test_labels),
            held_guide_rows=np.arange(len(test_labels), dtype=np.int64),
            source_test_model_sha256=np.asarray(
                FIT.sha256_file(directory / "source-test-model.npz")
            ),
            source_input_manifest_sha256=np.asarray(
                FIT.sha256_file(directory / "source-input-manifest.json")
            ),
            response_axis_sha256=np.asarray(response_contract.sha256),
        )
    holdout_result = FIT.validate_source_test_holdouts(
        directory / "source-test-holdouts.npz",
        expected_test_targets=splits["test"],
        expected_guide_targets=guide_targets,
        expected_model_sha256=FIT.sha256_file(directory / "source-test-model.npz"),
        expected_source_input_sha256=FIT.sha256_file(
            directory / "source-input-manifest.json"
        ),
        expected_response_axis_sha256=response_contract.sha256,
        expected_response_width=len(response),
    )
    bootstrap = FIT.target_cluster_bootstrap(
        holdout_result["per_guide_cosine_delta"],
        holdout_result["labels"],
        repetitions=20_000,
        seed=MODEL.SPLIT_SEED,
    )
    within = np.full(100, 0.005)
    cross = np.full(100, 0.004)
    shrink = np.asarray([holdout_result["shrink_metrics"]["cosine_delta"]])
    with (directory / "source-bootstrap.npz").open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(1),
            seed=np.asarray(MODEL.SPLIT_SEED),
            cosine_delta=bootstrap,
            source_test_holdouts_sha256=np.asarray(
                FIT.sha256_file(directory / "source-test-holdouts.npz")
            ),
        )
    identifiable = np.ones(len(guide_targets), dtype=bool)
    target_strata, strata_audit = FIT.build_split_local_cross_target_strata(
        null_features,
        guide_targets,
        identifiable,
        splits,
    )
    within_permutations = np.empty((100, len(guide_targets)), dtype=np.int64)
    cross_permutations = np.empty((100, len(guide_targets)), dtype=np.int64)
    null_selection = {"schema_version": 2}
    for family, values in (
        ("within_target", within),
        ("cross_target", cross),
    ):
        family_records = []
        for replicate, value in enumerate(values):
            seed = MODEL.EXP083.stable_seed(
                MODEL.SPLIT_SEED,
                f"exp109-{family}-null",
                replicate,
            )
            permutation = (
                FIT.within_target_permutation_indices(
                    guide_targets,
                    seed=seed,
                    row_mask=identifiable,
                    targets=development,
                )
                if family == "within_target"
                else FIT.cross_target_permutation_indices(
                    guide_targets,
                    target_strata,
                    seed=seed,
                    row_mask=identifiable,
                )
            )
            if family == "within_target":
                within_permutations[replicate] = permutation
            else:
                cross_permutations[replicate] = permutation
            family_records.append(
                {
                "replicate": replicate,
                "seed": seed,
                "permutation_indices_sha256": FIT.sha256_array(permutation),
                "selected_candidate": FIT.candidate_payload(winner),
                "selected_tune_metrics": selection_records[
                    list(MODEL.candidate_grid()).index(winner)
                ]["metrics"],
                "candidate_records": selection_records,
                "source_test_metrics": {
                    "targets": len(splits["test"]),
                    "guide_holdouts": len(test_labels),
                    "cosine_delta": float(value),
                    "nmse_delta": 0.0,
                },
                "permuted_features_sha256": FIT.sha256_array(
                    null_features[permutation]
                ),
            }
            )
        null_selection[family] = family_records
    null_selection["cross_target_strata"] = strata_audit
    with (directory / "source-nulls.npz").open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(2),
            within_target_cosine_delta=within,
            cross_target_cosine_delta=cross,
            within_target_permutation_indices=within_permutations,
            cross_target_permutation_indices=cross_permutations,
            direction_preserving_shrink_delta=shrink,
        )
    (directory / "source-null-selection.json").write_text(
        json.dumps(
            null_selection,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    response_corrections = np.zeros((len(panel), len(response)), dtype=np.float64)
    slot_corrections = np.zeros((len(panel), len(genes)), dtype=np.float64)
    diagnostics = MODEL.stable_score_extrapolation_diagnostics(
        final_model,
        development_features,
        singleton_features,
    )
    with (directory / "singleton-corrections.npz").open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(FIT.SINGLETON_CORRECTION_SCHEMA_VERSION),
            singleton_target_ids=np.asarray(centering),
            singleton_features=singleton_features,
            panel_singleton_features=panel_singleton_features,
            development_features=development_features,
            development_target_ids=development_ids,
            panel_target_ids=np.asarray(panel),
            panel_singleton_source_index=panel_index,
            response_corrections=response_corrections,
            source_slot_corrections=slot_corrections,
            source_test_model_sha256=np.asarray(
                FIT.sha256_file(directory / "source-test-model.npz")
            ),
            selected_model_sha256=np.asarray(
                FIT.sha256_file(directory / "selected-model.npz")
            ),
            source_input_manifest_sha256=np.asarray(
                FIT.sha256_file(directory / "source-input-manifest.json")
            ),
            response_axis_sha256=np.asarray(response_contract.sha256),
            centering_axis_sha256=np.asarray(centering_contract.sha256),
            extrapolation_metric_names=np.asarray(tuple(diagnostics)),
            extrapolation_metric_values=np.asarray(tuple(diagnostics.values())),
        )
    support = panel_index >= 0
    pm_targets, _ = FIT.canonical_axis("perturbmatch_targets", source_axes)
    strata = np.zeros(len(panel), dtype=np.int64)
    strata[support] = 1
    strata[
        support & np.asarray([target not in set(pm_targets) for target in panel])
    ] = 2
    permutations = FIT.deterministic_stratified_permutations(
        strata,
        support,
        repetitions=100,
        seed=FIT.HELD_LABEL_NULL_SEED,
    )
    with (directory / "held-label-nulls.npz").open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(FIT.HELD_LABEL_NULL_SCHEMA_VERSION),
            seed=np.asarray(FIT.HELD_LABEL_NULL_SEED),
            permutations=permutations,
            strata=strata,
            movable=support,
            panel_targets_sha256=np.asarray(FIT.sha256_lines(panel)),
            singleton_corrections_sha256=np.asarray(
                FIT.sha256_file(directory / "singleton-corrections.npz")
            ),
        )
    guard = FIT.validate_preregistration()["source_only_runner"]["memory_guard"]
    memory_receipt = {
        "schema_version": FIT.MEMORY_GUARD_RECEIPT_SCHEMA_VERSION,
        "experiment": FIT.EXPERIMENT,
        "status": "complete",
        "scope": guard["scope"],
        "guard_script": str(FIT.MEMORY_GUARD_PATH.relative_to(FIT.ROOT)),
        "guard_script_sha256": FIT.MEMORY_GUARD_SHA256,
        "root_pid": 100,
        "monitor_pid": 101,
        "max_rss_gib": guard["max_rss_gib"],
        "min_free_percent": guard["min_free_percent"],
        "poll_seconds": guard["poll_seconds"],
        "peak_rss_gib": 1.0,
        "samples": 1,
        "last_process_count": 1,
        "minimum_free_percent": 50,
        "elapsed_seconds": 1.0,
        "source_manifest_content_sha256": snapshot_receipt[
            "source_manifest_content_sha256"
        ],
        "source_cache_snapshot_sha256": FIT.sha256_file(
            directory / "source-cache-snapshot.json"
        ),
    }
    (directory / "memory-guard.json").write_text(
        json.dumps(memory_receipt, indent=2, sort_keys=True) + "\n"
    )
    source_test = FIT.source_test_decision(
        holdout_result["metrics"],
        bootstrap_samples=bootstrap,
        within_target_null=within,
        cross_target_null=cross,
        direction_preserving_shrink_delta=float(shrink[0]),
    )
    source_test["selected_candidate"] = FIT.candidate_payload(winner)
    source_test["final_model_invariants"] = (
        FIT.validate_preregistered_final_invariants(
            final_model,
            source_derived=derived,
            source_input_manifest=source_input,
        )
    )
    source_test["held_evaluation_authorized"] = False
    source_test["bindings"] = {
        "source_input_manifest_sha256": FIT.sha256_file(
            directory / "source-input-manifest.json"
        ),
        "source_cache_snapshot_sha256": FIT.sha256_file(
            directory / "source-cache-snapshot.json"
        ),
        "memory_guard_sha256": FIT.sha256_file(directory / "memory-guard.json"),
        "source_test_model_sha256": FIT.sha256_file(
            directory / "source-test-model.npz"
        ),
        "source_test_holdouts_sha256": FIT.sha256_file(
            directory / "source-test-holdouts.npz"
        ),
        "selected_model_sha256": FIT.sha256_file(directory / "selected-model.npz"),
        "tune_selection_sha256": FIT.sha256_file(directory / "tune-selection.json"),
        "source_bootstrap_sha256": FIT.sha256_file(directory / "source-bootstrap.npz"),
        "source_nulls_sha256": FIT.sha256_file(directory / "source-nulls.npz"),
        "source_null_selection_sha256": FIT.sha256_file(
            directory / "source-null-selection.json"
        ),
        "axes_sha256": FIT.sha256_file(directory / "axes.json"),
        "singleton_corrections_sha256": FIT.sha256_file(
            directory / "singleton-corrections.npz"
        ),
        "held_label_nulls_sha256": FIT.sha256_file(directory / "held-label-nulls.npz"),
        "selected_candidate": FIT.candidate_payload(winner),
        "split_hashes": {
            name: FIT.sha256_lines(splits[name]) for name in FIT.REQUIRED_SPLITS
        },
    }
    (directory / "source-test.json").write_text(
        json.dumps(source_test, indent=2, sort_keys=True) + "\n"
    )
    for filename in FIT.REQUIRED_FROZEN_ARTIFACTS:
        path = directory / filename
        if path.exists():
            continue
        with path.open("wb") as handle:
            np.savez_compressed(handle, marker=np.asarray([1]))
    manifest = FIT.build_freeze_manifest(
        directory,
        splits=splits,
        frozen_files=FIT.REQUIRED_FROZEN_ARTIFACTS,
        source_test=source_test,
        response_axis=response_axis,
        centering_target_ids=centering,
        centering_axis_contract=centering_contract,
    )
    manifest_path = directory / "freeze-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    authorization = {
        "authorized": True,
        "experiment": FIT.EXPERIMENT,
        "purpose": "held_context_evaluation_only",
        "no_held_truth_opened_during_fit": True,
        "freeze_manifest_sha256": FIT.sha256_file(manifest_path),
    }
    (directory / "held-evaluation-authorized.json").write_text(
        json.dumps(authorization, indent=2, sort_keys=True) + "\n"
    )
    return manifest, authorization


def test_freeze_authorization_rejects_failed_source_and_missing_artifacts() -> None:
    with tempfile.TemporaryDirectory() as raw:
        directory = Path(raw)
        manifest, authorization = _full_freeze_fixture(directory)
        try:
            EVALUATE.verify_frozen_authorization(directory)
        except RuntimeError as exc:
            assert "deliberately disabled" in str(exc)
        else:
            raise AssertionError("held evaluation unexpectedly became reachable")
        EVALUATE._verify_frozen_authorization_disabled(directory)

        failed = json.loads(json.dumps(manifest))
        failed["source_test"]["passes"] = False
        (directory / "freeze-manifest.json").write_text(
            json.dumps(failed, indent=2, sort_keys=True) + "\n"
        )
        authorization["freeze_manifest_sha256"] = FIT.sha256_file(
            directory / "freeze-manifest.json"
        )
        (directory / "held-evaluation-authorized.json").write_text(
            json.dumps(authorization, indent=2, sort_keys=True) + "\n"
        )
        try:
            EVALUATE._verify_frozen_authorization_disabled(directory)
        except ValueError as exc:
            assert "source test did not pass" in str(exc)
        else:
            raise AssertionError("failed source test authorized held evaluation")

        missing = json.loads(json.dumps(manifest))
        missing["frozen_artifacts"].pop("held-label-nulls.npz")
        (directory / "freeze-manifest.json").write_text(
            json.dumps(missing, indent=2, sort_keys=True) + "\n"
        )
        authorization["freeze_manifest_sha256"] = FIT.sha256_file(
            directory / "freeze-manifest.json"
        )
        (directory / "held-evaluation-authorized.json").write_text(
            json.dumps(authorization, indent=2, sort_keys=True) + "\n"
        )
        try:
            EVALUATE._verify_frozen_authorization_disabled(directory)
        except ValueError as exc:
            assert "artifact set" in str(exc)
        else:
            raise AssertionError("incomplete artifact set authorized held evaluation")


def _rewrite_npz(path: Path, values: dict[str, np.ndarray]) -> None:
    with path.open("wb") as handle:
        np.savez_compressed(handle, **values)


def test_hardened_artifacts_reject_every_known_forgery() -> None:
    with tempfile.TemporaryDirectory() as raw:
        directory = Path(raw)
        manifest, _ = _full_freeze_fixture(directory)
        splits = {
            name: tuple(manifest["splits"][name]["targets"])
            for name in FIT.REQUIRED_SPLITS
        }
        source_test = manifest["source_test"]

        forged = json.loads(json.dumps(source_test))
        forged["metrics"]["cosine_delta"] += 0.1
        try:
            FIT.validate_source_evidence(directory, forged, splits)
        except ValueError as exc:
            assert "metric was not reproduced" in str(exc)
        else:
            raise AssertionError("forged JSON source metrics were accepted")

        source_model_path = directory / "source-test-model.npz"
        original_source_model = source_model_path.read_bytes()
        source_model = FIT.load_frozen_model(source_model_path)
        FIT.save_frozen_model(
            source_model_path,
            dataclasses.replace(source_model, fit_targets_sha256="f" * 64),
            role="source_test_model",
        )
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "fit-target axis" in str(exc)
        else:
            raise AssertionError("wrong source-test training population was accepted")
        source_model_path.write_bytes(original_source_model)

        bootstrap_path = directory / "source-bootstrap.npz"
        original_bootstrap = bootstrap_path.read_bytes()
        with np.load(bootstrap_path, allow_pickle=False) as archive:
            bootstrap_values = {
                name: np.asarray(archive[name]) for name in archive.files
            }
        bootstrap_values["cosine_delta"] = bootstrap_values["cosine_delta"].copy()
        bootstrap_values["cosine_delta"][0] += 0.01
        _rewrite_npz(bootstrap_path, bootstrap_values)
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except AssertionError:
            pass
        else:
            raise AssertionError(
                "bootstrap not derived from frozen holdouts was accepted"
            )
        bootstrap_path.write_bytes(original_bootstrap)

        holdout_path = directory / "source-test-holdouts.npz"
        original_holdout = holdout_path.read_bytes()
        with np.load(holdout_path, allow_pickle=False) as archive:
            holdout_values = {name: np.asarray(archive[name]) for name in archive.files}
        holdout_values["recovered_nuisance"] = holdout_values[
            "recovered_nuisance"
        ].copy()
        holdout_values["recovered_nuisance"][0, 0] = 0.0
        _rewrite_npz(holdout_path, holdout_values)
        rebound = json.loads(json.dumps(source_test))
        rebound["bindings"]["source_test_holdouts_sha256"] = FIT.sha256_file(
            holdout_path
        )
        try:
            FIT.validate_source_evidence(directory, rebound, splits)
        except ValueError as exc:
            assert "metric was not reproduced" in str(exc)
        else:
            raise AssertionError("mutated per-guide outputs were accepted")
        holdout_path.write_bytes(original_holdout)

        source_input_path = directory / "source-input-manifest.json"
        original_source_input = source_input_path.read_text()
        bad_input = json.loads(original_source_input)
        bad_input["feature_cache"]["guide_axis_sha256"] = "0" * 64
        source_input_path.write_text(json.dumps(bad_input, indent=2) + "\n")
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "shared guide axis" in str(exc)
        else:
            raise AssertionError("reordered feature-cache guide axis was accepted")
        source_input_path.write_text(original_source_input)

        bad_input = json.loads(original_source_input)
        bad_input["eligibility"]["predicate"] = "nonempty_strings"
        source_input_path.write_text(json.dumps(bad_input, indent=2) + "\n")
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "normalized-sequence" in str(exc)
        else:
            raise AssertionError("weaker feature-eligibility semantics were accepted")
        source_input_path.write_text(original_source_input)

        snapshot_member = directory / "source-cache/feature/feature_names.txt"
        original_snapshot_member = snapshot_member.read_bytes()
        snapshot_member.write_bytes(original_snapshot_member + b"FORGED\n")
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "snapshot member hashes changed" in str(exc)
        else:
            raise AssertionError("mutated embedded source cache was accepted")
        snapshot_member.write_bytes(original_snapshot_member)

        before = FIT.source_cache_member_hashes(
            directory / "source-cache/guide",
            directory / "source-cache/feature",
            exact_members=True,
        )
        snapshot_member.write_bytes(original_snapshot_member + b"CHANGED\n")
        try:
            SOURCE_ONLY._finalize_source_cache_snapshot(
                guide_cache=directory / "source-cache/guide",
                feature_cache=directory / "source-cache/feature",
                staging=directory,
                original_manifest=json.loads(original_source_input),
                before=before,
            )
        except ValueError as exc:
            assert "changed during source computation" in str(exc)
        else:
            raise AssertionError("post-load cache mutation was accepted")
        snapshot_member.write_bytes(original_snapshot_member)

        memory_path = directory / "memory-guard.json"
        original_memory = memory_path.read_text()
        forged_memory = json.loads(original_memory)
        forged_memory["status"] = "limit_exceeded"
        memory_path.write_text(json.dumps(forged_memory, indent=2) + "\n")
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "did not complete" in str(exc)
        else:
            raise AssertionError("forged memory-guard receipt was accepted")
        memory_path.write_text(original_memory)

        source_nulls_path = directory / "source-nulls.npz"
        original_source_nulls = source_nulls_path.read_bytes()
        with np.load(source_nulls_path, allow_pickle=False) as archive:
            source_null_values = {
                name: np.asarray(archive[name]) for name in archive.files
            }
        source_null_values["cross_target_permutation_indices"] = source_null_values[
            "cross_target_permutation_indices"
        ].copy()
        source_null_values["cross_target_permutation_indices"][0, 0] = 0
        _rewrite_npz(source_nulls_path, source_null_values)
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except AssertionError:
            pass
        else:
            raise AssertionError("forged source-null assignment was accepted")
        source_nulls_path.write_bytes(original_source_nulls)

        null_selection_path = directory / "source-null-selection.json"
        original_null_selection = null_selection_path.read_text()
        forged_selection = json.loads(original_null_selection)
        forged_selection["cross_target"][0]["candidate_records"].pop()
        null_selection_path.write_text(
            json.dumps(forged_selection, indent=2, sort_keys=True) + "\n"
        )
        try:
            FIT.validate_source_evidence(directory, source_test, splits)
        except ValueError as exc:
            assert "candidate selection grid mismatch" in str(exc)
        else:
            raise AssertionError("incomplete null candidate grid was accepted")
        null_selection_path.write_text(original_null_selection)

        correction_path = directory / "singleton-corrections.npz"
        original_correction = correction_path.read_bytes()
        _rewrite_npz(correction_path, {"marker": np.asarray([1])})
        rebound = json.loads(json.dumps(source_test))
        rebound["bindings"]["singleton_corrections_sha256"] = FIT.sha256_file(
            correction_path
        )
        try:
            FIT.validate_source_evidence(directory, rebound, splits)
        except ValueError as exc:
            assert "singleton-correction field set" in str(exc)
        else:
            raise AssertionError("marker-only singleton corrections were accepted")
        correction_path.write_bytes(original_correction)

        with np.load(correction_path, allow_pickle=False) as archive:
            correction_values = {
                name: np.asarray(archive[name]) for name in archive.files
            }
        correction_values["extrapolation_metric_values"] = correction_values[
            "extrapolation_metric_values"
        ].copy()
        correction_values["extrapolation_metric_values"][0] += 1.0
        _rewrite_npz(correction_path, correction_values)
        rebound = json.loads(json.dumps(source_test))
        rebound["bindings"]["singleton_corrections_sha256"] = FIT.sha256_file(
            correction_path
        )
        try:
            FIT.validate_source_evidence(directory, rebound, splits)
        except AssertionError:
            pass
        else:
            raise AssertionError("forged extrapolation diagnostics were accepted")
        correction_path.write_bytes(original_correction)

        with np.load(correction_path, allow_pickle=False) as archive:
            correction_values = {
                name: np.asarray(archive[name]) for name in archive.files
            }
        correction_values["panel_singleton_features"] = correction_values[
            "panel_singleton_features"
        ].copy()
        correction_values["panel_singleton_features"][0, 0] += 1.0
        _rewrite_npz(correction_path, correction_values)
        rebound = json.loads(json.dumps(source_test))
        rebound["bindings"]["singleton_corrections_sha256"] = FIT.sha256_file(
            correction_path
        )
        try:
            FIT.validate_source_evidence(directory, rebound, splits)
        except (AssertionError, ValueError):
            pass
        else:
            raise AssertionError(
                "arbitrary panel-singleton feature matrix was accepted"
            )
        correction_path.write_bytes(original_correction)

        null_path = directory / "held-label-nulls.npz"
        with np.load(null_path, allow_pickle=False) as archive:
            null_values = {name: np.asarray(archive[name]) for name in archive.files}
        null_values["permutations"] = null_values["permutations"].copy()
        null_values["permutations"][0, 0] = null_values["permutations"][0, 1]
        _rewrite_npz(null_path, null_values)
        rebound = json.loads(json.dumps(source_test))
        rebound["bindings"]["held_label_nulls_sha256"] = FIT.sha256_file(null_path)
        try:
            FIT.validate_source_evidence(directory, rebound, splits)
        except AssertionError:
            pass
        else:
            raise AssertionError("non-bijective held-label null was accepted")


def test_zero_mean_stage_is_explicitly_pre_perturbmatch() -> None:
    protocol = FIT.validate_preregistration()
    assert (
        protocol["production_adapter"]["post_perturbmatch_zero_mean_claimed"] is False
    )
    correction = np.asarray([[1.0, -2.0], [-1.0, 2.0]])
    assert np.max(np.abs(correction.sum(axis=0))) == 0.0
    retained = MODEL.apply_k562_slot_retention(
        correction,
        np.asarray([True, False]),
        np.asarray([True, True]),
    )
    assert np.max(np.abs(retained.sum(axis=0))) > 0.0


def test_hct116_direction_is_not_a_promotion_input() -> None:
    metrics = {
        "source_test_passed": True,
        "source_bootstrap_passed": True,
        "source_nulls_passed": True,
        "correction_mean_abs_max": 0.0,
        "correction_outside_supported_singletons_abs_max": 0.0,
        "feature_nullspace_prediction_change_abs_max": 0.0,
        "excluded_coordinate_abs_max": 0.0,
        "eta_zero_exact_identity": True,
        "h1_raw_pds_delta": 0.003,
        "hct116_raw_pds_delta": 0.0,
        "hek293t_raw_pds_delta": 0.003,
        "h1_directional_credit_delta": 0.001,
        "hek293t_directional_credit_delta": 0.001,
        "macro_raw_pds_delta": 0.003,
        "singleton_macro_raw_pds_delta": 0.003,
        "above_frozen_label_null_p95": True,
        "exact_baseline_reconstruction": True,
        "unchanged_target_coverage": True,
        "hct116_directional_credit_delta": -999.0,
    }
    result = EVALUATE.promotion_decision(metrics)
    assert result["all_pass"] is True
    assert "hct116_directional" not in result["checks"]


def test_source_runner_has_no_held_or_submission_path() -> None:
    text = (HERE / "source_only.py").read_text()
    for forbidden in (
        "vcc submit",
        "cell-eval2",
        "anndata.read",
        "scanpy.read",
        "boto3",
        "kolf",
        "import evaluate",
    ):
        assert forbidden not in text


def test_memory_guard_wraps_the_complete_heavy_source_phase() -> None:
    text = (HERE / "source_only.py").read_text()
    run_body = text[text.index("def run_source_only") :]
    assert run_body.index("monitor, guard_paths = _start_memory_guard") < run_body.index(
        "initial = preflight"
    )
    assert run_body.index("initial = preflight") < run_body.index(
        "generate_or_verify_caches"
    )
    assert run_body.index("_run_source_nulls") < run_body.index(
        "_finalize_source_cache_snapshot"
    )
    assert run_body.index("_finalize_source_cache_snapshot") < run_body.index(
        "_finish_memory_guard"
    )
    assert run_body.index("_finish_memory_guard") < run_body.index("_publish")


def test_source_preflight_fails_closed_on_missing_verified_caches() -> None:
    parser = SOURCE_ONLY.build_parser()
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw)
        args = parser.parse_args(
            [
                "--preflight",
                "--cache-mode",
                "verify",
                "--guide-cache",
                str(root / "guide"),
                "--feature-cache",
                str(root / "feature"),
                "--output",
                str(root / "output"),
            ]
        )
        SOURCE_ONLY._resolve_paths(args)
        result = SOURCE_ONLY.preflight(args)
        assert result["status"] == "BLOCKED"
        assert result["held_evaluation_enabled"] is False
        assert result["held_outcomes_opened"] is False


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
