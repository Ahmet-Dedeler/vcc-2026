#!/usr/bin/env python3
"""Synthetic tests for experiment 083; no donor matrices or network access."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MODEL = load("exp083_model", "model.py")
EXTRACT = load("exp083_extract", "extract_guide_effects.py")
FEATURES = load("exp083_features", "build_seed_features.py")


def synthetic_latent_data():
    rng = np.random.default_rng(83)
    n_targets = 72
    guides_per_target = 3
    n_features = 14
    n_genes = 48
    rank = 4
    targets = tuple(
        f"T{target:03d}"
        for target in range(n_targets)
        for _ in range(guides_per_target)
    )
    theta = rng.normal(size=(n_targets, n_genes))
    features = rng.normal(size=(len(targets), n_features))
    left = rng.normal(size=(n_features, rank))
    right, _ = np.linalg.qr(rng.normal(size=(n_genes, rank)))
    nuisance = 1.6 * features @ left @ right.T
    effects = np.repeat(theta, guides_per_target, axis=0) + nuisance
    effects += rng.normal(scale=0.02, size=effects.shape)
    counts = rng.integers(50, 400, size=len(targets))
    identifiable = np.ones(len(targets), dtype=bool)
    genes = tuple(f"G{i}" for i in range(n_genes))
    return effects, features, targets, counts, identifiable, genes


def test_split_and_coordinate_contract() -> None:
    splits = MODEL.deterministic_target_split([f"T{i}" for i in range(40)], seed=83)
    MODEL.assert_split_contract(splits, ["P0", "P1"])
    assert set(splits["fit"]).isdisjoint(splits["tune"])
    assert set(splits["fit"]).isdisjoint(splits["test"])
    mask = MODEL.excluded_coordinate_mask(["A", "P0", "B", "P1"], ["P0", "P1"])
    np.testing.assert_array_equal(mask, [True, False, True, False])


def test_latent_model_predicts_held_guides() -> None:
    effects, features, targets, counts, identifiable, genes = synthetic_latent_data()
    splits = MODEL.deterministic_target_split(sorted(set(targets)), seed=83)
    candidate = MODEL.Candidate(
        rank=4, ridge=0.1, eta=1.0, weight_power=0.0, huber_c=2.0
    )
    fitted = MODEL.fit_latent_model(
        effects,
        features,
        targets,
        identifiable,
        splits["fit"],
        np.ones(len(genes), dtype=bool),
        candidate,
        seed=83,
    )
    truth, baseline, predicted, labels = MODEL.held_guide_predictions(
        fitted,
        effects,
        features,
        targets,
        counts,
        identifiable,
        splits["tune"],
        seed=83,
    )
    metrics = MODEL.guide_prediction_metrics(truth, baseline, predicted)
    assert len(labels) == len(splits["tune"])
    assert metrics["cosine_delta"] > 0.10
    assert metrics["nmse_delta"] < 0


def test_single_or_unidentifiable_target_is_bit_exact() -> None:
    effects, features, targets, counts, identifiable, genes = synthetic_latent_data()
    candidate = MODEL.Candidate(
        rank=4, ridge=1.0, eta=1.0, weight_power=0.0, huber_c=2.0
    )
    fitted = MODEL.fit_latent_model(
        effects,
        features,
        targets,
        identifiable,
        sorted(set(targets))[:40],
        np.ones(len(genes), dtype=bool),
        candidate,
        seed=83,
    )
    pooled = effects[0].copy()
    cleaned, audit = MODEL.infer_cleaned_target(
        fitted,
        effects[:1],
        features[:1],
        counts[:1],
        identifiable[:1],
        pooled,
    )
    assert not audit["modified"]
    assert cleaned.tobytes() == pooled.tobytes()
    cleaned, audit = MODEL.infer_cleaned_target(
        fitted,
        effects[:3],
        features[:3],
        counts[:3],
        np.zeros(3, dtype=bool),
        pooled,
    )
    assert not audit["modified"]
    assert cleaned.tobytes() == pooled.tobytes()


def test_no_whole_row_deletion_and_arm_d_scaling() -> None:
    raw = np.asarray([2.0, 4.0, -1.0, 3.0])
    shared = np.asarray([True, False, True, False])
    adjusted = MODEL.arm_d_k562_adjustment(
        raw,
        effect_scale=6.0,
        raw_median_norm=3.0,
        pm_shared_mask=shared,
        pm_available=True,
        companion_available=True,
    )
    # 6/3 source scaling, 1/2 PM blend on shared coordinates, 1/2
    # equal-source arm-D averaging.
    np.testing.assert_allclose(adjusted, [1.0, 4.0, -0.5, 3.0])
    absent = MODEL.arm_d_k562_adjustment(
        raw,
        effect_scale=6.0,
        raw_median_norm=3.0,
        pm_shared_mask=shared,
        pm_available=False,
        companion_available=True,
    )
    np.testing.assert_array_equal(absent, np.zeros_like(raw))


def test_permutation_nulls_preserve_strata_and_rows() -> None:
    features = np.arange(36, dtype=np.float64).reshape(12, 3)
    targets = ("A", "A", "B", "B", "C", "C", "D", "D", "E", "E", "F", "F")
    norms = {target: 1.0 for target in set(targets)}
    bundle = MODEL.stratified_feature_permutation(
        features, targets, norms, seed=83, bundle=True
    )
    rowwise = MODEL.stratified_feature_permutation(
        features, targets, norms, seed=84, bundle=False
    )
    assert sorted(map(tuple, bundle.tolist())) == sorted(map(tuple, features.tolist()))
    assert sorted(map(tuple, rowwise.tolist())) == sorted(map(tuple, features.tolist()))


def test_source_clusters_and_bootstrap_are_deterministic() -> None:
    source = np.asarray([[1, 0], [0.9, 0.1], [0, 1], [0.1, 0.9]], dtype=np.float64)
    clusters = MODEL.source_effect_clusters(source, k=1)
    assert len(np.unique(clusters)) == 2
    first = MODEL.clustered_bootstrap(
        np.asarray([0.1, 0.2, -0.1, 0.0]),
        clusters,
        repetitions=100,
        seed=83,
    )
    second = MODEL.clustered_bootstrap(
        np.asarray([0.1, 0.2, -0.1, 0.0]),
        clusters,
        repetitions=100,
        seed=83,
    )
    assert first == second


def test_grouped_sum_and_guide_parser() -> None:
    block = np.asarray([[1, 2], [3, 4], [5, 6], [7, 8]], dtype=np.float32)
    groups, sums = EXTRACT.grouped_sum_block(
        block, np.asarray([2, 1, 2, -1], dtype=np.int32)
    )
    np.testing.assert_array_equal(groups, [1, 2])
    np.testing.assert_allclose(sums, [[3, 4], [6, 8]])
    values, names = FEATURES.structural_features("ABC_+_100.23-P1|ABC_-_130.23-P1")
    assert values[names.index("pair_members")] == 2
    assert values[names.index("opposite_strand")] == 1
    assert values[names.index("parse_complete")] == 1


def test_low_cell_guide_does_not_change_target_reconstruction() -> None:
    targets = ("TARGET", "control")
    guides = ("target-high", "target-low", "control")
    target_codes = np.asarray([0] * 40 + [0] * 17 + [1] * 20, dtype=np.int32)
    guide_codes = np.asarray([0] * 40 + [1] * 17 + [2] * 20, dtype=np.int32)
    block = np.vstack(
        (
            np.tile([2.0, 0.0], (40, 1)),
            np.tile([0.0, 3.0], (17, 1)),
            np.tile([1.0, 1.0], (20, 1)),
        )
    ).astype(np.float32)

    guide_target, guide_cells = EXTRACT.assign_guide_targets(
        target_codes,
        guide_codes,
        len(guides),
        block_rows=11,
    )
    target_cells = EXTRACT.category_counts(target_codes, len(targets), block_rows=13)
    _, target_output_rows, output_target_cells = EXTRACT.compact_target_layout(
        targets,
        target_cells,
        ("TARGET",),
        control_code=1,
        min_cells=20,
    )
    valid_guides = np.flatnonzero(
        (guide_target >= 0) & (guide_target != 1) & (guide_cells >= 20)
    )
    np.testing.assert_array_equal(valid_guides, [0])
    np.testing.assert_array_equal(output_target_cells, [57])

    guide_groups, guide_sums = EXTRACT.grouped_sum_block(block, guide_codes)
    target_groups, target_sums = EXTRACT.grouped_sum_block(
        block,
        target_output_rows[target_codes],
    )
    np.testing.assert_array_equal(guide_groups, [0, 1, 2])
    np.testing.assert_array_equal(target_groups, [0])
    np.testing.assert_array_equal(guide_sums[valid_guides].sum(axis=0), [80.0, 0.0])
    np.testing.assert_array_equal(target_sums[0], [80.0, 51.0])

    control_sum = block[target_codes == 1].astype(np.float64).sum(axis=0)
    b_ctrl = np.log1p(EXTRACT.TARGET_SUM * control_sum / control_sum.sum())
    reconstructed = (
        np.log1p(
            EXTRACT.TARGET_SUM
            * target_sums[0]
            / float(target_sums[0].sum())
        )
        - b_ctrl
    )
    all_target_rows = block[target_codes == 0].astype(np.float64).sum(axis=0)
    expected = (
        np.log1p(
            EXTRACT.TARGET_SUM
            * all_target_rows
            / float(all_target_rows.sum())
        )
        - b_ctrl
    )
    filtered_only = (
        np.log1p(
            EXTRACT.TARGET_SUM
            * guide_sums[valid_guides].sum(axis=0)
            / float(guide_sums[valid_guides].sum())
        )
        - b_ctrl
    )
    np.testing.assert_array_equal(reconstructed, expected)
    assert not np.allclose(filtered_only, expected)


def test_no_experiment_080_dependency_or_submission_command() -> None:
    text = "\n".join(
        (HERE / name).read_text()
        for name in (
            "model.py",
            "inspect_source.py",
            "build_seed_features.py",
            "extract_guide_effects.py",
            "fit.py",
            "evaluate.py",
        )
    )
    forbidden = ("080-listwise", "vcc submit", "cell-eval2")
    for needle in forbidden:
        assert needle not in text


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
