#!/usr/bin/env python3
"""Synthetic and static contract tests for Experiment 145."""

from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_PATH = HERE / "run.py"


def load_run():
    spec = importlib.util.spec_from_file_location("experiment145_tested", RUN_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


M = load_run()


@contextmanager
def temporary_constants(**values):
    old = {name: getattr(M, name) for name in values}
    try:
        for name, value in values.items():
            setattr(M, name, value)
        yield
    finally:
        for name, value in old.items():
            setattr(M, name, value)


def score(value: float, size: int = 7) -> dict[str, dict[str, np.ndarray]]:
    row = np.full(size, value, dtype=np.float64)
    return {
        "k562_gwps": {"pds": row, "signed_credit": row},
        "k562_essential": {"pds": row, "signed_credit": row},
        "rpe1": {"pds": row, "signed_credit": row},
        "nadig_jurkat": {"pds": row, "signed_credit": row},
    }


def synthetic_cache():
    rng = np.random.default_rng(145)
    baseline = rng.normal(size=(4, 6)).astype(np.float32)
    transformed = rng.normal(size=(4, 6)).astype(np.float32)
    predictions = np.stack(
        [baseline, transformed, transformed * 0.9, transformed * 0.8, transformed * 0.7]
    )
    candidates = {
        "k562": {"baseline": baseline, "predictions": predictions},
        "non_k562": {"baseline": baseline, "predictions": predictions},
    }
    truth = {
        dataset: rng.normal(size=(4, 6)).astype(np.float64)
        for dataset in M.SCORED_DATASETS
    }
    parent = SimpleNamespace(
        DATASET_GROUP={
            "k562_gwps": "k562",
            "k562_essential": "k562",
            "rpe1": "non_k562",
            "nadig_jurkat": "non_k562",
        },
        signed_credit=lambda prediction, observed: (
            np.sum(np.sign(prediction) * np.sign(observed), axis=1)
            / prediction.shape[1]
        ),
    )
    folds = np.asarray([0, 0, 1, 1], dtype=np.int8)
    return M.build_score_cache(candidates, truth, parent, folds), folds


def test_protocol_matches_code_constants() -> None:
    value = json.loads((HERE / "preregistered_protocol.json").read_text())
    assert tuple(value["candidate"]["betas"]) == M.BETAS
    for group in M.GROUPS:
        assert (
            tuple(value["candidate"]["descending_trace_boundaries_one_based"][group])
            == M.TRACE_BOUNDARIES[group]
        )


def test_trace_boundaries_match_frozen_real_operators() -> None:
    for group in M.GROUPS:
        with np.load(
            M.EXP140_IDENTITY / f"whitener-{group}.npz", allow_pickle=False
        ) as data:
            assert M.trace_boundaries(data["evals"]) == M.TRACE_BOUNDARIES[group]


def test_frozen_plan_has_fail_closed_structure() -> None:
    with np.load(M.EXP142_PLAN, allow_pickle=False) as data:
        M.validate_plan_arrays(
            data["development_folds"],
            data["development_blocks"],
            data["bootstrap_block_draws"],
        )


def test_band_labels_preserve_exact_sizes() -> None:
    for group in M.GROUPS:
        labels = M.band_labels_ascending(M.TRACE_BOUNDARIES[group])
        b1, b2, b3, end = M.TRACE_BOUNDARIES[group]
        expected = np.asarray([b1, b2 - b1, b3 - b2, end - b3])
        np.testing.assert_array_equal(np.bincount(labels, minlength=4), expected)


def test_beta_zero_is_bit_exact() -> None:
    rng = np.random.default_rng(1)
    baseline = rng.normal(size=(9, M.RESPONSE_COUNT)).astype(np.float32)
    output = M.row_norm_preserving_spectral_transform(
        baseline,
        np.eye(M.RESPONSE_COUNT, dtype=np.float32),
        M.band_labels_ascending(M.TRACE_BOUNDARIES["k562"]),
        0.0,
    )
    assert np.array_equal(output, baseline)


def test_positive_beta_preserves_row_norm() -> None:
    rng = np.random.default_rng(2)
    dimension = 31
    q, _ = np.linalg.qr(rng.normal(size=(dimension, dimension)))
    baseline = rng.normal(size=(13, dimension)).astype(np.float32)
    old = M.RESPONSE_COUNT
    try:
        M.RESPONSE_COUNT = dimension
        labels = np.repeat(np.arange(4, dtype=np.uint8), [4, 6, 9, 12])
        output = M.row_norm_preserving_spectral_transform(
            baseline, q.astype(np.float32), labels, 0.75
        )
    finally:
        M.RESPONSE_COUNT = old
    np.testing.assert_allclose(
        np.linalg.norm(output, axis=1),
        np.linalg.norm(baseline, axis=1),
        rtol=2e-6,
        atol=2e-6,
    )


def test_tail_multipliers_are_monotone() -> None:
    old = M.RESPONSE_COUNT
    try:
        M.RESPONSE_COUNT = 8
        labels = np.asarray([3, 3, 2, 2, 1, 1, 0, 0], dtype=np.uint8)
        values = M.multipliers(labels, 0.5)
    finally:
        M.RESPONSE_COUNT = old
    assert np.all(np.diff(values) >= 0.0)
    assert values[0] < values[2] < values[4] == values[-1] == 1.0


def test_pds_square_similarity_uses_midrank() -> None:
    similarity = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.9, 0.2, 0.1],
        ]
    )
    np.testing.assert_allclose(
        M.pds_from_square_similarity(similarity), [1.0, 0.75, 0.0]
    )


def test_selector_has_identity_fallback_and_prefers_pds() -> None:
    scores = {
        0.0: score(0.0),
        0.25: score(-0.001),
        0.5: score(0.002),
        0.75: score(0.001),
        1.0: score(-0.002),
    }
    selected, _, eligible = M.select_beta(scores)
    assert selected == 0.5
    assert 0.0 in eligible
    assert 0.25 not in eligible


def test_observed_and_bootstrap_share_exact_tie_breaking() -> None:
    assert bool(M.prefer_candidate(1.0 + 0.5 * M.EPS, 0.0, 0.25, 1.0, 5.0, 0.0))
    assert not bool(M.prefer_candidate(1.0, 5.0, 0.25, 1.0, 5.0, 0.0))
    observed = M.prefer_candidate(
        np.asarray([1.0 + 0.5 * M.EPS, 1.0]),
        np.asarray([0.0, 5.0]),
        0.25,
        np.asarray([1.0, 1.0]),
        np.asarray([5.0, 5.0]),
        np.asarray([0.0, 0.0]),
    )
    np.testing.assert_array_equal(observed, [True, False])


def test_inner_gallery_is_not_a_subset_of_full_gallery_scores() -> None:
    similarity = np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.9, 0.1, 0.8, 0.7],
            [0.0, 0.0, 0.1, 0.8],
        ]
    )
    inner_indices = np.asarray([2, 3])
    inner = M.pds_from_square_similarity(
        similarity[np.ix_(inner_indices, inner_indices)]
    )
    full_subset = M.pds_from_square_similarity(similarity)[inner_indices]
    np.testing.assert_allclose(inner, [1.0, 1.0])
    np.testing.assert_allclose(full_subset, [2.0 / 3.0, 1.0])


def test_family_aggregation_weights_three_families_equally() -> None:
    values = {
        "k562_gwps": {
            "pds": np.full(2, 2.0),
            "signed_credit": np.full(2, 2.0),
        },
        "k562_essential": {
            "pds": np.full(2, 4.0),
            "signed_credit": np.full(2, 4.0),
        },
        "rpe1": {
            "pds": np.full(2, 6.0),
            "signed_credit": np.full(2, 6.0),
        },
        "nadig_jurkat": {
            "pds": np.full(2, 9.0),
            "signed_credit": np.full(2, 9.0),
        },
    }
    summary = M.weighted_summary(values)
    assert summary["families"]["pds"] == {"k562": 3.0, "rpe1": 6.0, "jurkat": 9.0}
    assert summary["macro"]["pds"] == 6.0


def test_all_blocks_once_bootstrap_reproduces_observed_nested_result() -> None:
    with temporary_constants(FOLD_COUNT=2):
        cache, folds = synthetic_cache()
        observed = M.observed_nested(cache, folds)
        blocks = np.arange(len(folds), dtype=np.int16)
        draws = blocks[None, :]
        _, arrays = M.bootstrap_reselected(cache, folds, blocks, draws)
    expected_indices = np.asarray(
        [M.BETAS.index(record["selected_beta"]) for record in observed["folds"]],
        dtype=np.uint8,
    )
    np.testing.assert_array_equal(
        arrays["selected_beta_index_by_fold"][0], expected_indices
    )
    np.testing.assert_allclose(
        arrays["pds"][0],
        observed["oof_summary"]["macro"]["pds"],
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        arrays["signed_credit"][0],
        observed["oof_summary"]["macro"]["signed_credit"],
        rtol=0.0,
        atol=1e-15,
    )


def test_null_plans_preserve_bins_and_outer_folds() -> None:
    folds = np.repeat(np.arange(M.FOLD_COUNT, dtype=np.int8), [123, 123, 123, 122, 122])
    plans = M.make_null_plans(folds)
    for group in M.GROUPS:
        observed = M.band_labels_ascending(M.TRACE_BOUNDARIES[group])
        expected = np.bincount(observed, minlength=4)
        draws = plans[f"spectral_band_labels__{group}"]
        assert draws.shape == (M.SPECTRAL_NULL_DRAWS, M.RESPONSE_COUNT)
        for draw in draws:
            np.testing.assert_array_equal(np.bincount(draw, minlength=4), expected)
    target = plans["target_alignment_permutations"]
    indices = np.arange(M.DEV_COUNT)
    assert target.shape == (M.TARGET_NULL_DRAWS, M.DEV_COUNT)
    assert np.all(target != indices[None, :])
    assert np.all(folds[target] == folds[None, :])


def test_block_count_matrix_counts_every_draw() -> None:
    draws = np.asarray([[0, 0, 2], [1, 2, 2]], dtype=np.int16)
    counts = M.block_count_matrix(draws, 3)
    np.testing.assert_array_equal(counts, [[2, 0, 1], [0, 1, 2]])


def test_command_surface_has_no_unsafe_stage() -> None:
    source = RUN_PATH.read_text()
    assert 'choices=("prepare", "develop", "audit")' in source
    for forbidden in (
        "load_test_truth",
        "load_confirmation_truth",
        "vcc submit",
        "write_submission",
    ):
        assert forbidden not in source


def test_development_locks_before_slice_only_truth_load() -> None:
    source = inspect.getsource(M.develop)
    assert source.index("create_outcome_access_lock") < source.index(
        "load_development_truth_slice_only"
    )


def test_no_full_truth_container_hash_or_read() -> None:
    validation = inspect.getsource(M.validate_development_authorization)
    loader = inspect.getsource(M.load_development_truth_slice_only)
    assert "sha256_file(path)" not in validation
    assert "load_effect_rows" in loader
    assert "np.load(" not in loader
    develop = inspect.getsource(M.develop)
    assert "load_development_truth(0.0)" not in develop
    assert "process_tree_memory_guard" in develop


def test_actual_candidate_artifacts_are_identity_norm_and_manifest_valid() -> None:
    if not M.PREOUTCOME_DIR.exists():
        return
    M.validate_directory_manifest(M.PREOUTCOME_DIR)
    for group in M.GROUPS:
        with np.load(
            M.PREOUTCOME_DIR / f"candidates-{group}.npz", allow_pickle=False
        ) as data:
            baseline = np.asarray(data["baseline"])
            predictions = np.asarray(data["predictions"])
        assert np.array_equal(predictions[0], baseline)
        old_norm = np.linalg.norm(baseline, axis=1)
        for candidate in predictions[1:]:
            np.testing.assert_allclose(
                np.linalg.norm(candidate, axis=1),
                old_norm,
                rtol=2e-6,
                atol=2e-6,
            )


def test_all_safety_flags_are_required_false() -> None:
    required = {
        "untouched_410_test_authorized": False,
        "reserved_confirmation_expression_authorized": False,
        "count_generation_authorized": False,
        "cell_eval2_authorized": False,
        "challenge_build_authorized": False,
        "submission_authorized": False,
        "null_stage_authorized": False,
    }
    M.require_safety_flags(required)
    for key in required:
        broken = dict(required)
        broken[key] = True
        try:
            M.require_safety_flags(broken)
        except M.ContractError:
            pass
        else:
            raise AssertionError(f"unsafe flag accepted: {key}")


def main() -> int:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
