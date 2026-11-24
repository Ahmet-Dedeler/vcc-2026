#!/usr/bin/env python3
"""Synthetic tests for experiment 084; no donor data or cloud access."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


STUDY = load_module("exp084_study_test", HERE / "study.py")
RUN = STUDY.load_exp080(require_pin=True)


def test_default_fit_count_is_exact_and_conditional_refits_are_explicit() -> None:
    grid = STUDY.Grid()
    counts = STUDY.fit_counts(grid)
    assert grid.candidates == 108
    assert counts["mandatory_grid_fits"] == 43_632
    assert counts["mandatory_grid_fits_by_rank"] == {
        "4": 14_544,
        "8": 14_544,
        "16": 14_544,
    }
    assert counts["mandatory_optimizer_epochs"] == 5_235_840
    assert counts["mandatory_qr_calls"] == 87_264
    assert counts["mandatory_svd_calls"] == 43_632
    assert counts["minimum_total_fits"] == 43_632
    assert counts["maximum_if_pre_gate_fails"] == 43_834
    assert counts["maximum_if_pre_gate_passes"] == 44_137


def test_default_source_level_matmul_accounting_matches_hand_counts() -> None:
    account = STUDY.default_flop_accounting(STUDY.Dimensions(), STUDY.Grid())
    assert account["one_inner_fit_preparation"]["total"] == 4_569_154_560
    assert account["one_epoch_by_rank"] == {
        "4": 34_107_392,
        "8": 35_053_568,
        "16": 36_995_072,
    }
    grid = account["mandatory_grid"]
    assert grid["preparation_flops"] == 199_361_351_761_920
    assert grid["optimizer_flops"] == 185_271_999_528_960
    assert grid["model_fit_flops"] == 384_633_351_290_880
    assert grid["current_candidate_evaluation_flops"] == 121_815_517_593_600
    assert grid["fit_plus_evaluation_flops"] == 506_448_868_884_480
    cache = account["safe_preparation_cache"]
    assert cache["unique_inner_problems"] == 4
    assert cache["basis_svd_calls_after_cache"] == 4
    assert cache["basis_qr_calls_after_cache"] == 8


def test_exact_naive_dot_and_conventional_counts_are_distinguished() -> None:
    assert STUDY.conventional_matmul_flops(2, 3, 4) == 48
    assert STUDY.exact_naive_dot_flops(2, 3, 4) == 40


def test_cached_fit_is_byte_identical_for_real_and_null_labels() -> None:
    rng = np.random.default_rng(840)
    n, d, c = 40, 31, 23
    query = rng.normal(size=(n, d))
    truth = query + rng.normal(scale=0.4, size=(n, d))
    truths = {"teacher": truth}
    rows = np.arange(n, dtype=np.int64)
    correction = np.arange(c, dtype=np.int64)
    basis_seed = 841
    prepared = STUDY.prepare_fit_problem(
        RUN,
        query,
        truths,
        correction,
        rows,
        basis_dim=12,
        basis_seed=basis_seed,
    )
    for rank in (4, 8):
        candidate = RUN.Candidate(rank, 0.10, 0.07, 1e-3)
        for permuted in (False, True):
            permutation_seed = 842 if permuted else None
            baseline = RUN.fit_adapter(
                query,
                truths,
                correction,
                candidate,
                rows,
                basis_dim=12,
                epochs=4,
                learning_rate=0.03,
                gradient_clip=5.0,
                seed=843 + rank,
                basis_seed=basis_seed,
                permutation_seed=permutation_seed,
                permute_labels=permuted,
            )
            cached = STUDY.fit_prepared(
                RUN,
                prepared,
                candidate,
                epochs=4,
                learning_rate=0.03,
                gradient_clip=5.0,
                seed=843 + rank,
                permutation_seed=permutation_seed,
                permute_labels=permuted,
            )
            assert np.array_equal(baseline.basis, cached.basis)
            assert np.array_equal(baseline.u, cached.u)
            assert np.array_equal(baseline.v, cached.v)
            assert baseline.fit_audit == cached.fit_audit


def test_cached_problem_is_candidate_invariant_but_models_are_not_collapsed() -> None:
    rng = np.random.default_rng(844)
    n, d, c = 30, 26, 19
    query = rng.normal(size=(n, d))
    truth = rng.normal(size=(n, d))
    rows = np.arange(n, dtype=np.int64)
    correction = np.arange(c, dtype=np.int64)
    prepared = STUDY.prepare_fit_problem(
        RUN,
        query,
        {"teacher": truth},
        correction,
        rows,
        basis_dim=10,
        basis_seed=845,
    )
    first = RUN.Candidate(4, 0.02, 0.03, 1e-4)
    second = RUN.Candidate(4, 0.20, 0.15, 1e-2)
    first_model = STUDY.fit_prepared(
        RUN,
        prepared,
        first,
        epochs=3,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=846,
    )
    second_model = STUDY.fit_prepared(
        RUN,
        prepared,
        second,
        epochs=3,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=846,
    )
    assert first_model.basis is second_model.basis
    assert not np.array_equal(first_model.u, second_model.u)
    assert not np.array_equal(first_model.v, second_model.v)


def test_latent_evaluation_preserves_pds_on_representative_data() -> None:
    rng = np.random.default_rng(847)
    n, d, c = 45, 29, 21
    query = rng.normal(size=(n, d))
    truth = query + rng.normal(scale=0.5, size=(n, d))
    rows = np.arange(15, dtype=np.int64)
    correction = np.arange(c, dtype=np.int64)
    candidate = RUN.Candidate(4, 0.10, 0.07, 1e-3)
    model = RUN.fit_adapter(
        query,
        {"teacher": truth},
        correction,
        candidate,
        np.arange(n, dtype=np.int64),
        basis_dim=12,
        epochs=4,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=848,
        basis_seed=849,
    )
    current = RUN.evaluate_adapter(
        query,
        truth,
        correction,
        model,
        rows,
        rows,
    )
    latent = STUDY.latent_evaluation(
        RUN,
        query,
        truth,
        correction,
        model,
        rows,
        rows,
    )
    assert np.array_equal(
        current["per_target_baseline"], latent["per_target_baseline"]
    )
    assert np.array_equal(
        current["per_target_adapted"], latent["per_target_adapted"]
    )
    np.testing.assert_allclose(
        current["angles_degrees"], latent["angles_degrees"], atol=2e-7
    )


def test_memory_accounting_cache_payload_is_small() -> None:
    memory = STUDY.array_memory_accounting(STUDY.Dimensions())
    assert (
        memory["named_array_bytes"]["persistent_3_queries_plus_2_truths"]
        == 161_088_000
    )
    assert memory["cache_payload_mib"]["one_inner_problem"] < 4.1
    assert (
        memory["cache_payload_mib"]["all_4_inner_plus_2_pretest_plus_3_final"]
        < 58
    )


def test_report_output_is_confined_to_experiment_084() -> None:
    parser = STUDY.argparse.ArgumentParser()
    # The actual confinement is enforced in main; this guards the invariant
    # used by the implementation and documents the intended root.
    assert STUDY.HERE.name == "084-listwise-runtime-study"
    assert STUDY.EXP080.name == "080-listwise-metric-adapter"
    assert STUDY.HERE != STUDY.EXP080
    assert parser is not None


def main() -> int:
    tests = [
        test_default_fit_count_is_exact_and_conditional_refits_are_explicit,
        test_default_source_level_matmul_accounting_matches_hand_counts,
        test_exact_naive_dot_and_conventional_counts_are_distinguished,
        test_cached_fit_is_byte_identical_for_real_and_null_labels,
        test_cached_problem_is_candidate_invariant_but_models_are_not_collapsed,
        test_latent_evaluation_preserves_pds_on_representative_data,
        test_memory_accounting_cache_payload_is_small,
        test_report_output_is_confined_to_experiment_084,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed; peak RSS={STUDY.peak_rss_gib():.3f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
