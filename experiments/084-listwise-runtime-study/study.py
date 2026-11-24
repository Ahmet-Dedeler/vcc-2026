#!/usr/bin/env python3
"""Read-only runtime accounting and bounded benchmark for experiment 080.

This module never edits experiment 080. It imports the pinned implementation,
builds a semantics-equivalent prepared-fit path in this directory, and writes
only the requested JSON report.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import resource
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

# Keep the representative run bounded and reproducible. These must be set
# before NumPy is imported.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXP080 = REPO / "experiments" / "080-listwise-metric-adapter"
RUN_PATH = EXP080 / "run.py"
PINNED_RUN_SHA256 = "5680481fcedac491f37a6f0e3f31a3d33869885c6fdc44bd761532286fa7fbae"
SNAPSHOT_UTC = "2026-08-23T07:39:00Z"


@dataclass(frozen=True)
class Dimensions:
    training_targets: int = 600
    fit_targets: int = 360
    tune_targets: int = 120
    test_targets: int = 120
    master_genes: int = 7312
    train_score_genes: int = 6712
    correction_genes: int = 5582
    basis_dim: int = 64
    randomized_basis_oversample: int = 8
    teachers_inner: int = 1
    teachers_refit: int = 2


@dataclass(frozen=True)
class Grid:
    ranks: tuple[int, ...] = (4, 8, 16)
    betas: tuple[float, ...] = (0.02, 0.05, 0.10, 0.20)
    temperatures: tuple[float, ...] = (0.03, 0.07, 0.15)
    regularizations: tuple[float, ...] = (1e-4, 1e-3, 1e-2)
    epochs: int = 120
    null_replicates: int = 100
    development_folds: int = 2
    held_teacher_roles: int = 2
    final_folds: int = 3

    @property
    def candidates(self) -> int:
        return (
            len(self.ranks)
            * len(self.betas)
            * len(self.temperatures)
            * len(self.regularizations)
        )

    @property
    def candidates_per_rank(self) -> int:
        return (
            len(self.betas)
            * len(self.temperatures)
            * len(self.regularizations)
        )

    @property
    def selection_repetitions(self) -> int:
        return 1 + self.null_replicates

    @property
    def cells_per_selection(self) -> int:
        return self.development_folds * self.held_teacher_roles

    @property
    def mandatory_grid_fits(self) -> int:
        return (
            self.selection_repetitions
            * self.cells_per_selection
            * self.candidates
        )


@dataclass(frozen=True)
class PreparedPair:
    teacher: str
    base_scores: np.ndarray
    truth_projection: np.ndarray
    raw_truth_norms: np.ndarray


@dataclass(frozen=True)
class PreparedFit:
    basis: np.ndarray
    query_projection: np.ndarray
    pairs: tuple[PreparedPair, ...]
    rows: int
    basis_dim: int


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def load_exp080(*, require_pin: bool = True):
    observed = sha256_file(RUN_PATH)
    if require_pin and observed != PINNED_RUN_SHA256:
        raise RuntimeError(
            "experiment 080 changed after the runtime-study snapshot: "
            f"{observed} != {PINNED_RUN_SHA256}"
        )
    spec = importlib.util.spec_from_file_location("exp080_runtime_study", RUN_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parser_defaults(run) -> dict[str, object]:
    wanted = {
        "ranks",
        "betas",
        "temperatures",
        "regularizations",
        "basis_dim",
        "epochs",
        "null_replicates",
        "bootstrap_replicates",
        "learning_rate",
        "gradient_clip",
    }
    return {
        action.dest: action.default
        for action in run.build_parser()._actions
        if action.dest in wanted
    }


def conventional_matmul_flops(rows: int, inner: int, columns: int) -> int:
    """Conventional dense-GEMM count, where one multiply and one add are 2 FLOPs."""
    return 2 * rows * inner * columns


def exact_naive_dot_flops(rows: int, inner: int, columns: int) -> int:
    """Exact scalar multiply/add count for naive dot products."""
    return rows * columns * (2 * inner - 1)


def fit_preparation_matmul_flops(
    *,
    rows: int,
    score_genes: int,
    correction_genes: int,
    basis_dim: int,
    teachers: int,
    oversample: int = 8,
) -> dict[str, int]:
    width = basis_dim + oversample
    basis_rows = rows * (1 + teachers)
    calls = {
        "basis_values_times_omega": conventional_matmul_flops(
            basis_rows, correction_genes, width
        ),
        "basis_values_transpose_times_q": conventional_matmul_flops(
            correction_genes, basis_rows, width
        ),
        "basis_values_times_power_projection": conventional_matmul_flops(
            basis_rows, correction_genes, width
        ),
        "basis_q_transpose_times_values": conventional_matmul_flops(
            width, basis_rows, correction_genes
        ),
        "query_projection": conventional_matmul_flops(
            rows, correction_genes, basis_dim
        ),
        "truth_projections": teachers
        * conventional_matmul_flops(rows, correction_genes, basis_dim),
        "base_scores": teachers
        * conventional_matmul_flops(rows, score_genes, rows),
    }
    return {**calls, "total": sum(calls.values())}


def epoch_matmul_flops(
    *, rows: int, basis_dim: int, rank: int, teachers: int
) -> dict[str, int]:
    # Per teacher, _softmax_cross_entropy_gradient has five n-k-r products
    # and two n-n-k products. The regularizer has four k-r-r products.
    calls = {
        "pair_low_rank_products": teachers * 10 * rows * basis_dim * rank,
        "pair_gallery_products": teachers * 4 * rows * rows * basis_dim,
        "regularizer_products": 8 * basis_dim * rank * rank,
    }
    return {**calls, "total": sum(calls.values())}


def current_selection_evaluation_matmul_flops(
    *,
    query_rows_total: int,
    evaluated_rows: int,
    gallery_rows: int,
    score_genes: int,
    correction_genes: int,
    basis_dim: int,
    rank: int,
) -> int:
    # transform_query runs on every query row, then both baseline and adapted
    # score matrices are formed on the selected rows.
    return (
        2 * conventional_matmul_flops(
            query_rows_total, correction_genes, basis_dim
        )
        + 4 * query_rows_total * basis_dim * rank
        + 2
        * conventional_matmul_flops(evaluated_rows, score_genes, gallery_rows)
    )


def fit_counts(grid: Grid, selected_nonidentity: int | None = None) -> dict[str, object]:
    if selected_nonidentity is not None and not (
        0 <= selected_nonidentity <= grid.selection_repetitions
    ):
        raise ValueError("selected_nonidentity is outside [0, selection_repetitions]")
    mandatory = grid.mandatory_grid_fits
    formula = (
        f"{mandatory} + 2*S + 3*I(pre_gate_passed)*S, "
        f"where 0<=S<={grid.selection_repetitions}"
    )
    result: dict[str, object] = {
        "candidate_count_excluding_identity": grid.candidates,
        "candidate_count_including_identity": grid.candidates + 1,
        "selection_repetitions_real_plus_null": grid.selection_repetitions,
        "selection_cells_per_repetition": grid.cells_per_selection,
        "mandatory_grid_fits": mandatory,
        "mandatory_grid_fits_by_rank": {
            str(rank): (
                grid.selection_repetitions
                * grid.cells_per_selection
                * grid.candidates_per_rank
            )
            for rank in grid.ranks
        },
        "mandatory_optimizer_epochs": mandatory * grid.epochs,
        "mandatory_qr_calls": mandatory * 2,
        "mandatory_svd_calls": mandatory,
        "fit_count_formula": formula,
        "minimum_total_fits": mandatory,
        "maximum_if_pre_gate_fails": mandatory
        + grid.development_folds * grid.selection_repetitions,
        "maximum_if_pre_gate_passes": mandatory
        + (grid.development_folds + grid.final_folds)
        * grid.selection_repetitions,
        "null_share_of_grid_fits": grid.null_replicates
        / grid.selection_repetitions,
    }
    if selected_nonidentity is not None:
        result["for_selected_nonidentity_count"] = {
            "S": selected_nonidentity,
            "pretest_refits": grid.development_folds * selected_nonidentity,
            "total_if_pre_gate_fails": mandatory
            + grid.development_folds * selected_nonidentity,
            "final_refits_if_pre_gate_passes": grid.final_folds
            * selected_nonidentity,
            "total_if_pre_gate_passes": mandatory
            + (grid.development_folds + grid.final_folds)
            * selected_nonidentity,
        }
    return result


def default_flop_accounting(dim: Dimensions, grid: Grid) -> dict[str, object]:
    prep = fit_preparation_matmul_flops(
        rows=dim.fit_targets,
        score_genes=dim.train_score_genes,
        correction_genes=dim.correction_genes,
        basis_dim=dim.basis_dim,
        teachers=dim.teachers_inner,
        oversample=dim.randomized_basis_oversample,
    )
    per_rank_epoch = {
        rank: epoch_matmul_flops(
            rows=dim.fit_targets,
            basis_dim=dim.basis_dim,
            rank=rank,
            teachers=dim.teachers_inner,
        )["total"]
        for rank in grid.ranks
    }
    fits_per_rank = (
        grid.selection_repetitions
        * grid.cells_per_selection
        * grid.candidates_per_rank
    )
    preparation_total = prep["total"] * grid.mandatory_grid_fits
    optimization_total = sum(
        fits_per_rank * grid.epochs * per_rank_epoch[rank] for rank in grid.ranks
    )
    evaluation_per_rank = {
        rank: current_selection_evaluation_matmul_flops(
            query_rows_total=dim.training_targets,
            evaluated_rows=dim.tune_targets,
            gallery_rows=dim.training_targets,
            score_genes=dim.train_score_genes,
            correction_genes=dim.correction_genes,
            basis_dim=dim.basis_dim,
            rank=rank,
        )
        for rank in grid.ranks
    }
    evaluation_total = sum(
        fits_per_rank * evaluation_per_rank[rank] for rank in grid.ranks
    )
    cached_unique_inner_problems = grid.cells_per_selection
    cached_preparation_total = prep["total"] * cached_unique_inner_problems
    return {
        "convention": (
            "Exact count of source-level dense @ operations using the standard "
            "2*m*n*k convention. QR, SVD, exp, sqrt, divisions, comparisons, "
            "copies, and BLAS implementation details are reported separately."
        ),
        "one_inner_fit_preparation": prep,
        "one_epoch_by_rank": {str(k): v for k, v in per_rank_epoch.items()},
        "one_120_epoch_fit_by_rank": {
            str(rank): prep["total"] + grid.epochs * per_rank_epoch[rank]
            for rank in grid.ranks
        },
        "mandatory_grid": {
            "preparation_flops": preparation_total,
            "optimizer_flops": optimization_total,
            "model_fit_flops": preparation_total + optimization_total,
            "current_candidate_evaluation_flops": evaluation_total,
            "fit_plus_evaluation_flops": (
                preparation_total + optimization_total + evaluation_total
            ),
            "current_basis_qr_calls": grid.mandatory_grid_fits * 2,
            "current_basis_svd_calls": grid.mandatory_grid_fits,
        },
        "safe_preparation_cache": {
            "unique_inner_problems": cached_unique_inner_problems,
            "preparation_flops_after_cache": cached_preparation_total,
            "preparation_flops_avoided": preparation_total
            - cached_preparation_total,
            "basis_qr_calls_after_cache": cached_unique_inner_problems * 2,
            "basis_svd_calls_after_cache": cached_unique_inner_problems,
            "optimizer_flops_unchanged": optimization_total,
            "selection_and_null_repetitions_unchanged": True,
        },
        "current_selection_evaluation_per_candidate_by_rank": {
            str(k): v for k, v in evaluation_per_rank.items()
        },
    }


def array_memory_accounting(dim: Dimensions) -> dict[str, object]:
    f64 = 8
    n = dim.fit_targets
    d = dim.train_score_genes
    c = dim.correction_genes
    k = dim.basis_dim
    w = k + dim.randomized_basis_oversample
    named = {
        "persistent_3_queries_plus_2_truths": 5
        * dim.training_targets
        * d
        * f64,
        "one_inner_normalized_query_plus_truth": 2 * n * d * f64,
        "one_inner_basis_source": 2 * n * c * f64,
        "one_basis_omega": c * w * f64,
        "one_basis_projected": 2 * n * w * f64,
        "one_basis_c_by_w_temporary": c * w * f64,
        "one_basis": c * k * f64,
        "one_base_score_matrix": n * n * f64,
        "two_projection_matrices": 2 * n * k * f64,
    }
    one_cached_inner = (
        c * k + n * n + 2 * n * k
    ) * f64
    pre_n = dim.fit_targets + dim.tune_targets
    one_cached_pretest = (
        c * k
        + dim.teachers_refit * pre_n * pre_n
        + (1 + dim.teachers_refit) * pre_n * k
    ) * f64
    final_n = dim.training_targets
    one_cached_final = (
        c * k
        + dim.teachers_refit * final_n * final_n
        + (1 + dim.teachers_refit) * final_n * k
    ) * f64
    all_problem_payloads = (
        4 * one_cached_inner + 2 * one_cached_pretest + 3 * one_cached_final
    )
    return {
        "dtype": "float64",
        "named_array_bytes": named,
        "named_array_mib": {
            key: value / 2**20 for key, value in named.items()
        },
        "named_preparation_lower_bound_bytes": sum(named.values()),
        "named_preparation_lower_bound_mib": sum(named.values()) / 2**20,
        "cache_payload_bytes": {
            "one_inner_problem": one_cached_inner,
            "one_pretest_problem": one_cached_pretest,
            "one_final_problem": one_cached_final,
            "all_4_inner_plus_2_pretest_plus_3_final": all_problem_payloads,
        },
        "cache_payload_mib": {
            "one_inner_problem": one_cached_inner / 2**20,
            "one_pretest_problem": one_cached_pretest / 2**20,
            "one_final_problem": one_cached_final / 2**20,
            "all_4_inner_plus_2_pretest_plus_3_final": all_problem_payloads / 2**20,
        },
        "scope_note": (
            "Exact NumPy payload bytes for named arrays only. Process RSS also "
            "contains Python, NumPy, allocator, and BLAS/LAPACK workspaces and "
            "is measured by the bounded benchmark."
        ),
    }


def prepare_fit_problem(
    run,
    query: np.ndarray,
    truths: Mapping[str, np.ndarray],
    correction_in_score: np.ndarray,
    rows: np.ndarray,
    *,
    basis_dim: int,
    basis_seed: int,
) -> PreparedFit:
    """Compute the candidate- and null-invariant portion of fit_adapter."""
    q = run.row_normalize(np.asarray(query[rows], dtype=np.float64))
    normalized_truths = {
        name: run.row_normalize(np.asarray(values[rows], dtype=np.float64))
        for name, values in truths.items()
    }
    correction = np.asarray(correction_in_score, dtype=np.int64)
    basis_rank = min(
        int(basis_dim),
        len(correction),
        len(rows) * (1 + len(normalized_truths)),
    )
    basis_source = np.vstack(
        [
            q[:, correction],
            *(truth[:, correction] for truth in normalized_truths.values()),
        ]
    )
    basis = run.randomized_right_basis(basis_source, basis_rank, basis_seed)
    query_projection = q[:, correction] @ basis
    pairs = tuple(
        PreparedPair(
            teacher=teacher,
            base_scores=q @ truth.T,
            truth_projection=truth[:, correction] @ basis,
            raw_truth_norms=np.linalg.norm(
                np.asarray(truths[teacher][rows], dtype=np.float64), axis=1
            ),
        )
        for teacher, truth in normalized_truths.items()
    )
    return PreparedFit(
        basis=basis,
        query_projection=query_projection,
        pairs=pairs,
        rows=len(rows),
        basis_dim=basis_rank,
    )


def labels_from_norms(run, norms: np.ndarray, *, seed: int) -> np.ndarray:
    order = np.argsort(np.asarray(norms, dtype=np.float64), kind="stable")
    labels = np.arange(len(norms), dtype=np.int64)
    rng = np.random.default_rng(seed)
    for bucket in np.array_split(order, min(10, len(order))):
        labels[bucket] = rng.permutation(bucket)
    if sorted(labels.tolist()) != list(range(len(labels))):
        raise AssertionError("cached norm-stratified permutation is not bijective")
    return labels


def fit_prepared(
    run,
    prepared: PreparedFit,
    candidate,
    *,
    epochs: int,
    learning_rate: float,
    gradient_clip: float,
    seed: int,
    permutation_seed: int | None = None,
    permute_labels: bool = False,
):
    """Run the unchanged optimizer over cached, immutable training pairs."""
    pairs = []
    permutation_hashes: dict[str, str] = {}
    for item in prepared.pairs:
        labels = (
            labels_from_norms(
                run,
                item.raw_truth_norms,
                seed=run.stable_seed(
                    seed if permutation_seed is None else permutation_seed,
                    item.teacher,
                    "permutation",
                ),
            )
            if permute_labels
            else np.arange(prepared.rows, dtype=np.int64)
        )
        permutation_hashes[item.teacher] = run.sha256_array(labels)
        pairs.append(
            run.TrainingPair(
                base_scores=item.base_scores,
                query_projection=prepared.query_projection,
                truth_projection=item.truth_projection,
                labels=labels,
            )
        )

    rng = np.random.default_rng(seed)
    initial = rng.standard_normal((prepared.basis_dim, candidate.rank))
    v, _ = np.linalg.qr(initial, mode="reduced")
    u = np.zeros((prepared.basis_dim, candidate.rank), dtype=np.float64)
    first_u = np.zeros_like(u)
    second_u = np.zeros_like(u)
    first_v = np.zeros_like(v)
    second_v = np.zeros_like(v)
    beta1, beta2 = 0.9, 0.999
    losses: list[float] = []
    for step in range(1, epochs + 1):
        loss = 0.0
        gradient_u = np.zeros_like(u)
        gradient_v = np.zeros_like(v)
        for pair in pairs:
            pair_loss, pair_gu, pair_gv = run._softmax_cross_entropy_gradient(
                pair,
                u,
                v,
                beta=candidate.beta,
                temperature=candidate.temperature,
            )
            loss += pair_loss
            gradient_u += pair_gu
            gradient_v += pair_gv
        loss /= len(pairs)
        gradient_u /= len(pairs)
        gradient_v /= len(pairs)

        gram_u = u.T @ u
        gram_v = v.T @ v
        residual_frobenius_sq = float(np.sum(gram_u * gram_v))
        loss += 0.5 * candidate.regularization * residual_frobenius_sq
        gradient_u += candidate.regularization * (u @ gram_v)
        gradient_v += candidate.regularization * (v @ gram_u)
        factor_regularization = candidate.regularization * 1e-3
        loss += 0.5 * factor_regularization * (
            float(np.sum(u * u)) + float(np.sum(v * v))
        )
        gradient_u += factor_regularization * u
        gradient_v += factor_regularization * v

        gradient_norm = math.sqrt(
            float(np.sum(gradient_u * gradient_u) + np.sum(gradient_v * gradient_v))
        )
        if gradient_norm > gradient_clip:
            scale = gradient_clip / gradient_norm
            gradient_u *= scale
            gradient_v *= scale
        first_u = beta1 * first_u + (1 - beta1) * gradient_u
        second_u = beta2 * second_u + (1 - beta2) * (gradient_u * gradient_u)
        first_v = beta1 * first_v + (1 - beta1) * gradient_v
        second_v = beta2 * second_v + (1 - beta2) * (gradient_v * gradient_v)
        correction1 = 1 - beta1**step
        correction2 = 1 - beta2**step
        u -= (
            learning_rate
            * (first_u / correction1)
            / (np.sqrt(second_u / correction2) + 1e-8)
        )
        v -= (
            learning_rate
            * (first_v / correction1)
            / (np.sqrt(second_v / correction2) + 1e-8)
        )
        if (
            not np.isfinite(loss)
            or not np.isfinite(u).all()
            or not np.isfinite(v).all()
        ):
            raise FloatingPointError(f"{candidate.key}: non-finite optimizer state")
        losses.append(loss)

    product_frobenius = math.sqrt(float(np.sum((u.T @ u) * (v.T @ v))))
    return run.Adapter(
        candidate=candidate,
        basis=prepared.basis,
        u=u,
        v=v,
        fit_audit={
            "rows": int(prepared.rows),
            "teachers": [pair.teacher for pair in prepared.pairs],
            "basis_dim": int(prepared.basis_dim),
            "epochs": int(epochs),
            "learning_rate": float(learning_rate),
            "gradient_clip": float(gradient_clip),
            "first_loss": float(losses[0]),
            "last_loss": float(losses[-1]),
            "minimum_loss": float(min(losses)),
            "residual_product_frobenius": product_frobenius,
            "labels_permuted_within_norm_deciles": bool(permute_labels),
            "label_permutation_sha256": permutation_hashes,
        },
    )


def latent_evaluation(
    run,
    query: np.ndarray,
    truth: np.ndarray,
    correction_in_score: np.ndarray,
    adapter,
    query_rows: np.ndarray,
    positive_columns: np.ndarray,
) -> dict[str, np.ndarray | float]:
    """Algebraically exact fast path; not recommended without real tie audit."""
    rows = np.asarray(query_rows, dtype=np.int64)
    q = run.row_normalize(np.asarray(query[rows], dtype=np.float64))
    t = run.row_normalize(np.asarray(truth, dtype=np.float64))
    baseline = q @ t.T
    x = q[:, correction_in_score] @ adapter.basis
    truth_projection = t[:, correction_in_score] @ adapter.basis
    residual = (x @ adapter.v) @ adapter.u.T
    cross = np.sum(x * residual, axis=1)
    residual_norm_sq = np.sum(residual * residual, axis=1)
    denominator = np.sqrt(
        np.maximum(
            1.0
            + 2.0 * adapter.candidate.beta * cross
            + adapter.candidate.beta**2 * residual_norm_sq,
            1e-18,
        )
    )
    adapted_scores = (
        baseline
        + adapter.candidate.beta * (residual @ truth_projection.T)
    ) / denominator[:, None]
    baseline_pds = run.pds_per_query_from_scores(baseline, positive_columns)
    adapted_pds = run.pds_per_query_from_scores(adapted_scores, positive_columns)
    cosine = np.clip(
        (1.0 + adapter.candidate.beta * cross) / denominator, -1.0, 1.0
    )
    return {
        "per_target_baseline": baseline_pds,
        "per_target_adapted": adapted_pds,
        "angles_degrees": np.degrees(np.arccos(cosine)),
        "delta": float(np.mean(adapted_pds - baseline_pds)),
    }


def timed(callable_):
    start = time.perf_counter()
    value = callable_()
    return value, time.perf_counter() - start


def run_benchmark(
    run,
    *,
    dim: Dimensions,
    epochs: int,
    seed: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    n = dim.fit_targets
    d = dim.train_score_genes
    c = dim.correction_genes
    query = rng.normal(size=(n, d))
    truth = query + rng.normal(scale=0.35, size=(n, d))
    truths = {"teacher": truth}
    rows = np.arange(n, dtype=np.int64)
    correction = np.arange(c, dtype=np.int64)
    basis_seed = run.stable_seed(8012026, "inner-basis", "outer", "held")
    candidates = (
        run.Candidate(4, 0.02, 0.03, 1e-4),
        run.Candidate(8, 0.10, 0.07, 1e-3),
        run.Candidate(16, 0.20, 0.15, 1e-2),
    )

    prepared, prepare_seconds = timed(
        lambda: prepare_fit_problem(
            run,
            query,
            truths,
            correction,
            rows,
            basis_dim=dim.basis_dim,
            basis_seed=basis_seed,
        )
    )
    baseline_seconds: dict[str, float] = {}
    cached_seconds: dict[str, float] = {}
    optimizer_seconds_per_epoch: dict[str, float] = {}
    equivalence: dict[str, object] = {}
    cached_models = {}
    for candidate in candidates:
        optimizer_seed = run.stable_seed(
            8012026, "inner-init", "outer", "held", candidate.rank
        )
        baseline, baseline_elapsed = timed(
            lambda candidate=candidate, optimizer_seed=optimizer_seed: run.fit_adapter(
                query,
                truths,
                correction,
                candidate,
                rows,
                basis_dim=dim.basis_dim,
                epochs=epochs,
                learning_rate=0.03,
                gradient_clip=5.0,
                seed=optimizer_seed,
                basis_seed=basis_seed,
            )
        )
        cached, cached_elapsed = timed(
            lambda candidate=candidate, optimizer_seed=optimizer_seed: fit_prepared(
                run,
                prepared,
                candidate,
                epochs=epochs,
                learning_rate=0.03,
                gradient_clip=5.0,
                seed=optimizer_seed,
            )
        )
        key = str(candidate.rank)
        baseline_seconds[key] = baseline_elapsed
        cached_seconds[key] = cached_elapsed
        optimizer_seconds_per_epoch[key] = cached_elapsed / epochs
        equivalence[key] = {
            "basis_array_equal": bool(np.array_equal(baseline.basis, cached.basis)),
            "u_array_equal": bool(np.array_equal(baseline.u, cached.u)),
            "v_array_equal": bool(np.array_equal(baseline.v, cached.v)),
            "fit_audit_equal": baseline.fit_audit == cached.fit_audit,
            "max_abs_u": float(np.max(np.abs(baseline.u - cached.u))),
            "max_abs_v": float(np.max(np.abs(baseline.v - cached.v))),
        }
        cached_models[candidate.rank] = cached

    # Verify the cached-label path too. This is one representative null, while
    # the accounting preserves all 100 complete selection repetitions.
    candidate = candidates[1]
    optimizer_seed = run.stable_seed(
        8012026, "inner-init", "outer", "held", candidate.rank
    )
    permutation_seed = run.stable_seed(8022026, "inner", 17, "outer", "held")
    null_baseline = run.fit_adapter(
        query,
        truths,
        correction,
        candidate,
        rows,
        basis_dim=dim.basis_dim,
        epochs=epochs,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=optimizer_seed,
        basis_seed=basis_seed,
        permutation_seed=permutation_seed,
        permute_labels=True,
    )
    null_cached = fit_prepared(
        run,
        prepared,
        candidate,
        epochs=epochs,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=optimizer_seed,
        permutation_seed=permutation_seed,
        permute_labels=True,
    )
    null_equivalence = {
        "basis_array_equal": bool(
            np.array_equal(null_baseline.basis, null_cached.basis)
        ),
        "u_array_equal": bool(np.array_equal(null_baseline.u, null_cached.u)),
        "v_array_equal": bool(np.array_equal(null_baseline.v, null_cached.v)),
        "fit_audit_equal": null_baseline.fit_audit == null_cached.fit_audit,
        "label_hashes_equal": (
            null_baseline.fit_audit["label_permutation_sha256"]
            == null_cached.fit_audit["label_permutation_sha256"]
        ),
    }

    # One full-shape selection evaluation. The faster latent form is measured
    # and audited, but intentionally not the recommended first change.
    query_eval = rng.normal(size=(dim.training_targets, d))
    truth_eval = query_eval + rng.normal(
        scale=0.35, size=(dim.training_targets, d)
    )
    tune_rows = np.arange(dim.tune_targets, dtype=np.int64)
    positives = tune_rows.copy()
    eval_model = cached_models[16]
    current_eval, current_eval_seconds = timed(
        lambda: run.evaluate_adapter(
            query_eval,
            truth_eval,
            correction,
            eval_model,
            tune_rows,
            positives,
        )
    )
    latent_eval, latent_eval_seconds = timed(
        lambda: latent_evaluation(
            run,
            query_eval,
            truth_eval,
            correction,
            eval_model,
            tune_rows,
            positives,
        )
    )
    eval_equivalence = {
        "baseline_pds_array_equal": bool(
            np.array_equal(
                current_eval["per_target_baseline"],
                latent_eval["per_target_baseline"],
            )
        ),
        "adapted_pds_array_equal": bool(
            np.array_equal(
                current_eval["per_target_adapted"],
                latent_eval["per_target_adapted"],
            )
        ),
        "max_abs_angle_degrees": float(
            np.max(
                np.abs(
                    current_eval["angles_degrees"]
                    - latent_eval["angles_degrees"]
                )
            )
        ),
        "delta_abs_difference": float(
            abs(float(current_eval["delta"]) - float(latent_eval["delta"]))
        ),
    }

    grid = Grid()
    model_current_seconds = (
        grid.mandatory_grid_fits * prepare_seconds
        + sum(
            (
                grid.selection_repetitions
                * grid.cells_per_selection
                * grid.candidates_per_rank
            )
            * grid.epochs
            * optimizer_seconds_per_epoch[str(rank)]
            for rank in grid.ranks
        )
    )
    model_cached_seconds = (
        grid.cells_per_selection * prepare_seconds
        + sum(
            (
                grid.selection_repetitions
                * grid.cells_per_selection
                * grid.candidates_per_rank
            )
            * grid.epochs
            * optimizer_seconds_per_epoch[str(rank)]
            for rank in grid.ranks
        )
    )
    unchanged_eval_seconds = grid.mandatory_grid_fits * current_eval_seconds
    return {
        "kind": "synthetic full-dimension, one-inner-cell representative",
        "shape": {
            "fit_rows": n,
            "score_genes": d,
            "correction_genes": c,
            "basis_dim": dim.basis_dim,
            "teachers": 1,
            "timed_epochs_per_fit": epochs,
        },
        "threads": {
            name: os.environ.get(name)
            for name in (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
        "prepare_once_seconds": prepare_seconds,
        "baseline_fit_seconds_by_rank": baseline_seconds,
        "cached_fit_seconds_by_rank": cached_seconds,
        "cached_optimizer_seconds_per_epoch_by_rank": optimizer_seconds_per_epoch,
        "fit_equivalence": equivalence,
        "null_fit_equivalence": null_equivalence,
        "selection_evaluation": {
            "current_seconds": current_eval_seconds,
            "latent_seconds": latent_eval_seconds,
            "speedup": current_eval_seconds / latent_eval_seconds,
            "equivalence": eval_equivalence,
            "recommendation_status": (
                "defer: algebraically equivalent and PDS-identical here, but "
                "requires exhaustive real-data rank/tie equivalence before use"
            ),
        },
        "local_single_thread_extrapolation": {
            "warning": (
                "A ratio-oriented extrapolation from synthetic local timings, "
                "not an AWS wall-clock forecast. It excludes staging, I/O, "
                "pretest/final conditional refits, and bootstrap."
            ),
            "current_grid_model_seconds": model_current_seconds,
            "cached_grid_model_seconds": model_cached_seconds,
            "safe_cache_model_speedup": model_current_seconds
            / model_cached_seconds,
            "current_grid_evaluation_seconds_if_unchanged": unchanged_eval_seconds,
            "current_fit_plus_eval_seconds": model_current_seconds
            + unchanged_eval_seconds,
            "safe_cache_fit_plus_unchanged_eval_seconds": model_cached_seconds
            + unchanged_eval_seconds,
            "safe_cache_end_to_end_grid_speedup": (
                (model_current_seconds + unchanged_eval_seconds)
                / (model_cached_seconds + unchanged_eval_seconds)
            ),
        },
        "peak_rss_gib": peak_rss_gib(),
    }


def build_report(*, benchmark_epochs: int, benchmark_seed: int) -> dict[str, object]:
    run = load_exp080(require_pin=True)
    dim = Dimensions()
    grid = Grid()
    defaults = parser_defaults(run)
    expected_defaults = {
        "ranks": list(grid.ranks),
        "betas": list(grid.betas),
        "temperatures": list(grid.temperatures),
        "regularizations": list(grid.regularizations),
        "basis_dim": dim.basis_dim,
        "epochs": grid.epochs,
        "null_replicates": grid.null_replicates,
        "bootstrap_replicates": 20_000,
        "learning_rate": 0.03,
        "gradient_clip": 5.0,
    }
    if defaults != expected_defaults:
        raise RuntimeError(f"080 defaults drifted: {defaults} != {expected_defaults}")
    benchmark = run_benchmark(
        run,
        dim=dim,
        epochs=benchmark_epochs,
        seed=benchmark_seed,
    )
    if benchmark["peak_rss_gib"] >= 1.0:
        raise MemoryError(
            f"benchmark peak RSS {benchmark['peak_rss_gib']:.3f} GiB is not <1 GiB"
        )
    return {
        "study": "084-listwise-runtime-study",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "read_only_target": str(EXP080),
        "snapshot": {
            "snapshot_utc": SNAPSHOT_UTC,
            "run_py": str(RUN_PATH),
            "run_py_sha256": sha256_file(RUN_PATH),
            "pinned_run_py_sha256": PINNED_RUN_SHA256,
            "line_anchors": {
                "randomized_right_basis": "run.py:269-294",
                "fit_adapter": "run.py:357-526",
                "complete_candidate_selection": "run.py:1133-1226",
                "evaluate_selected_pretest": "run.py:1229-1298",
                "fit_and_freeze_final_models": "run.py:1345-1390",
                "null_loop": "run.py:1515-1543",
                "defaults": "run.py:2036-2072",
            },
        },
        "evidence": {
            "dimensions": asdict(dim),
            "dimension_sources": {
                "training_targets": (
                    "artifact_manifest truth-cache meta.json and train_perts.txt"
                ),
                "master_genes": (
                    "S3-pinned master_genes.txt SHA-256 "
                    "dedeefb0d286fb0e715ef2a5175320067408b0445f591e39acd87a6560db168f"
                ),
                "correction_mask": (
                    "S3-pinned train_mask.npy SHA-256 "
                    "0762c03532e38590d50335f695065dbf13c09d928b6be599e831bc7e4b02220f; "
                    "shape 7312, 5582 true, 1730 false"
                ),
                "train_score_genes": "7312 master genes - 600 training targets = 6712",
            },
            "parser_defaults": defaults,
        },
        "fit_counts": fit_counts(grid),
        "flops": default_flop_accounting(dim, grid),
        "memory": array_memory_accounting(dim),
        "semantic_interactions": {
            "basis_depends_on": [
                "outer-role query",
                "held-teacher fit donor set",
                "fit rows",
                "correction axis",
                "basis_dim",
                "basis_seed",
            ],
            "basis_does_not_depend_on": [
                "rank",
                "beta",
                "temperature",
                "regularization",
                "real versus null labels",
                "null replicate",
            ],
            "normalized_arrays_and_training_pair_do_not_depend_on": [
                "rank",
                "beta",
                "temperature",
                "regularization",
                "real versus null labels",
                "null replicate",
            ],
            "optimizer_initialization_depends_on": [
                "outer role",
                "held teacher",
                "rank",
            ],
            "optimizer_trajectory_depends_on": [
                "rank",
                "beta",
                "temperature",
                "regularization",
                "labels/null replicate",
            ],
            "null_semantics_that_must_not_change": (
                "Each of 100 nulls must still train all 108 candidates in all "
                "four inner cells, select independently, refit its own winner, "
                "and participate as that selected null in the pretest/final null."
            ),
        },
        "recommendation": {
            "smallest_safe_change": (
                "Introduce an immutable PreparedFit cache keyed by protocol "
                "cell (stage, fold/outer role, held-teacher set, row split, "
                "basis seed). Reuse its exact basis, normalized projections, "
                "base score matrices, and raw response norms across candidates "
                "and all real/null repetitions. Keep the existing optimizer, "
                "candidate loop order, evaluate_adapter, choose_candidate, and "
                "all 100 complete null selections unchanged."
            ),
            "why_safe": (
                "The cached objects are mathematically and byte-for-byte "
                "candidate/null invariant in the pinned implementation. The "
                "representative benchmark checks exact basis/U/V/audit equality "
                "for ranks 4/8/16 and for a permuted-label null."
            ),
            "staged_benchmark": [
                "First patch only the inner four selection problem caches.",
                "On a guarded worker, run one real selection and nulls 0 and 1 both ways.",
                "Require byte-identical selected keys, cell deltas, angles, label hashes, U/V, and selection JSON before scaling to 100 nulls.",
                "Only afterward consider cached pretest/final problems.",
                "Defer latent-space evaluation and candidate batching until exhaustive real-data rank/tie equivalence is proved.",
            ],
        },
        "benchmark": benchmark,
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=HERE / "benchmark-results.json"
    )
    parser.add_argument("--benchmark-epochs", type=int, default=3)
    parser.add_argument("--benchmark-seed", type=int, default=8402026)
    args = parser.parse_args()
    if args.benchmark_epochs < 1 or args.benchmark_epochs > 10:
        raise SystemExit("--benchmark-epochs must be in [1,10]")
    if args.output.resolve().parent != HERE.resolve():
        raise SystemExit("--output must stay inside experiment 084")
    report = build_report(
        benchmark_epochs=args.benchmark_epochs,
        benchmark_seed=args.benchmark_seed,
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "peak_rss_gib": report["benchmark"]["peak_rss_gib"],
        "mandatory_grid_fits": report["fit_counts"]["mandatory_grid_fits"],
        "safe_cache_model_speedup": report["benchmark"][
            "local_single_thread_extrapolation"
        ]["safe_cache_model_speedup"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
