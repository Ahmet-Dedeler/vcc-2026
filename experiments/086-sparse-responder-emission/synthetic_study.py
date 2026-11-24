#!/usr/bin/env python3
"""Bounded full-shape and small DE diagnostics for sparse responders."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, mannwhitneyu


HERE = Path(__file__).resolve().parent


def load_emitter():
    spec = importlib.util.spec_from_file_location(
        "exp086_emitter_study", HERE / "emitter.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


EMITTER = load_emitter()


def exact_quasi_from_totals(
    totals: np.ndarray,
    *,
    cells: int,
    spread: float,
    seed: int,
) -> sp.csr_matrix:
    """Small reference implementation of experiment-077 exact-total quasi."""
    totals = np.asarray(totals, dtype=np.int64)
    rng = np.random.default_rng(seed)
    random_totals = np.floor(totals * float(spread)).astype(np.int64)
    even_totals = totals - random_totals
    block = np.zeros((cells, len(totals)), dtype=np.int64)
    nonzero = np.flatnonzero(even_totals)
    if nonzero.size:
        base, remainder = np.divmod(even_totals[nonzero], cells)
        offsets = np.concatenate(([0], np.cumsum(remainder)[:-1])) % cells
        position = (
            np.arange(cells, dtype=np.int64)[:, None] - offsets[None, :]
        ) % cells
        block[:, nonzero] = base[None, :] + (position < remainder[None, :])
    random_genes = np.flatnonzero(random_totals)
    if random_genes.size:
        draw = rng.poisson(
            np.broadcast_to(
                random_totals[random_genes] / cells,
                (cells, len(random_genes)),
            )
        )
        residual = random_totals[random_genes] - draw.sum(axis=0)
        for local in np.flatnonzero(residual):
            difference = int(residual[local])
            if difference > 0:
                destinations = rng.integers(0, cells, difference)
                np.add.at(draw[:, local], destinations, 1)
            else:
                for _ in range(-difference):
                    occupied = np.flatnonzero(draw[:, local] > 0)
                    draw[rng.choice(occupied), local] -= 1
        block[:, random_genes] += draw
    result = sp.csr_matrix(block)
    result.eliminate_zeros()
    return result


def hamilton_totals(profile: np.ndarray, grand_total: int) -> np.ndarray:
    weights = np.asarray(profile, dtype=np.float64)
    weights = weights / weights.sum()
    ideal = weights * int(grand_total)
    totals = np.floor(ideal).astype(np.int64)
    left = int(grand_total) - int(totals.sum())
    order = np.argsort(-(ideal - totals), kind="stable")
    totals[order[:left]] += 1
    return totals


def de_table(
    perturbation: sp.csr_matrix,
    control: sp.csr_matrix,
    *,
    gate: np.ndarray,
) -> dict[str, np.ndarray]:
    def cpm(matrix: sp.csr_matrix) -> np.ndarray:
        values = matrix.toarray().astype(np.float64)
        libraries = values.sum(axis=1)
        return values * (1_000_000.0 / libraries)[:, None]

    pert = cpm(perturbation)
    ctrl = cpm(control)
    columns = np.flatnonzero(gate)
    test = mannwhitneyu(
        pert[:, columns],
        ctrl[:, columns],
        alternative="two-sided",
        axis=0,
        method="asymptotic",
    )
    pvalues = np.nan_to_num(np.asarray(test.pvalue), nan=1.0)
    adjusted = false_discovery_control(pvalues, method="bh")
    lfc = np.log2(
        (pert[:, columns].mean(axis=0) + 1e-9) / (ctrl[:, columns].mean(axis=0) + 1e-9)
    )
    return {
        "columns": columns,
        "p": pvalues,
        "padj": adjusted,
        "lfc": lfc,
        "significant": adjusted < 0.05,
    }


def compare_de(
    predicted: dict[str, np.ndarray],
    real: dict[str, np.ndarray],
    *,
    active: np.ndarray,
) -> dict[str, float | int]:
    pred_sig = predicted["significant"]
    real_sig = real["significant"]
    intersection = pred_sig & real_sig
    union = pred_sig | real_sig
    return {
        "predicted_significant": int(pred_sig.sum()),
        "real_significant": int(real_sig.sum()),
        "jaccard": float(intersection.sum() / union.sum()) if union.any() else 1.0,
        "direction_accuracy_intersection": (
            float(
                np.mean(
                    np.sign(predicted["lfc"][intersection])
                    == np.sign(real["lfc"][intersection])
                )
            )
            if intersection.any()
            else 0.0
        ),
        "false_positive_calls_outside_planted_effect": int(
            (pred_sig & ~active[predicted["columns"]]).sum()
        ),
    }


def heterogeneous_controls(
    rng: np.random.Generator,
    *,
    cells: int,
    genes: int,
    median_library: int,
) -> sp.csr_matrix:
    """A bounded zero-inflated control fixture with realistic rank skew."""
    baseline = rng.lognormal(mean=-1.0, sigma=1.5, size=genes)
    baseline /= baseline.sum()
    libraries = rng.integers(
        int(0.70 * median_library),
        int(1.30 * median_library) + 1,
        size=cells,
    )
    rows = []
    for library in libraries:
        weights = baseline * rng.gamma(shape=0.35, scale=1.0 / 0.35, size=genes)
        weights /= weights.sum()
        rows.append(rng.multinomial(int(library), weights))
    return sp.csr_matrix(np.asarray(rows, dtype=np.int32))


def small_de_study(seed: int) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    genes = 512
    emitted_cells = 400
    control_cells = 2_000
    control = heterogeneous_controls(
        rng,
        cells=control_cells,
        genes=genes,
        median_library=4_000,
    )

    true_responder_fraction = 0.40
    real_template_indices = rng.choice(control_cells, emitted_cells, replace=False)
    templates = control[real_template_indices]
    baseline_columns = np.asarray(templates.sum(axis=0), dtype=np.int64).ravel()
    grand_total = int(baseline_columns.sum())
    control_profile = baseline_columns.astype(np.float64) / grand_total

    # A true response changes exactly 40 genes: mass is transferred from twenty
    # down genes to twenty up genes, leaving every other aggregate margin fixed.
    true_profile = control_profile.copy()
    transferable = min(0.02, 0.45 * float(true_profile[20:40].sum()))
    down_weights = true_profile[20:40] / true_profile[20:40].sum()
    up_weights = true_profile[:20] / true_profile[:20].sum()
    true_profile[20:40] -= transferable * down_weights
    true_profile[:20] += transferable * up_weights
    true_totals = hamilton_totals(true_profile, grand_total)
    real, _ = EMITTER.emit_sparse_responders(
        templates,
        column_totals=true_totals,
        responder_fraction=true_responder_fraction,
        key="synthetic-real",
        seed=seed + 50,
        control_reference_totals=baseline_columns,
        scored_gene_mask=np.ones(genes, dtype=bool),
    )

    active = np.zeros(genes, dtype=bool)
    active[:40] = True

    # Null-margin control: both layouts have the selected templates' exact
    # columns, but quasi flattens their ranks while sparse emission is an exact
    # no-op. This isolates the source of quasi's false-positive DE calls.
    null_quasi = exact_quasi_from_totals(
        baseline_columns,
        cells=emitted_cells,
        spread=0.25,
        seed=seed + 1,
    )
    null_sparse, _ = EMITTER.emit_sparse_responders(
        templates,
        column_totals=baseline_columns,
        responder_fraction=0.40,
        key="synthetic-null",
        seed=seed + 2,
        control_reference_totals=baseline_columns,
        scored_gene_mask=np.ones(genes, dtype=bool),
    )

    # Low-alpha model: correct direction but only 40% of the aggregate change.
    model_profile = control_profile + 0.40 * (true_profile - control_profile)
    totals = hamilton_totals(model_profile, grand_total)

    quasi = exact_quasi_from_totals(
        totals, cells=emitted_cells, spread=0.25, seed=seed + 3
    )
    sparse, diagnostics = EMITTER.emit_sparse_responders(
        templates,
        column_totals=totals,
        responder_fraction=0.40,
        key="synthetic-de",
        seed=seed + 4,
        control_reference_totals=baseline_columns,
        scored_gene_mask=np.ones(genes, dtype=bool),
    )
    gate = EMITTER.arithmetic_mean_cpm(control) > 5.0
    null_quasi_de = de_table(null_quasi, control, gate=gate)
    null_sparse_de = de_table(null_sparse, control, gate=gate)
    real_de = de_table(real, control, gate=gate)
    quasi_de = de_table(quasi, control, gate=gate)
    sparse_de = de_table(sparse, control, gate=gate)
    quasi_bulk = EMITTER.scorer_pseudobulk(quasi)
    sparse_bulk = EMITTER.scorer_pseudobulk(sparse)
    return {
        "scenario": (
            "synthetic truth has 40% responders; model uses 40% of its aggregate "
            "effect; zero-inflated controls expose distribution-only false calls"
        ),
        "planted_effect_genes": int(active.sum()),
        "pseudobulk_array_equal": bool(np.array_equal(quasi_bulk, sparse_bulk)),
        "max_abs_pseudobulk_difference": float(
            np.max(np.abs(quasi_bulk - sparse_bulk))
        ),
        "null_margin_control": {
            "pseudobulk_array_equal": bool(
                np.array_equal(
                    EMITTER.scorer_pseudobulk(null_quasi),
                    EMITTER.scorer_pseudobulk(null_sparse),
                )
            ),
            "quasi_significant": int(null_quasi_de["significant"].sum()),
            "sparse_significant": int(null_sparse_de["significant"].sum()),
        },
        "low_alpha_quasi": {
            **compare_de(quasi_de, real_de, active=active),
            "jackknife_correction": EMITTER.jackknife_correction(quasi),
        },
        "low_alpha_sparse_responder": {
            **compare_de(sparse_de, real_de, active=active),
            "jackknife_correction": EMITTER.jackknife_correction(sparse),
            "diagnostics": diagnostics.to_dict(),
        },
    }


def full_shape_memory_study(seed: int) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    cells, genes = 400, 18_533
    # Sparse but full-dimensional representative control templates.
    dense = rng.poisson(1.4, size=(cells, genes)).astype(np.int16)
    dense[rng.random((cells, genes)) < 0.72] = 0
    for row in np.flatnonzero(dense.sum(axis=1) == 0):
        dense[row, 0] = 1
    template = sp.csr_matrix(dense, dtype=np.int32)
    control_totals = np.asarray(template.sum(axis=0), dtype=np.int64).ravel()
    multipliers = np.ones(genes)
    multipliers[:500] = np.linspace(0.65, 0.95, 500)
    multipliers[500:1_000] = np.linspace(1.05, 1.45, 500)
    weighted = control_totals.astype(np.float64) * multipliers
    ideal = weighted * (float(control_totals.sum()) / weighted.sum())
    target = np.floor(ideal).astype(np.int64)
    left = int(control_totals.sum() - target.sum())
    order = np.argsort(-(ideal - target), kind="stable")
    target[order[:left]] += 1
    scored = control_totals > np.quantile(control_totals, 0.20)
    started = time.perf_counter()
    emitted, diagnostics = EMITTER.emit_sparse_responders(
        template,
        column_totals=target,
        responder_fraction=0.50,
        key="full-shape",
        seed=seed + 1,
        control_reference_totals=control_totals,
        scored_gene_mask=scored,
    )
    elapsed = time.perf_counter() - started
    return {
        "shape": list(emitted.shape),
        "elapsed_seconds": elapsed,
        "peak_rss_gib": EMITTER.process_peak_rss_gib(),
        "column_totals_exact": bool(
            np.array_equal(
                np.asarray(emitted.sum(axis=0), dtype=np.int64).ravel(),
                target,
            )
        ),
        "diagnostics": diagnostics.to_dict(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "synthetic-results.json")
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()
    if args.output.resolve().parent != HERE.resolve():
        raise SystemExit("--output must stay inside experiment 086")
    result = {
        "study": "086-sparse-responder-emission",
        "small_de": small_de_study(args.seed),
        "full_shape_memory": full_shape_memory_study(args.seed + 100),
    }
    result["peak_rss_gib"] = EMITTER.process_peak_rss_gib()
    if result["peak_rss_gib"] >= 1.0:
        raise MemoryError(
            f"synthetic study exceeded 1 GiB: {result['peak_rss_gib']:.3f}"
        )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
