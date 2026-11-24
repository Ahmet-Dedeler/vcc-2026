"""Bounded real-data diagnostic for the fixed-margin template emitter.

Reads one 400-cell quasi-emitted perturbation block from experiment 040 and
compares it with a fixed-margin re-emission based on real Orion-HCT116 controls.
Only backed slices are read; this script does not invoke cell-eval2.
"""

from __future__ import annotations

import argparse
import json
import math
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import (
    false_discovery_control,
    ks_2samp,
    mannwhitneyu,
    spearmanr,
)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from template_emitter import (  # noqa: E402
    arithmetic_mean_cpm,
    emit_template_fixed_margins,
    integer_column_totals,
    scorer_pseudobulk,
    template_scaled_cpm_target,
)


CONTROL = "non-targeting"
PERT_COL = "target_gene"
GIB = 1024**3
QUANTILES = np.array([0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])


def peak_rss_bytes() -> int:
    """ru_maxrss is bytes on macOS and KiB on Linux."""
    rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss if sys.platform == "darwin" else rss * 1024


def memory_checkpoint(stage: str, limit_bytes: int) -> None:
    rss = peak_rss_bytes()
    print(f"[rss] {stage}: {rss / GIB:.3f} GiB", flush=True)
    if rss >= limit_bytes:
        raise MemoryError(
            f"peak RSS {rss / GIB:.3f} GiB exceeded limit "
            f"{limit_bytes / GIB:.3f} GiB at {stage}"
        )


def contiguous_rows(indices: np.ndarray) -> slice | np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size and np.all(np.diff(indices) == 1):
        return slice(int(indices[0]), int(indices[-1]) + 1)
    return np.sort(indices)


def read_rows(backed: ad.AnnData, indices: np.ndarray) -> sp.csr_matrix:
    """Read only requested rows from a backed CSR dataset."""
    x = backed.X[contiguous_rows(indices)]
    return x.tocsr() if not sp.isspmatrix_csr(x) else x


def align_real_rows(
    real_rows: sp.csr_matrix,
    real_cols_in_pred_order: np.ndarray,
    pred_positions: np.ndarray,
    n_pred_genes: int,
) -> sp.csr_matrix:
    """Map a real-axis CSR slice onto the experiment-040 prediction axis."""
    common = real_rows[:, real_cols_in_pred_order].tocoo()
    return sp.csr_matrix(
        (
            common.data,
            (common.row, pred_positions[common.col]),
        ),
        shape=(real_rows.shape[0], n_pred_genes),
    )


def _row_scaled_column_sum(
    x: sp.csr_matrix,
    row_denominators: np.ndarray,
    scale: float = 1_000_000.0,
) -> np.ndarray:
    x = x.tocsr().astype(np.float64)
    denom = np.asarray(row_denominators, dtype=np.float64)
    if denom.shape != (x.shape[0],) or np.any(denom <= 0):
        raise ValueError("row denominators must be positive and match matrix rows")
    row_of_nnz = np.repeat(np.arange(x.shape[0]), np.diff(x.indptr))
    x.data *= scale / denom[row_of_nnz]
    return np.asarray(x.sum(axis=0), dtype=np.float64).ravel()


def stream_control_summaries(
    real: ad.AnnData,
    control_indices: np.ndarray,
    real_cols_in_pred_order: np.ndarray,
    *,
    chunk_rows: int,
    rss_limit: int,
) -> dict[str, np.ndarray | float]:
    """One bounded pass over all real controls."""
    n = control_indices.size
    genes = real_cols_in_pred_order.size
    libs = np.empty(n, dtype=np.int64)
    raw_common = np.zeros(genes, dtype=np.float64)
    cpm_sum_common = np.zeros(genes, dtype=np.float64)
    detected_common = np.zeros(genes, dtype=np.int64)
    total_all = 0.0

    for start in range(0, n, chunk_rows):
        stop = min(n, start + chunk_rows)
        idx = control_indices[start:stop]
        full = read_rows(real, idx)
        chunk_libs = np.asarray(full.sum(axis=1), dtype=np.int64).ravel()
        common = full[:, real_cols_in_pred_order].tocsr()
        libs[start:stop] = chunk_libs
        total_all += float(chunk_libs.sum())
        raw_common += np.asarray(common.sum(axis=0), dtype=np.float64).ravel()
        cpm_sum_common += _row_scaled_column_sum(common, chunk_libs)
        detected_common += np.diff(common.tocsc().indptr)
        if start % (8 * chunk_rows) == 0:
            memory_checkpoint(f"streamed controls {stop}/{n}", rss_limit)

    return {
        "library_sizes": libs,
        "raw_common": raw_common,
        "mean_cpm_common": cpm_sum_common / n,
        "detection_rate_common": detected_common / n,
        "total_counts_all_genes": total_all,
    }


def stratified_positions(values: np.ndarray, n: int) -> np.ndarray:
    """Select approximately equal-mass quantiles without replacement."""
    if n > values.size:
        raise ValueError(f"cannot select {n} rows from {values.size}")
    order = np.argsort(values, kind="stable")
    pos = np.floor((np.arange(n) + 0.5) * values.size / n).astype(np.int64)
    picked = order[pos]
    if np.unique(picked).size != n:
        raise RuntimeError("quantile selection unexpectedly produced duplicate rows")
    return picked


def quantile_summary(values: np.ndarray) -> dict[str, float]:
    q = np.quantile(np.asarray(values, dtype=np.float64), QUANTILES)
    return {f"q{int(100 * p):02d}": float(v) for p, v in zip(QUANTILES, q, strict=True)}


def distribution_summary(values: np.ndarray) -> dict[str, object]:
    x = np.asarray(values, dtype=np.float64)
    mean = float(x.mean())
    return {
        "mean": mean,
        "median": float(np.median(x)),
        "cv": float(x.std(ddof=1) / mean) if mean else 0.0,
        "quantiles": quantile_summary(x),
        "relative_to_median_quantiles": quantile_summary(x / np.median(x)),
    }


def safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    value = float(spearmanr(a, b).statistic)
    return value if math.isfinite(value) else 0.0


def support_jaccard(a: np.ndarray, b: np.ndarray, threshold: float) -> float:
    sa = np.asarray(a) > threshold
    sb = np.asarray(b) > threshold
    union = int(np.count_nonzero(sa | sb))
    return 1.0 if union == 0 else float(np.count_nonzero(sa & sb) / union)


def detection_summary(
    matrix: sp.csr_matrix,
    real_target_rate: np.ndarray,
    real_control_rate: np.ndarray,
) -> dict[str, object]:
    rate = np.diff(matrix.tocsc().indptr) / matrix.shape[0]
    return {
        "mean_detected_genes_per_cell": float(np.diff(matrix.indptr).mean()),
        "detected_genes_per_cell_quantiles": quantile_summary(np.diff(matrix.indptr)),
        "gene_detection_rate_spearman_vs_real_target": safe_spearman(
            rate, real_target_rate
        ),
        "gene_detection_rate_mae_vs_real_target": float(
            np.mean(np.abs(rate - real_target_rate))
        ),
        "gene_detection_rate_spearman_vs_control": safe_spearman(
            rate, real_control_rate
        ),
        "gene_detection_rate_mae_vs_control": float(
            np.mean(np.abs(rate - real_control_rate))
        ),
        "detected_gene_jaccard_vs_real_target_any": support_jaccard(
            rate, real_target_rate, 0.0
        ),
        "detected_gene_jaccard_vs_real_target_5pct": support_jaccard(
            rate, real_target_rate, 0.05
        ),
    }


def cellwise_support_jaccard(a: sp.csr_matrix, b: sp.csr_matrix) -> dict[str, float]:
    aa = a.astype(bool)
    bb = b.astype(bool)
    inter = np.asarray(aa.multiply(bb).sum(axis=1)).ravel()
    union = np.diff(aa.indptr) + np.diff(bb.indptr) - inter
    score = np.divide(
        inter,
        union,
        out=np.ones_like(inter, dtype=np.float64),
        where=union > 0,
    )
    return {
        "mean": float(score.mean()),
        "median": float(np.median(score)),
        "q05": float(np.quantile(score, 0.05)),
        "q95": float(np.quantile(score, 0.95)),
    }


def mean_cpm(
    matrix: sp.csr_matrix,
    full_library_sizes: np.ndarray | None = None,
) -> np.ndarray:
    if full_library_sizes is None:
        return arithmetic_mean_cpm(matrix)
    return _row_scaled_column_sum(matrix, full_library_sizes) / matrix.shape[0]


def cpm_comparison(
    predicted: np.ndarray,
    *,
    nominal: np.ndarray,
    template_target: np.ndarray,
    real_target: np.ndarray,
    real_control: np.ndarray,
    gate: np.ndarray,
) -> dict[str, float]:
    eps = 1e-9

    def rel_l1(a, b):
        denom = float(np.abs(b).sum())
        return float(np.abs(a - b).sum() / denom) if denom else 0.0

    real_lfc = np.log2((real_target + eps) / (real_control + eps))
    pred_lfc = np.log2((predicted + eps) / (real_control + eps))
    denom = float(np.abs(real_lfc[gate]).sum())
    return {
        "relative_l1_vs_nominal_column_cpm": rel_l1(predicted, nominal),
        "relative_l1_vs_template_scaled_cpm": rel_l1(predicted, template_target),
        "relative_l1_vs_real_target_cpm": rel_l1(predicted, real_target),
        "mean_cpm_spearman_vs_real_target": safe_spearman(predicted, real_target),
        "lfc_spearman_vs_real": safe_spearman(pred_lfc[gate], real_lfc[gate]),
        "lfc_nmae_proxy_vs_real": (
            float(np.abs(pred_lfc[gate] - real_lfc[gate]).sum() / denom)
            if denom
            else 0.0
        ),
    }


def mwu_against_control(
    arms: dict[str, tuple[sp.csr_matrix, np.ndarray]],
    control: sp.csr_matrix,
    control_libs: np.ndarray,
    gate_common: np.ndarray,
    *,
    gene_chunk: int,
    rss_limit: int,
) -> dict[str, dict[str, np.ndarray]]:
    """Chunked two-sided Mann-Whitney U against one sampled real-control pool."""
    gate_cols = np.flatnonzero(gate_common)
    out = {
        name: {
            "p": np.ones(gate_cols.size, dtype=np.float64),
            "rbc": np.zeros(gate_cols.size, dtype=np.float64),
        }
        for name in arms
    }
    ctrl_scale = 1_000_000.0 / np.asarray(control_libs, dtype=np.float64)
    arm_scales = {
        name: 1_000_000.0 / np.asarray(libs, dtype=np.float64)
        for name, (_, libs) in arms.items()
    }

    for start in range(0, gate_cols.size, gene_chunk):
        stop = min(gate_cols.size, start + gene_chunk)
        cols = gate_cols[start:stop]
        ctrl = control[:, cols].toarray().astype(np.float64, copy=False)
        ctrl *= ctrl_scale[:, None]
        for name, (matrix, _) in arms.items():
            values = matrix[:, cols].toarray().astype(np.float64, copy=False)
            values *= arm_scales[name][:, None]
            test = mannwhitneyu(
                values,
                ctrl,
                alternative="two-sided",
                axis=0,
                method="asymptotic",
            )
            out[name]["p"][start:stop] = np.asarray(test.pvalue, dtype=np.float64)
            u = np.asarray(test.statistic, dtype=np.float64)
            out[name]["rbc"][start:stop] = (
                2.0 * u / (matrix.shape[0] * control.shape[0]) - 1.0
            )
        if start % (8 * gene_chunk) == 0:
            memory_checkpoint(
                f"MWU genes {stop}/{gate_cols.size}", rss_limit
            )
    return out


def de_rank_summary(
    stats: dict[str, dict[str, np.ndarray]],
    mean_cpms: dict[str, np.ndarray],
    control_mean_cpm: np.ndarray,
    gate_common: np.ndarray,
) -> dict[str, dict[str, object]]:
    gate_cols = np.flatnonzero(gate_common)
    real = stats["real"]
    real_sign = np.sign(
        mean_cpms["real"][gate_cols] - control_mean_cpm[gate_cols]
    )
    real_score = real_sign * -np.log10(np.clip(real["p"], 1e-300, 1.0))
    real_order = np.lexsort((-np.abs(real["rbc"]), real["p"]))
    real_adj = false_discovery_control(real["p"], method="bh")
    real_sig = real_adj < 0.05

    result: dict[str, dict[str, object]] = {}
    for name in ("quasi", "template"):
        pred = stats[name]
        pred_sign = np.sign(
            mean_cpms[name][gate_cols] - control_mean_cpm[gate_cols]
        )
        pred_score = pred_sign * -np.log10(np.clip(pred["p"], 1e-300, 1.0))
        pred_order = np.lexsort((-np.abs(pred["rbc"]), pred["p"]))
        pred_adj = false_discovery_control(pred["p"], method="bh")
        pred_sig = pred_adj < 0.05
        inter = pred_sig & real_sig
        union = pred_sig | real_sig

        top: dict[str, dict[str, float | int]] = {}
        for k in (100, 500, 1000):
            kk = min(k, gate_cols.size)
            rt = real_order[:kk]
            pt = pred_order[:kk]
            overlap = np.intersect1d(rt, pt).size
            top[str(k)] = {
                "k_used": int(kk),
                "set_overlap_fraction": float(overlap / kk),
                "direction_accuracy_on_real_top": float(
                    np.mean(pred_sign[rt] == real_sign[rt])
                ),
                "direction_accuracy_on_pred_top": float(
                    np.mean(pred_sign[pt] == real_sign[pt])
                ),
            }

        result[name] = {
            "signed_neglog10p_spearman_vs_real": safe_spearman(
                pred_score, real_score
            ),
            "rank_biserial_spearman_vs_real": safe_spearman(
                pred["rbc"], real["rbc"]
            ),
            "bh05_n_sig": int(pred_sig.sum()),
            "real_bh05_n_sig": int(real_sig.sum()),
            "bh05_jaccard_vs_real": (
                float(inter.sum() / union.sum()) if union.any() else 1.0
            ),
            "bh05_recall_vs_real": (
                float(inter.sum() / real_sig.sum()) if real_sig.any() else 1.0
            ),
            "bh05_direction_accuracy_in_intersection": (
                float(np.mean(pred_sign[inter] == real_sign[inter]))
                if inter.any()
                else 0.0
            ),
            "top_k": top,
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prediction",
        default=(
            "experiments/040-union-gamma/out/runs/"
            "g0.75_q1.0_a1.08_spread0.25/pred.h5ad"
        ),
    )
    parser.add_argument("--real", default="data/harness/orion-hct116.h5ad")
    parser.add_argument(
        "--output",
        default="experiments/043-template-emitter/real_data_results.json",
    )
    parser.add_argument("--template-cells", type=int, default=400)
    parser.add_argument("--de-control-cells", type=int, default=2000)
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--gene-chunk", type=int, default=128)
    parser.add_argument("--seed", type=int, default=4300)
    parser.add_argument("--max-rss-gib", type=float, default=3.8)
    args = parser.parse_args()

    started = time.time()
    rss_limit = int(args.max_rss_gib * GIB)
    prediction_path = Path(args.prediction)
    real_path = Path(args.real)
    output_path = Path(args.output)
    memory_checkpoint("start", rss_limit)

    pred = ad.read_h5ad(prediction_path, backed="r")
    real = ad.read_h5ad(real_path, backed="r")
    try:
        pred_genes = np.asarray(pred.var_names, dtype=str)
        real_genes = np.asarray(real.var_names, dtype=str)
        real_lookup = {gene: i for i, gene in enumerate(real_genes)}
        real_map = np.array([real_lookup.get(gene, -1) for gene in pred_genes])
        pred_positions = np.flatnonzero(real_map >= 0)
        real_cols = real_map[pred_positions]

        pred_labels = pred.obs[PERT_COL].astype(str).to_numpy()
        real_labels = real.obs[PERT_COL].astype(str).to_numpy()
        pred_targets = set(pred_labels) - {CONTROL}
        real_targets = set(real_labels) - {CONTROL}
        overlap = sorted(pred_targets & real_targets)
        if not overlap:
            raise RuntimeError("prediction and HCT116 reference share no perturbations")

        real_control_idx = np.flatnonzero(real_labels == CONTROL)
        control = stream_control_summaries(
            real,
            real_control_idx,
            real_cols,
            chunk_rows=args.chunk_rows,
            rss_limit=rss_limit,
        )
        ctrl_bulk = np.log1p(
            50_000.0
            * np.asarray(control["raw_common"])
            / float(control["total_counts_all_genes"])
        )
        memory_checkpoint("control summaries", rss_limit)

        candidate_norms: dict[str, float] = {}
        for target in overlap:
            idx = np.flatnonzero(real_labels == target)
            full = read_rows(real, idx)
            total = float(full.sum())
            common_cols = np.asarray(
                full[:, real_cols].sum(axis=0), dtype=np.float64
            ).ravel()
            bulk = np.log1p(50_000.0 * common_cols / total)
            keep = pred_genes[pred_positions] != target
            candidate_norms[target] = float(
                np.linalg.norm((bulk - ctrl_bulk)[keep]) / np.sqrt(keep.sum())
            )
        target = max(candidate_norms, key=candidate_norms.get)
        print(
            f"selected {target} from overlap {overlap}; "
            f"real bulk RMS effects={candidate_norms}",
            flush=True,
        )

        pred_target_idx = np.flatnonzero(pred_labels == target)
        real_target_idx = np.flatnonzero(real_labels == target)
        if pred_target_idx.size != 400:
            raise RuntimeError(
                f"expected 400 predicted cells for {target}, got {pred_target_idx.size}"
            )
        quasi = read_rows(pred, pred_target_idx)
        if not np.array_equal(quasi.data, np.rint(quasi.data)):
            raise RuntimeError("experiment-040 block is not integer-valued")
        quasi.data = np.rint(quasi.data).astype(np.int64)
        quasi.eliminate_zeros()
        quasi_cols = np.asarray(quasi.sum(axis=0), dtype=np.int64).ravel()
        quasi_total = int(quasi_cols.sum())
        memory_checkpoint("prediction block", rss_limit)

        control_libs = np.asarray(control["library_sizes"], dtype=np.int64)
        template_local = stratified_positions(control_libs, args.template_cells)
        de_control_local = stratified_positions(control_libs, args.de_control_cells)
        template_global = np.sort(real_control_idx[template_local])
        de_control_global = np.sort(real_control_idx[de_control_local])

        template_full = read_rows(real, template_global)
        template_full_libs = np.asarray(
            template_full.sum(axis=1), dtype=np.int64
        ).ravel()
        template_aligned = align_real_rows(
            template_full, real_cols, pred_positions, pred.n_vars
        )
        row_margins = integer_column_totals(
            template_full_libs, quasi_total, seed=args.seed + 1
        )
        template, emitter_diag = emit_template_fixed_margins(
            template_aligned,
            column_totals=quasi_cols,
            row_margins=row_margins,
            template_library_sizes=template_full_libs,
            seed=args.seed,
        )
        memory_checkpoint("template emission", rss_limit)

        # Load only the selected real target and the 2,000-cell control sample.
        real_target_full = read_rows(real, real_target_idx)
        real_target_libs = np.asarray(
            real_target_full.sum(axis=1), dtype=np.int64
        ).ravel()
        real_target_common = real_target_full[:, real_cols].tocsr()

        de_control_full = read_rows(real, de_control_global)
        de_control_libs = np.asarray(
            de_control_full.sum(axis=1), dtype=np.int64
        ).ravel()
        de_control_common = de_control_full[:, real_cols].tocsr()
        memory_checkpoint("real target and control sample", rss_limit)

        quasi_common = quasi[:, pred_positions].tocsr()
        template_common = template[:, pred_positions].tocsr()
        template_source_common = template_aligned[:, pred_positions].tocsr()
        real_target_rate = (
            np.diff(real_target_common.tocsc().indptr) / real_target_common.shape[0]
        )
        real_control_rate = np.asarray(control["detection_rate_common"])

        quasi_libs = np.asarray(quasi.sum(axis=1), dtype=np.int64).ravel()
        template_libs = np.asarray(
            template.sum(axis=1), dtype=np.int64
        ).ravel()
        library = {
            "hct116_control_all": distribution_summary(control_libs),
            "hct116_control_template_original": distribution_summary(
                template_full_libs
            ),
            "quasi": distribution_summary(quasi_libs),
            "template_fixed_margin": distribution_summary(template_libs),
            "shape_ks_vs_hct116_control": {
                "quasi": float(
                    ks_2samp(
                        quasi_libs / np.median(quasi_libs),
                        control_libs / np.median(control_libs),
                    ).statistic
                ),
                "template_fixed_margin": float(
                    ks_2samp(
                        template_libs / np.median(template_libs),
                        control_libs / np.median(control_libs),
                    ).statistic
                ),
            },
        }

        detection = {
            "real_target": {
                "mean_detected_genes_per_cell": float(
                    np.diff(real_target_common.indptr).mean()
                ),
                "detected_genes_per_cell_quantiles": quantile_summary(
                    np.diff(real_target_common.indptr)
                ),
            },
            "real_control_all": {
                "mean_detected_genes_per_cell": float(
                    np.sum(real_control_rate)
                ),
            },
            "quasi": detection_summary(
                quasi_common, real_target_rate, real_control_rate
            ),
            "template_fixed_margin": detection_summary(
                template_common, real_target_rate, real_control_rate
            ),
            "template_cellwise_support_jaccard_to_source": (
                cellwise_support_jaccard(
                    template_common, template_source_common
                )
            ),
        }

        quasi_mean = mean_cpm(quasi_common)
        template_mean = mean_cpm(template_common)
        real_target_mean = mean_cpm(real_target_common, real_target_libs)
        real_control_mean = np.asarray(control["mean_cpm_common"])
        nominal = (
            1_000_000.0
            * quasi_cols[pred_positions].astype(np.float64)
            / quasi_total
        )
        template_target_full = template_scaled_cpm_target(
            template_aligned,
            quasi_cols,
            row_margins,
            template_library_sizes=template_full_libs,
        )
        template_target = template_target_full[pred_positions]
        gate = real_control_mean > 5.0
        target_common_position = np.flatnonzero(
            pred_genes[pred_positions] == target
        )
        if target_common_position.size:
            gate[target_common_position[0]] = False

        cpm = {
            "quasi": cpm_comparison(
                quasi_mean,
                nominal=nominal,
                template_target=template_target,
                real_target=real_target_mean,
                real_control=real_control_mean,
                gate=gate,
            ),
            "template_fixed_margin": cpm_comparison(
                template_mean,
                nominal=nominal,
                template_target=template_target,
                real_target=real_target_mean,
                real_control=real_control_mean,
                gate=gate,
            ),
            "tested_genes_control_cpm_gt_5_excluding_target": int(gate.sum()),
        }

        stats = mwu_against_control(
            {
                "real": (real_target_common, real_target_libs),
                "quasi": (quasi_common, quasi_libs),
                "template": (template_common, template_libs),
            },
            de_control_common,
            de_control_libs,
            gate,
            gene_chunk=args.gene_chunk,
            rss_limit=rss_limit,
        )
        de = de_rank_summary(
            stats,
            {
                "real": real_target_mean,
                "quasi": quasi_mean,
                "template": template_mean,
            },
            real_control_mean,
            gate,
        )
        memory_checkpoint("DE rank diagnostics", rss_limit)

        quasi_bulk = scorer_pseudobulk(quasi)
        template_bulk = scorer_pseudobulk(template)
        got_template_cols = np.asarray(
            template.sum(axis=0), dtype=np.int64
        ).ravel()
        pseudobulk = {
            "column_totals_exact": bool(
                np.array_equal(quasi_cols, got_template_cols)
            ),
            "max_abs_column_total_difference": int(
                np.max(np.abs(quasi_cols - got_template_cols), initial=0)
            ),
            "max_abs_bulk_lognorm_difference": float(
                np.max(np.abs(quasi_bulk - template_bulk), initial=0.0)
            ),
            "array_equal_bulk_lognorm": bool(
                np.array_equal(quasi_bulk, template_bulk)
            ),
        }

        result = {
            "date": "2026-08-22",
            "prediction": str(prediction_path),
            "real": str(real_path),
            "prediction_shape": list(pred.shape),
            "real_shape": list(real.shape),
            "prediction_target_count": len(pred_targets),
            "real_target_count": len(real_targets),
            "overlap_targets": overlap,
            "common_gene_count": int(pred_positions.size),
            "prediction_only_genes": pred_genes[real_map < 0].tolist(),
            "selected_target": target,
            "selection_rule": "largest HCT116 real bulk-log effect RMS among overlap targets",
            "selection_candidate_effect_rms": candidate_norms,
            "cells": {
                "quasi": int(quasi.shape[0]),
                "template": int(template.shape[0]),
                "real_target": int(real_target_common.shape[0]),
                "real_control_all": int(real_control_idx.size),
                "real_control_mwu_sample": int(de_control_common.shape[0]),
            },
            "pseudobulk": pseudobulk,
            "library_sizes": library,
            "detection_and_zeros": detection,
            "arithmetic_mean_cpm": cpm,
            "simple_mwu_rank_diagnostic": de,
            "template_emitter_diagnostics": emitter_diag.to_dict(),
            "peak_rss_gib": peak_rss_bytes() / GIB,
            "elapsed_seconds": time.time() - started,
            "notes": [
                "No cell-eval2 command or API was invoked.",
                "MWU diagnostics use 2,000 library-quantile-stratified real controls, "
                "not the official full-control scorer.",
                "The experiment-040 artifact has 85 targets and 18,080 genes, while "
                "Orion-HCT116 has 132 targets and 18,533 genes; diagnostics use their "
                "18,077 common genes.",
                "Fixed column totals force the template row depths to the quasi block's "
                "grand total. Their HCT116 library-size distribution shape is preserved "
                "up to one global scale and integer rounding.",
            ],
        }
        result["peak_rss_gib"] = peak_rss_bytes() / GIB
        output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    finally:
        pred.file.close()
        real.file.close()

    memory_checkpoint("done", rss_limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
