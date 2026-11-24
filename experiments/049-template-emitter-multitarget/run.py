#!/usr/bin/env python3
"""Bounded multi-target diagnostic for experiment 043's template emitter.

This compares the existing experiment-040 quasi-emitted cell blocks with
fixed-margin re-emissions for MED15, MED25, and PLAGL2.  Both arms have exactly
the same per-gene column totals, so their scorer pseudobulk is identical.

The experiment is intentionally not a model-validation result: the prediction
artifact has the H1-test target/layout, while the real perturbation cells and
control templates come from Orion-HCT116.  Only backed sparse slices are read,
and cell-eval2 is never imported or invoked.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.stats import ks_2samp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPPORT_DIR = ROOT / "experiments" / "043-template-emitter"
sys.path.insert(0, str(SUPPORT_DIR))

from real_data_diagnostic import (  # noqa: E402
    CONTROL,
    GIB,
    PERT_COL,
    align_real_rows,
    cpm_comparison,
    de_rank_summary,
    detection_summary,
    mean_cpm,
    memory_checkpoint,
    mwu_against_control,
    peak_rss_bytes,
    read_rows,
    stratified_positions,
    stream_control_summaries,
)
from template_emitter import (  # noqa: E402
    emit_template_fixed_margins,
    integer_column_totals,
    scorer_pseudobulk,
    template_scaled_cpm_target,
)


DEFAULT_TARGETS = ("MED15", "MED25", "PLAGL2")
DEFAULT_TARGET_SEED_OFFSETS = {"MED25": 0, "MED15": 100, "PLAGL2": 200}
ARM_NAMES = ("quasi", "template")
PRIMARY_METRICS = {
    "library_shape_ks_vs_hct116_control": "lower",
    "detection_rate_mae_vs_real_target": "lower",
    "de_set_jaccard_bh05": "higher",
    "direction_accuracy_bh05_intersection": "higher",
    "top500_set_overlap_fraction": "higher",
    "lfc_nmae_proxy_vs_real": "lower",
}


def arm_summary(
    *,
    library_ks: float,
    detection: dict[str, object],
    cpm: dict[str, float],
    de: dict[str, object],
) -> dict[str, float | int]:
    """Extract the requested comparable metrics for one emission arm."""
    top500 = de["top_k"]["500"]
    return {
        "library_shape_ks_vs_hct116_control": float(library_ks),
        "detection_rate_mae_vs_real_target": float(
            detection["gene_detection_rate_mae_vs_real_target"]
        ),
        "de_set_jaccard_bh05": float(de["bh05_jaccard_vs_real"]),
        "direction_accuracy_bh05_intersection": float(
            de["bh05_direction_accuracy_in_intersection"]
        ),
        "top500_set_overlap_fraction": float(top500["set_overlap_fraction"]),
        "lfc_nmae_proxy_vs_real": float(cpm["lfc_nmae_proxy_vs_real"]),
        "predicted_bh05_n_sig": int(de["bh05_n_sig"]),
        "real_bh05_n_sig": int(de["real_bh05_n_sig"]),
    }


def macro_summary(
    per_target: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Macro-average requested metrics and count target-level wins."""
    arms: dict[str, dict[str, float]] = {}
    for arm in ARM_NAMES:
        arms[arm] = {
            metric: float(
                np.mean(
                    [
                        float(per_target[target]["arms"][arm][metric])
                        for target in per_target
                    ]
                )
            )
            for metric in PRIMARY_METRICS
        }

    comparison: dict[str, dict[str, float | int | str]] = {}
    for metric, direction in PRIMARY_METRICS.items():
        raw_delta = arms["template"][metric] - arms["quasi"][metric]
        favorable_delta = raw_delta if direction == "higher" else -raw_delta
        wins = 0
        ties = 0
        for target in per_target:
            q = float(per_target[target]["arms"]["quasi"][metric])
            t = float(per_target[target]["arms"]["template"][metric])
            if np.isclose(q, t, rtol=0.0, atol=1e-12):
                ties += 1
            elif (direction == "higher" and t > q) or (
                direction == "lower" and t < q
            ):
                wins += 1
        comparison[metric] = {
            "preferred": direction,
            "template_minus_quasi": float(raw_delta),
            "favorable_delta_for_template": float(favorable_delta),
            "template_target_wins": wins,
            "ties": ties,
            "targets": len(per_target),
        }

    return {"arms": arms, "comparison": comparison}


def target_seed(base_seed: int, target: str, position: int) -> int:
    """Keep MED25 bitwise-comparable with experiment 043's seed."""
    return base_seed + DEFAULT_TARGET_SEED_OFFSETS.get(target, 1000 + 100 * position)


def evaluate_target(
    *,
    target: str,
    pred: ad.AnnData,
    real: ad.AnnData,
    pred_labels: np.ndarray,
    real_labels: np.ndarray,
    pred_genes: np.ndarray,
    pred_positions: np.ndarray,
    real_cols: np.ndarray,
    control: dict[str, np.ndarray | float],
    template_aligned,
    template_full_libs: np.ndarray,
    de_control_common,
    de_control_libs: np.ndarray,
    seed: int,
    gene_chunk: int,
    rss_limit: int,
) -> dict[str, object]:
    """Evaluate one target, releasing all target-sized arrays on return."""
    pred_target_idx = np.flatnonzero(pred_labels == target)
    real_target_idx = np.flatnonzero(real_labels == target)
    if pred_target_idx.size != 400:
        raise RuntimeError(
            f"expected 400 predicted cells for {target}, got {pred_target_idx.size}"
        )
    if real_target_idx.size == 0:
        raise RuntimeError(f"real reference has no rows for {target}")

    quasi = read_rows(pred, pred_target_idx)
    if not np.array_equal(quasi.data, np.rint(quasi.data)):
        raise RuntimeError(f"{target}: prediction block is not integer-valued")
    quasi.data = np.rint(quasi.data).astype(np.int64)
    quasi.eliminate_zeros()
    quasi_cols = np.asarray(quasi.sum(axis=0), dtype=np.int64).ravel()
    quasi_total = int(quasi_cols.sum())

    row_margins = integer_column_totals(
        template_full_libs, quasi_total, seed=seed + 1
    )
    template, emitter_diag = emit_template_fixed_margins(
        template_aligned,
        column_totals=quasi_cols,
        row_margins=row_margins,
        template_library_sizes=template_full_libs,
        seed=seed,
    )
    memory_checkpoint(f"{target}: template emission", rss_limit)

    real_target_full = read_rows(real, real_target_idx)
    real_target_libs = np.asarray(
        real_target_full.sum(axis=1), dtype=np.int64
    ).ravel()
    real_target_common = real_target_full[:, real_cols].tocsr()

    quasi_common = quasi[:, pred_positions].tocsr()
    template_common = template[:, pred_positions].tocsr()
    real_target_rate = (
        np.diff(real_target_common.tocsc().indptr) / real_target_common.shape[0]
    )
    real_control_rate = np.asarray(control["detection_rate_common"])

    quasi_libs = np.asarray(quasi.sum(axis=1), dtype=np.int64).ravel()
    template_libs = np.asarray(template.sum(axis=1), dtype=np.int64).ravel()
    control_libs = np.asarray(control["library_sizes"], dtype=np.int64)
    library_ks = {
        "quasi": float(
            ks_2samp(
                quasi_libs / np.median(quasi_libs),
                control_libs / np.median(control_libs),
            ).statistic
        ),
        "template": float(
            ks_2samp(
                template_libs / np.median(template_libs),
                control_libs / np.median(control_libs),
            ).statistic
        ),
    }

    detection = {
        "quasi": detection_summary(
            quasi_common, real_target_rate, real_control_rate
        ),
        "template": detection_summary(
            template_common, real_target_rate, real_control_rate
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
        "template": cpm_comparison(
            template_mean,
            nominal=nominal,
            template_target=template_target,
            real_target=real_target_mean,
            real_control=real_control_mean,
            gate=gate,
        ),
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
        gene_chunk=gene_chunk,
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

    got_template_cols = np.asarray(
        template.sum(axis=0), dtype=np.int64
    ).ravel()
    quasi_bulk = scorer_pseudobulk(quasi)
    template_bulk = scorer_pseudobulk(template)
    pseudobulk = {
        "column_totals_exact": bool(
            np.array_equal(quasi_cols, got_template_cols)
        ),
        "max_abs_column_total_difference": int(
            np.max(np.abs(quasi_cols - got_template_cols), initial=0)
        ),
        "bulk_lognorm_array_equal": bool(
            np.array_equal(quasi_bulk, template_bulk)
        ),
        "max_abs_bulk_lognorm_difference": float(
            np.max(np.abs(quasi_bulk - template_bulk), initial=0.0)
        ),
    }
    if not pseudobulk["column_totals_exact"]:
        raise RuntimeError(f"{target}: template changed column totals")
    if not pseudobulk["bulk_lognorm_array_equal"]:
        raise RuntimeError(f"{target}: template changed scorer pseudobulk")

    arms = {
        arm: arm_summary(
            library_ks=library_ks[arm],
            detection=detection[arm],
            cpm=cpm[arm],
            de=de[arm],
        )
        for arm in ARM_NAMES
    }
    result = {
        "target": target,
        "cells": {
            "quasi": int(quasi.shape[0]),
            "template": int(template.shape[0]),
            "real_target": int(real_target_common.shape[0]),
            "real_control_mwu_sample": int(de_control_common.shape[0]),
        },
        "tested_genes_control_cpm_gt_5_excluding_target": int(gate.sum()),
        "pseudobulk": pseudobulk,
        "arms": arms,
        "template_emitter_diagnostics": emitter_diag.to_dict(),
    }
    memory_checkpoint(f"{target}: completed", rss_limit)
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
        default="experiments/049-template-emitter-multitarget/results.json",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=list(DEFAULT_TARGETS),
        help="Bounded target list; defaults to the complete three-target overlap.",
    )
    parser.add_argument("--template-cells", type=int, default=400)
    parser.add_argument("--de-control-cells", type=int, default=2000)
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--gene-chunk", type=int, default=128)
    parser.add_argument("--seed", type=int, default=4300)
    parser.add_argument("--max-rss-gib", type=float, default=2.8)
    args = parser.parse_args()

    if len(args.targets) > 5:
        parser.error("bounded experiment permits at most five targets")
    targets = tuple(dict.fromkeys(args.targets))
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
        missing = [t for t in targets if t not in overlap]
        if missing:
            raise RuntimeError(
                f"requested targets missing from prediction/reference overlap: {missing}"
            )

        real_control_idx = np.flatnonzero(real_labels == CONTROL)
        control = stream_control_summaries(
            real,
            real_control_idx,
            real_cols,
            chunk_rows=args.chunk_rows,
            rss_limit=rss_limit,
        )
        control_libs = np.asarray(control["library_sizes"], dtype=np.int64)
        template_local = stratified_positions(control_libs, args.template_cells)
        de_control_local = stratified_positions(
            control_libs, args.de_control_cells
        )
        template_global = np.sort(real_control_idx[template_local])
        de_control_global = np.sort(real_control_idx[de_control_local])

        template_full = read_rows(real, template_global)
        template_full_libs = np.asarray(
            template_full.sum(axis=1), dtype=np.int64
        ).ravel()
        template_aligned = align_real_rows(
            template_full, real_cols, pred_positions, pred.n_vars
        )
        de_control_full = read_rows(real, de_control_global)
        de_control_libs = np.asarray(
            de_control_full.sum(axis=1), dtype=np.int64
        ).ravel()
        de_control_common = de_control_full[:, real_cols].tocsr()
        del template_full, de_control_full
        gc.collect()
        memory_checkpoint("shared controls loaded", rss_limit)

        per_target: dict[str, dict[str, object]] = {}
        for offset, target in enumerate(targets):
            print(f"[target] {target}", flush=True)
            per_target[target] = evaluate_target(
                target=target,
                pred=pred,
                real=real,
                pred_labels=pred_labels,
                real_labels=real_labels,
                pred_genes=pred_genes,
                pred_positions=pred_positions,
                real_cols=real_cols,
                control=control,
                template_aligned=template_aligned,
                template_full_libs=template_full_libs,
                de_control_common=de_control_common,
                de_control_libs=de_control_libs,
                seed=target_seed(args.seed, target, offset),
                gene_chunk=args.gene_chunk,
                rss_limit=rss_limit,
            )
            gc.collect()

        result = {
            "date": "2026-08-22",
            "experiment": "049-template-emitter-multitarget",
            "prediction": str(prediction_path),
            "real": str(real_path),
            "prediction_shape": list(pred.shape),
            "real_shape": list(real.shape),
            "prediction_target_count": len(pred_targets),
            "real_target_count": len(real_targets),
            "complete_overlap_targets": overlap,
            "evaluated_targets": list(targets),
            "common_gene_count": int(pred_positions.size),
            "cells": {
                "real_control_all": int(real_control_idx.size),
                "control_templates": int(args.template_cells),
                "real_control_mwu_sample": int(args.de_control_cells),
            },
            "primary_metric_directions": PRIMARY_METRICS,
            "per_target": per_target,
            "macro": macro_summary(per_target),
            "all_pseudobulk_exact": bool(
                all(
                    target_result["pseudobulk"]["column_totals_exact"]
                    and target_result["pseudobulk"]["bulk_lognorm_array_equal"]
                    for target_result in per_target.values()
                )
            ),
            "peak_rss_gib": peak_rss_bytes() / GIB,
            "elapsed_seconds": time.time() - started,
            "limitations": [
                "The prediction artifact has an H1-test-style 85-target, "
                "18,080-gene layout, while the real perturbations and control "
                "templates are Orion-HCT116 with 132 targets and 18,533 genes.",
                "This is an emission-layout diagnostic on the 18,077 common genes, "
                "not a clean model or cross-context transfer validation.",
                "The experiment-040 artifact inherits the quarantined H1 target-mask "
                "mismatch; comparisons remain informative only because quasi and "
                "template reuse identical per-target column totals.",
                "The DE diagnostic is a two-sided SciPy Mann-Whitney U test against "
                "2,000 library-quantile-stratified HCT116 controls, not cell-eval2.",
                "No cell-eval2 invocation and no challenge submission occurred.",
            ],
        }
        result["peak_rss_gib"] = peak_rss_bytes() / GIB
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(result["macro"], indent=2, sort_keys=True), flush=True)
    finally:
        pred.file.close()
        real.file.close()

    memory_checkpoint("done", rss_limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
