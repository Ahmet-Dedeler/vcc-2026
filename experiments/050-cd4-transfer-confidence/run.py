"""Audit CD4-only metadata as a target/condition transfer-confidence signal.

The public supplementary table contains direct K562-versus-CD4 effect
correlations for a subset of targets.  We use those correlations only as a
held-out biological diagnostic: candidate confidence features themselves must
come entirely from CD4 (effect breadth, donor reproducibility, cell count, and
on-target QC).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")


def safe_spearman(x: pd.Series, y: pd.Series) -> dict[str, float | int]:
    frame = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 3 or frame["x"].nunique() < 2:
        return {"n": len(frame), "rho": float("nan"), "pvalue": float("nan")}
    result = spearmanr(frame["x"], frame["y"])
    return {
        "n": len(frame),
        "rho": float(result.statistic),
        "pvalue": float(result.pvalue),
    }


def bootstrap_mean_delta(
    selected: np.ndarray,
    baseline: np.ndarray,
    *,
    seed: int,
    replicates: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    delta = np.asarray(selected, dtype=np.float64) - np.asarray(
        baseline, dtype=np.float64
    )
    draws = rng.choice(delta, size=(replicates, len(delta)), replace=True).mean(axis=1)
    return {
        "mean": float(delta.mean()),
        "ci95_low": float(np.quantile(draws, 0.025)),
        "ci95_high": float(np.quantile(draws, 0.975)),
        "positive_probability": float(np.mean(draws > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--comparison",
        default=ROOT / "data/external/cd4t/K562_comparison.suppl_table.csv",
        type=Path,
    )
    parser.add_argument(
        "--de-stats",
        default=ROOT / "data/external/cd4t/DE_stats.suppl_table.csv",
        type=Path,
    )
    parser.add_argument(
        "--panel",
        default=ROOT / "data/challenge/pert_counts.csv",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default=Path(__file__).resolve().parent / "results.json",
        type=Path,
    )
    parser.add_argument(
        "--donor-root",
        default=ROOT / "data/donors",
        type=Path,
    )
    parser.add_argument("--bootstrap", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=50)
    args = parser.parse_args()

    panel = set(pd.read_csv(args.panel)["target_gene"].dropna().astype(str))
    panel_order = pd.read_csv(args.panel)["target_gene"].dropna().astype(str).tolist()
    comparison = pd.read_csv(args.comparison)
    comparison = comparison[
        comparison["target_contrast_gene_name"].astype(str).isin(panel)
    ].copy()
    comparison = comparison[comparison["condition"].isin(CONDITIONS)].copy()

    de_stats = pd.read_csv(args.de_stats)
    de_stats = de_stats[
        de_stats["target_contrast_gene_name"].astype(str).isin(panel)
    ].copy()
    de_stats = de_stats[de_stats["culture_condition"].isin(CONDITIONS)].copy()

    merged = comparison.merge(
        de_stats[
            [
                "target_contrast_gene_name",
                "culture_condition",
                "ontarget_significant",
                "offtarget_flag",
                "n_cells_target",
                "n_downstream",
                "n_total_de_genes",
            ]
        ],
        left_on=["target_contrast_gene_name", "condition"],
        right_on=["target_contrast_gene_name", "culture_condition"],
        how="left",
        validate="one_to_one",
    )
    own_de_column = {
        "Rest": "n_degs_MASH_Rest",
        "Stim8hr": "n_degs_MASH_Stim8hr",
        "Stim48hr": "n_degs_MASH_Stim48hr",
    }
    merged["own_mash_de"] = [
        row[own_de_column[row["condition"]]] for _, row in merged.iterrows()
    ]

    condition_summary: dict[str, dict[str, float | int]] = {}
    feature_correlations: dict[str, dict[str, dict[str, float | int]]] = {}
    features = (
        "own_mash_de",
        "n_downstream",
        "n_total_de_genes",
        "donor_correlation_mean",
        "n_cells_target",
    )
    for condition in CONDITIONS:
        subset = merged[merged["condition"] == condition]
        transfer = subset["logfc_pearson_r"].dropna().to_numpy(dtype=np.float64)
        condition_summary[condition] = {
            "n": len(transfer),
            "mean_transfer_r": float(transfer.mean()),
            "median_transfer_r": float(np.median(transfer)),
            "fraction_positive": float(np.mean(transfer > 0)),
            "fraction_above_0_2": float(np.mean(transfer > 0.2)),
        }
        feature_correlations[condition] = {
            feature: safe_spearman(
                np.log1p(subset[feature])
                if feature in {"own_mash_de", "n_downstream", "n_total_de_genes"}
                else subset[feature],
                subset["logfc_pearson_r"],
            )
            for feature in features
        }

    transfer = merged.pivot(
        index="target_contrast_gene_name",
        columns="condition",
        values="logfc_pearson_r",
    ).dropna()
    breadth = merged.pivot(
        index="target_contrast_gene_name",
        columns="condition",
        values="own_mash_de",
    ).reindex(transfer.index)

    fixed_best_condition = max(
        CONDITIONS, key=lambda condition: float(transfer[condition].mean())
    )
    breadth_choice = breadth.idxmax(axis=1)
    breadth_selected = np.asarray(
        [
            transfer.loc[target, breadth_choice.loc[target]]
            for target in transfer.index
        ],
        dtype=np.float64,
    )
    fixed = transfer[fixed_best_condition].to_numpy(dtype=np.float64)
    oracle = transfer.max(axis=1).to_numpy(dtype=np.float64)

    condition_pairwise = {
        f"{left}__{right}": safe_spearman(transfer[left], transfer[right])
        for i, left in enumerate(CONDITIONS)
        for right in CONDITIONS[i + 1 :]
    }
    choice_counts = {
        condition: int((breadth_choice == condition).sum())
        for condition in CONDITIONS
    }
    donor_names = ("k562_gwps", "orion_hct116", "orion_hek293t")
    donor_targets = {
        name: set((args.donor_root / name / "perts.txt").read_text().splitlines())
        for name in donor_names
    }
    cd4_targets = set(de_stats["target_contrast_gene_name"].astype(str))
    existing_coverage = {
        target: sum(target in donor_targets[name] for name in donor_names)
        for target in panel_order
    }
    with_cd4_coverage = {
        target: count + int(target in cd4_targets)
        for target, count in existing_coverage.items()
    }

    result = {
        "protocol": {
            "panel_targets_in_comparison": int(
                merged["target_contrast_gene_name"].nunique()
            ),
            "rows": len(merged),
            "conditions": list(CONDITIONS),
            "confidence_features_use_cd4_only": True,
            "k562_correlation_used_only_as_diagnostic_truth": True,
        },
        "condition_summary": condition_summary,
        "feature_correlations": feature_correlations,
        "condition_pairwise_transfer_correlation": condition_pairwise,
        "breadth_condition_choice": {
            "fixed_best_condition": fixed_best_condition,
            "choice_counts": choice_counts,
            "fixed_mean_transfer_r": float(fixed.mean()),
            "breadth_selected_mean_transfer_r": float(breadth_selected.mean()),
            "oracle_mean_transfer_r": float(oracle.mean()),
            "breadth_vs_fixed": bootstrap_mean_delta(
                breadth_selected,
                fixed,
                seed=args.seed,
                replicates=args.bootstrap,
            ),
            "oracle_vs_fixed": bootstrap_mean_delta(
                oracle,
                fixed,
                seed=args.seed + 1,
                replicates=args.bootstrap,
            ),
        },
        "panel_coverage": {
            "existing_count_histogram": {
                str(count): int(sum(value == count for value in existing_coverage.values()))
                for count in sorted(set(existing_coverage.values()))
            },
            "with_cd4_count_histogram": {
                str(count): int(sum(value == count for value in with_cd4_coverage.values()))
                for count in sorted(set(with_cd4_coverage.values()))
            },
            "cd4_present_k562_missing": [
                target
                for target in panel_order
                if target in cd4_targets and target not in donor_targets["k562_gwps"]
            ],
            "missing_from_all_existing": [
                target
                for target in panel_order
                if not any(target in donor_targets[name] for name in donor_names)
            ],
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
