#!/usr/bin/env python3
"""Collect paired quasi/template cell-eval2 means from experiment 053."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

LOWER_IS_BETTER = {
    "de_wilcoxon_lfc_nmae",
    "expr_distance_unbiased",
    "expr_mse_unbiased",
    "expr_mse_unbiased_capped",
    "expr_mse_unbiased_capped_norm",
}


def read_mean(path: Path) -> dict[str, float]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle)
        for row in rows:
            if row.get("statistic") == "mean":
                result = {}
                for key, value in row.items():
                    if key == "statistic" or value in (None, "", "NaN"):
                        continue
                    number = float(value)
                    if math.isfinite(number):
                        result[key] = number
                return result
    raise ValueError(f"{path} has no mean row")


def collect(root: Path) -> dict[str, object]:
    paired: dict[str, object] = {}
    for tag_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        quasi_path = tag_dir / "quasi-agg_results.csv"
        template_path = tag_dir / "template-run" / "agg_results.csv"
        if not quasi_path.is_file() or not template_path.is_file():
            continue
        quasi = read_mean(quasi_path)
        template = read_mean(template_path)
        common = sorted(set(quasi) & set(template))
        comparison = {}
        for metric in common:
            delta = template[metric] - quasi[metric]
            comparison[metric] = {
                "quasi": quasi[metric],
                "template": template[metric],
                "template_minus_quasi": delta,
                "favorable_delta_for_template": (
                    -delta if metric in LOWER_IS_BETTER else delta
                ),
                "preferred": ("lower" if metric in LOWER_IS_BETTER else "higher"),
            }
        paired[tag_dir.name] = comparison
    return {
        "paired_gamma_runs": len(paired),
        "runs": paired,
        "note": (
            "Raw metric comparison only. Reference-scaled overall scores should "
            "be computed with the same real bundle before choosing a candidate."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = collect(args.root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
