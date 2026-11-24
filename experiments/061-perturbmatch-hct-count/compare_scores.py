#!/usr/bin/env python3
"""Compare experiment-061 candidate and existing clean HCT baseline scores."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

METRICS = (
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
    "de_wilcoxon_lfc_nmae",
)
LOWER_IS_BETTER = {
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
}
# Midpoints of the three official 2026 context anchors, matching experiment 031/040.
ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def read_mean(path: Path) -> dict[str, float]:
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("statistic") != "mean":
                continue
            values: dict[str, float] = {}
            for metric in METRICS:
                value = float(row[metric])
                if not math.isfinite(value):
                    raise ValueError(f"{path}: non-finite {metric}")
                values[metric] = value
            return values
    raise ValueError(f"{path}: no mean row")


def scaled(metric: str, value: float) -> float:
    baseline, replicate = ANCHORS[metric]
    result = (value - baseline) / (replicate - baseline)
    if metric == "expr_mse_unbiased_capped_norm":
        return max(0.0, min(1.0, result))
    return result


def compare(baseline_path: Path, candidate_path: Path) -> dict[str, object]:
    baseline_raw = read_mean(baseline_path)
    candidate_raw = read_mean(candidate_path)
    metrics: dict[str, object] = {}
    baseline_scaled_values = []
    candidate_scaled_values = []
    for metric in METRICS:
        baseline_scaled = scaled(metric, baseline_raw[metric])
        candidate_scaled = scaled(metric, candidate_raw[metric])
        baseline_scaled_values.append(baseline_scaled)
        candidate_scaled_values.append(candidate_scaled)
        raw_delta = candidate_raw[metric] - baseline_raw[metric]
        metrics[metric] = {
            "preferred": "lower" if metric in LOWER_IS_BETTER else "higher",
            "baseline_raw": baseline_raw[metric],
            "candidate_raw": candidate_raw[metric],
            "raw_candidate_minus_baseline": raw_delta,
            "raw_favorable_delta": (
                -raw_delta if metric in LOWER_IS_BETTER else raw_delta
            ),
            "baseline_scaled": baseline_scaled,
            "candidate_scaled": candidate_scaled,
            "scaled_candidate_minus_baseline": candidate_scaled - baseline_scaled,
        }
    baseline_avg6 = sum(baseline_scaled_values) / len(METRICS)
    candidate_avg6 = sum(candidate_scaled_values) / len(METRICS)
    return {
        "experiment": "061-perturbmatch-hct-count",
        "status": "scored",
        "comparison": "candidate minus existing clean gamma-0.65 HCT baseline",
        "scaling": (
            "reference scaling uses midpoint of official 2026 A/B/C baseline "
            "and replicate anchors, exactly as experiments 031/040"
        ),
        "metrics": metrics,
        "avg6": {
            "baseline": baseline_avg6,
            "candidate": candidate_avg6,
            "candidate_minus_baseline": candidate_avg6 - baseline_avg6,
        },
        "score_files": {
            "baseline": {
                "path": str(baseline_path),
                "sha256": file_sha256(baseline_path),
            },
            "candidate": {
                "path": str(candidate_path),
                "sha256": file_sha256(candidate_path),
            },
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = compare(args.baseline, args.candidate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

