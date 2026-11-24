#!/usr/bin/env python3
"""Compare experiment-077 candidates with the frozen experiment-040 baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}


def scaled(raw: dict[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for metric, (baseline, replicate) in ANCHORS.items():
        value = (raw[metric] - baseline) / (replicate - baseline)
        if metric == "expr_mse_unbiased_capped_norm":
            value = max(0.0, min(1.0, value))
        result[metric] = float(value)
    result["avg6"] = sum(result.values()) / 6
    return result


def read_aggregate(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path).set_index("statistic")
    return {metric: float(frame.loc["mean", metric]) for metric in ANCHORS}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    baseline_raw = read_aggregate(args.baseline)
    baseline_scaled = scaled(baseline_raw)
    runs = json.loads(args.results.read_text())
    if len(runs) != 2 or any("error" in run for run in runs):
        raise ValueError("expected exactly two successful candidate runs")

    candidates = []
    for run in sorted(runs, key=lambda item: item["tag"]):
        raw = {metric: float(run["raw"][metric]) for metric in ANCHORS}
        score = scaled(raw)
        candidates.append(
            {
                "tag": run["tag"],
                "raw": raw,
                "scaled": score,
                "delta_avg6": score["avg6"] - baseline_scaled["avg6"],
                "delta_pds_raw": raw["pds_cosine"] - baseline_raw["pds_cosine"],
                "positive_mse": raw["expr_mse_unbiased_capped_norm"] < 0.989,
            }
        )

    for candidate in candidates:
        candidate["passes_individual_gate"] = bool(
            candidate["positive_mse"]
            and candidate["delta_pds_raw"] >= -0.003
            and candidate["delta_avg6"] >= 0.003
        )
    adjacent_positive = all(item["delta_avg6"] > 0 for item in candidates)
    selected = max(candidates, key=lambda item: item["scaled"]["avg6"])
    decision = bool(selected["passes_individual_gate"] and adjacent_positive)
    report = {
        "baseline": {
            "raw": baseline_raw,
            "scaled": baseline_scaled,
        },
        "candidates": candidates,
        "adjacent_alphas_both_positive": adjacent_positive,
        "selected": selected["tag"],
        "passes_promotion_gate": decision,
        "submission_authorized": False,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
