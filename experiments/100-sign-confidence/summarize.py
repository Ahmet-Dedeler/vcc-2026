#!/usr/bin/env python3
"""Aggregate experiment-100 surface results with target-block uncertainty."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SURFACES = ("h1", "hct116", "hek293t")
N_BOOT = 2_000


def seed(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], "little")


def macro_boot(values: dict[str, np.ndarray], label: str) -> dict[str, float]:
    rng = np.random.default_rng(seed(label))
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        means = []
        for surface in SURFACES:
            x = values[surface]
            means.append(x[rng.integers(0, len(x), len(x))].mean())
        boot[i] = np.mean(means)
    estimate = float(np.mean([x.mean() for x in values.values()]))
    return {
        "estimate": estimate,
        "lo": float(np.percentile(boot, 2.5)),
        "hi": float(np.percentile(boot, 97.5)),
    }


def main() -> None:
    metrics = {
        surface: json.loads((WORK / surface / "metrics.json").read_text())
        for surface in SURFACES
    }
    rows = {
        surface: np.load(WORK / surface / "per_target_metrics.npz")
        for surface in SURFACES
    }
    comparisons = {
        "rank8_precision_minus_agreement": {
            s: rows[s]["rank8_precision"] - rows[s]["agreement_precision"]
            for s in SURFACES
        },
        "rank16_precision_minus_agreement": {
            s: rows[s]["rank16_precision"] - rows[s]["agreement_precision"]
            for s in SURFACES
        },
        "rank8_linear_credit_minus_identity": {
            s: rows[s]["rank8_credit"] - rows[s]["identity_credit"]
            for s in SURFACES
        },
        "rank8_linear_credit_minus_agreement": {
            s: rows[s]["rank8_credit"] - rows[s]["agreement_credit"]
            for s in SURFACES
        },
        "rank16_linear_credit_minus_agreement": {
            s: rows[s]["rank16_credit"] - rows[s]["agreement_credit"]
            for s in SURFACES
        },
        "rank16_linear_credit_minus_identity": {
            s: rows[s]["rank16_credit"] - rows[s]["identity_credit"]
            for s in SURFACES
        },
        "rank8_pds_minus_identity": {
            s: rows[s]["pds_rank8"] - rows[s]["pds_identity"]
            for s in SURFACES
        },
        "rank16_pds_minus_identity": {
            s: rows[s]["pds_rank16"] - rows[s]["pds_identity"]
            for s in SURFACES
        },
        "rank8_pds_minus_agreement": {
            s: rows[s]["pds_rank8"] - rows[s]["pds_agreement"]
            for s in SURFACES
        },
        "rank16_pds_minus_agreement": {
            s: rows[s]["pds_rank16"] - rows[s]["pds_agreement"]
            for s in SURFACES
        },
    }
    bootstrap = {
        name: macro_boot(values, name) for name, values in comparisons.items()
    }
    observed = bootstrap["rank8_precision_minus_agreement"]["estimate"]
    null = np.mean(
        [
            np.asarray(
                metrics[s]["null_rank8_precision_gain_over_agreement"]["values"],
                dtype=np.float64,
            )
            for s in SURFACES
        ],
        axis=0,
    )
    null_p95 = float(np.percentile(null, 95))
    summary = {
        "surfaces": {
            s: {
                "targets_total": metrics[s]["targets_total"],
                "targets_de_scored": metrics[s]["targets_de_scored"],
                "weighted_sign_precision": metrics[s]["weighted_sign_precision"],
                "linear_nmae_proxy": metrics[s]["linear_nmae_proxy"],
                "pds": metrics[s]["pds"],
            }
            for s in SURFACES
        },
        "macro_bootstrap": bootstrap,
        "rank8_null": {
            "observed_precision_gain_over_agreement": observed,
            "null_values": null.tolist(),
            "null_p50": float(np.median(null)),
            "null_p95": null_p95,
            "empirical_p": float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        },
        "decision": {
            "weighted_precision_nonnegative_each_surface": all(
                metrics[s]["bootstrap"]["rank8_precision_minus_agreement"][
                    "estimate"
                ]
                >= 0
                for s in SURFACES
            ),
            "linear_proxy_nonnegative_each_surface": all(
                metrics[s]["bootstrap"]["rank8_linear_credit_minus_identity"][
                    "estimate"
                ]
                >= 0
                for s in SURFACES
            ),
            "pds_nonnegative_each_surface": all(
                metrics[s]["bootstrap"]["rank8_pds_minus_identity"]["estimate"]
                >= 0
                for s in SURFACES
            ),
            "macro_precision_ci_lower_positive": (
                bootstrap["rank8_precision_minus_agreement"]["lo"] > 0
            ),
            "observed_above_null_p95": observed > null_p95,
            "rank16_nonnegative_each_surface": all(
                metrics[s]["bootstrap"]["rank16_precision_minus_agreement"][
                    "estimate"
                ]
                >= 0
                for s in SURFACES
            ),
        },
    }
    summary["decision"]["go"] = all(summary["decision"].values())
    (WORK / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
