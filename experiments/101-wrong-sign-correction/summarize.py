#!/usr/bin/env python3
"""Aggregate experiment-101 results and write RESULTS.md."""

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
    for draw in range(N_BOOT):
        means = []
        for surface in SURFACES:
            x = np.asarray(values[surface], dtype=np.float64)
            means.append(x[rng.integers(0, len(x), len(x))].mean())
        boot[draw] = np.mean(means)
    estimate = float(np.mean([np.mean(x) for x in values.values()]))
    return {
        "estimate": estimate,
        "lo": float(np.percentile(boot, 2.5)),
        "hi": float(np.percentile(boot, 97.5)),
    }


def macro_ratio_boot(
    numerator: dict[str, np.ndarray],
    denominator: dict[str, np.ndarray],
    label: str,
) -> dict[str, float | None]:
    estimate_parts = []
    for surface in SURFACES:
        d = float(denominator[surface].sum())
        estimate_parts.append(float(numerator[surface].sum() / d) if d > 0 else np.nan)
    estimate = float(np.nanmean(estimate_parts))
    rng = np.random.default_rng(seed(label))
    boot = []
    for _ in range(N_BOOT):
        parts = []
        for surface in SURFACES:
            n = numerator[surface]
            d = denominator[surface]
            idx = rng.integers(0, len(n), len(n))
            total = float(d[idx].sum())
            if total > 0:
                parts.append(float(n[idx].sum() / total))
        if parts:
            boot.append(float(np.mean(parts)))
    if not boot:
        return {"estimate": None, "lo": None, "hi": None}
    return {
        "estimate": estimate,
        "lo": float(np.percentile(boot, 2.5)),
        "hi": float(np.percentile(boot, 97.5)),
    }


def fmt(value: float | None, digits: int = 5) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def main() -> None:
    metrics = {
        surface: json.loads((WORK / surface / "metrics.json").read_text())
        for surface in SURFACES
    }
    rows = {
        surface: np.load(WORK / surface / "per_target_metrics.npz", allow_pickle=False)
        for surface in SURFACES
    }
    comparisons = {
        "flip_credit_minus_identity": {
            s: rows[s]["flip_credit"] - rows[s]["identity_credit"] for s in SURFACES
        },
        "zero_credit_minus_identity": {
            s: rows[s]["zero_credit"] - rows[s]["identity_credit"] for s in SURFACES
        },
        "flip_credit_minus_zero": {
            s: rows[s]["flip_credit"] - rows[s]["zero_credit"] for s in SURFACES
        },
        "flip_pds_minus_identity": {
            s: rows[s]["pds_flip"] - rows[s]["pds_identity"] for s in SURFACES
        },
        "zero_pds_minus_identity": {
            s: rows[s]["pds_zero"] - rows[s]["pds_identity"] for s in SURFACES
        },
    }
    macro_bootstrap = {
        name: macro_boot(values, name) for name, values in comparisons.items()
    }
    selected_tail_surfaces = [
        surface
        for surface in SURFACES
        if float(rows[surface]["tail_selected_weight"].sum()) > 0
    ]
    tail_error = macro_ratio_boot(
        {s: rows[s]["tail_wrong_weight"] for s in SURFACES},
        {s: rows[s]["tail_selected_weight"] for s in SURFACES},
        "selected-surface-tail-error",
    )
    null = np.mean(
        [
            np.asarray(
                metrics[s]["permutation_null"]["linear_credit_gain"],
                dtype=np.float64,
            )
            for s in SURFACES
        ],
        axis=0,
    )
    observed = macro_bootstrap["flip_credit_minus_identity"]["estimate"]
    null_p95 = float(np.percentile(null, 95, method="higher"))
    decision = {
        "linear_gain_at_least_0.003_each_surface": all(
            metrics[s]["linear_credit_gain"]["flip_minus_identity"] >= 0.003
            for s in SURFACES
        ),
        "macro_linear_bootstrap_lower_positive": (
            macro_bootstrap["flip_credit_minus_identity"]["lo"] > 0
        ),
        "tail_error_above_0.60_each_surface": all(
            (
                metrics[s]["flipped_tail"]["weighted_error_rate"]["estimate"]
                is not None
                and metrics[s]["flipped_tail"]["weighted_error_rate"]["estimate"] > 0.60
            )
            for s in SURFACES
        ),
        "pds_nonnegative_each_surface": all(
            metrics[s]["pds_delta"]["flip_minus_identity"] >= 0 for s in SURFACES
        ),
        "observed_above_permutation_p95": observed > null_p95,
        "flip_beats_zero_macro": (
            macro_bootstrap["flip_credit_minus_zero"]["estimate"] > 0
        ),
    }
    decision["go"] = all(decision.values())
    summary = {
        "surfaces": {
            surface: {
                "selected": metrics[surface]["selected"],
                "targets_total": metrics[surface]["targets_total"],
                "targets_de_scored": metrics[surface]["targets_de_scored"],
                "linear_nmae_proxy": metrics[surface]["linear_nmae_proxy"],
                "linear_credit_gain": metrics[surface]["linear_credit_gain"],
                "pds": metrics[surface]["pds"],
                "pds_delta": metrics[surface]["pds_delta"],
                "flipped_tail": metrics[surface]["flipped_tail"],
            }
            for surface in SURFACES
        },
        "macro_bootstrap": macro_bootstrap,
        "selected_surface_tail_error": {
            **tail_error,
            "surfaces": selected_tail_surfaces,
            "n_surfaces": len(selected_tail_surfaces),
        },
        "permutation_null": {
            "observed_macro_linear_credit_gain": observed,
            "values": null.tolist(),
            "p50": float(np.median(null)),
            "p95": null_p95,
            "empirical_p": float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        },
        "decision": decision,
    }
    json_dump = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    (WORK / "summary.json").write_text(json_dump)

    verdict = "GO" if decision["go"] else "NO-GO"
    lines = [
        "# 101 — cross-fitted selective wrong-sign correction",
        "",
        f"**Verdict: {verdict}.**",
        "",
        "The experiment used the preregistered rank-16 signed ridge grid, "
        "strict outer-donor holdout, four whole-target folds, 2,000 "
        "target-block bootstrap draws, and 20 complete target-label "
        "permutation nulls. No count-level scorer, AWS worker, challenge "
        "artifact, submission, or commit was used.",
        "",
        "## Selected configurations and outer results",
        "",
        "| surface | selected | linear NMAE identity | flip | zero | flip credit gain | PDS delta | tail error | selected pairs |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for surface in SURFACES:
        item = metrics[surface]
        selected = item["selected"]
        selected_text = (
            "identity"
            if selected == "identity"
            else (
                f"r={selected['ridge']:g}, tau={selected['tau']:.2f}, "
                f"lambda={selected['lambda']:g}"
            )
        )
        tail = item["flipped_tail"]["weighted_error_rate"]["estimate"]
        lines.append(
            f"| {surface.upper()} | {selected_text} | "
            f"{item['linear_nmae_proxy']['identity']:.5f} | "
            f"{item['linear_nmae_proxy']['flip']:.5f} | "
            f"{item['linear_nmae_proxy']['zero']:.5f} | "
            f"{item['linear_credit_gain']['flip_minus_identity']:+.5f} | "
            f"{item['pds_delta']['flip_minus_identity']:+.5f} | "
            f"{fmt(tail)} | "
            f"{item['flipped_tail']['selected_pairs']} |"
        )
    lines.extend(
        [
            "",
            "## Macro uncertainty and null",
            "",
            f"- Flip linear-credit gain: "
            f"{observed:+.5f} "
            f"[{macro_bootstrap['flip_credit_minus_identity']['lo']:+.5f}, "
            f"{macro_bootstrap['flip_credit_minus_identity']['hi']:+.5f}].",
            f"- Zero-tail linear-credit gain: "
            f"{macro_bootstrap['zero_credit_minus_identity']['estimate']:+.5f}.",
            f"- Flip minus zero: "
            f"{macro_bootstrap['flip_credit_minus_zero']['estimate']:+.5f} "
            f"[{macro_bootstrap['flip_credit_minus_zero']['lo']:+.5f}, "
            f"{macro_bootstrap['flip_credit_minus_zero']['hi']:+.5f}].",
            f"- Flip PDS delta: "
            f"{macro_bootstrap['flip_pds_minus_identity']['estimate']:+.5f} "
            f"[{macro_bootstrap['flip_pds_minus_identity']['lo']:+.5f}, "
            f"{macro_bootstrap['flip_pds_minus_identity']['hi']:+.5f}].",
            f"- Selected-surface weighted tail error "
            f"({len(selected_tail_surfaces)}/{len(SURFACES)} surfaces: "
            f"{', '.join(selected_tail_surfaces)}): "
            f"{fmt(tail_error['estimate'])} "
            f"[{fmt(tail_error['lo'])}, {fmt(tail_error['hi'])}].",
            f"- Permutation null p50 / p95: "
            f"{np.median(null):+.5f} / {null_p95:+.5f}; "
            f"empirical one-sided p={summary['permutation_null']['empirical_p']:.4f}.",
            "",
            "## Preregistered decision",
            "",
            "| condition | result |",
            "| --- | --- |",
        ]
    )
    labels = {
        "linear_gain_at_least_0.003_each_surface": "linear gain >= 0.003 on every surface",
        "macro_linear_bootstrap_lower_positive": "macro bootstrap lower bound above zero",
        "tail_error_above_0.60_each_surface": "tail error > 60% on every surface",
        "pds_nonnegative_each_surface": "PDS nonnegative on every surface",
        "observed_above_permutation_p95": "observed gain above permutation p95",
        "flip_beats_zero_macro": "flip beats matched zero-tail control",
    }
    for key, label in labels.items():
        lines.append(f"| {label} | {'PASS' if decision[key] else '**FAIL**'} |")
    lines.extend(
        [
            f"| overall | **{verdict}** |",
            "",
            "## Interpretation",
            "",
            "The proposed mechanism is not supported. H1 selected only 26 "
            "pairs across 10 targets, and just 44.8% of their metric weight "
            "was actually wrong-sign. In other words, the selected tail was "
            "more often correct than wrong, so reversal was not a reliable "
            "repair. HCT116 and HEK293T found no nested PDS-safe candidate "
            "with positive signed-credit gain and therefore selected identity.",
            "",
            "The tiny positive H1 result is unresolved: the macro gain is "
            "roughly 270 times smaller than the preregistered 0.003 "
            "per-surface threshold, its bootstrap interval crosses zero, and "
            "it falls below the 20-draw permutation-null p95. Although flip "
            "is numerically above the matched zero control, that difference "
            "is also tiny and its interval crosses zero.",
            "",
            "See `work/summary.json` for the complete machine-readable result "
            "and each surface's `selection.json`, `freeze-audit.json`, "
            "`metrics.json`, and `per_target_metrics.npz` for the full "
            "selection, isolation, bootstrap, and null evidence.",
            "",
            "Inner selection intentionally used an unwhitened, "
            "median-norm-scaled base because some available precomputed "
            "whiteners contain the outer held donor. Outer evaluation retained "
            "experiment 100's exact surface-specific whitened base.",
            "",
            "## Provenance limitation",
            "",
            "The freeze audit recorded hashes for predictions, selections, "
            "and null-model arrays, but the scorer did not enforce those "
            "hashes before opening held truth and did not freeze source-code "
            "or whitener hashes. `run.py` was changed after H1/HCT116 scoring "
            "to repair the HEK293T identity-null crash. The change did not "
            "alter the selected H1 model or either identity result, and every "
            "stored number reproduces from the frozen arrays, but code "
            "immutability across all three score stages cannot be proven "
            "retroactively. This weakens provenance, not the NO-GO: the "
            "corrected conservative null threshold is larger and the result "
            "still fails every substantive promotion gate.",
            "",
            "Each permutation draw also creates one deterministic mapping per "
            "whole-target fold after that fold's targets are removed, rather "
            "than one global mapping reused across folds. This is leakage-safe "
            "but should be made explicit in any reuse.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "uv run python experiments/101-wrong-sign-correction/test_run.py",
            "uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py prepare",
            "uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py freeze --surface h1",
            "uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py score --surface h1",
            "# Repeat freeze/score sequentially for hct116 and hek293t.",
            "uv run python experiments/101-wrong-sign-correction/summarize.py",
            "```",
            "",
        ]
    )
    (HERE / "RESULTS.md").write_text("\n".join(lines))
    print(json_dump)


if __name__ == "__main__":
    main()
