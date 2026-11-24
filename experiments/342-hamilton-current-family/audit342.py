#!/usr/bin/env python3
"""Independent replay for 342 (optional, after run342.py run).

Rescore every arm from its saved sufficient statistics, re-derive every contrast
interval with a dense gene-wise weighted MSE budget instead of the report's Gram
shortcut, and re-check the 0/1 totals difference and mass identity between arms
from the saved per-gene totals alone. Does not regenerate counts.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import run323 as R            # noqa: E402
import evaluate323 as E       # noqa: E402
import report342 as P         # noqa: E402
from cell_eval2.metrics.delta import PRED_TRACE_CAP_K  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means         # noqa: E402


def dense_mse(weights, p, t, ctrl, pjk, rjk, cjk, own):
    """Direct gene-wise weighted sums, no Gram/target-patch shortcut (337's audit form)."""
    n, g = p.shape
    a = weights.astype(float) / (g - 1)
    total = a.sum(1)
    sums = a @ p
    denominators = np.broadcast_to(total[:, None], sums.shape).copy()
    for i, col in enumerate(own):
        sums[:, col] -= a[:, i] * p[i, col]
        denominators[:, col] -= a[:, i]
    centering = np.divide(sums * sums, denominators, out=np.zeros_like(sums),
                          where=denominators > 0).sum(1)
    square = p * p
    square[np.arange(n), own] = 0
    budget = np.maximum(0, a @ square.sum(1) - centering)
    claim = a @ np.minimum(pjk, PRED_TRACE_CAP_K * rjk)
    rho = np.minimum(1, np.divide(budget, claim, out=np.ones_like(claim), where=claim > 0))
    error = (p - t) ** 2
    error[np.arange(n), own] = 0
    energy = (t - ctrl) ** 2
    energy[np.arange(n), own] = 0
    num = a @ (error.sum(1) - rjk) - rho * claim
    den = a @ (energy.sum(1) - rjk - cjk)
    if np.any(den <= 0):
        raise ValueError("nonpositive denominator")
    return num / den


def main() -> None:
    out = WORK / "independent-audit.json"
    if out.exists():
        raise FileExistsError(out)
    result = json.loads((WORK / "results.json").read_text())
    r, c, _ = E.load_reference()
    reference = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = reference[:-1], reference[-1]
    lookup = {g: i for i, g in enumerate(r["genes"])}
    own = np.array([lookup[t] for t in r["targets"]])
    targets = list(r["targets"])
    members = list(E.MEMBERS)
    anchors = json.loads((BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])
    weights = P.resample_counts(len(targets))
    raw_by_arm, overall, replay, sums_by_arm = {}, {}, {}, {}
    for name in P.ORDER:
        record = json.loads((WORK / f"metrics-{name}.json").read_text())
        path = WORK / f"stats-{name}.npz"
        if R.sha(path) != record["statistics_sha256"]:
            raise ValueError(f"statistics changed: {name}")
        with np.load(path, allow_pickle=False) as z:
            stats = {k: z[k] for k in z.files}
        recomputed, rows, _ = E.score(r, c, stats["sums"], stats["jk"], stats["p"],
                                      stats["q"], stats["lfc"])
        replay[name] = {"six_metric_max_abs": max(abs(recomputed[m] - record["raw"][m])
                                                  for m in members)}
        if replay[name]["six_metric_max_abs"] > 1e-12:
            raise AssertionError(f"metric replay drift: {name}")
        saved = list(csv.DictReader((WORK / f"per-target-{name}.csv").open()))
        index = {(x["perturbation"], x["metric"]): float(x["value"]) for x in saved}
        differences = []
        for row in rows:
            old, new = index[row["perturbation"], row["metric"]], row["value"]
            if np.isnan(old) and np.isnan(new):
                continue
            differences.append(abs(old - new))
        replay[name]["per_target_max_abs"] = max(differences)
        if replay[name]["per_target_max_abs"] > 1e-12:
            raise AssertionError(f"per-target replay drift: {name}")
        p = bulk_lognorm_means(stats["sums"].astype(float), E.TS)
        raw = np.empty((len(weights), len(members)))
        for start in range(0, len(weights), 64):
            raw[start:start + 64, 1] = dense_mse(weights[start:start + 64], p, truth, control,
                                                 stats["jk"], r["jk"][:-1], r["jk"][-1], own)
        for j, m in enumerate(members):
            if m == "mse":
                continue
            values = np.array([index.get((str(t), E.MEMBERS[m]), np.nan) for t in targets])
            valid = np.isfinite(values)
            raw[:, j] = (weights[:, valid] @ values[valid]) / weights[:, valid].sum(1)
        if np.max(np.abs(raw[0] - [record["raw"][m] for m in members])) > 1e-9:
            raise AssertionError(f"point replay drift: {name}")
        scaled = (raw - base) / (replicate - base)
        scaled[:, 1] = np.clip(scaled[:, 1], 0, 1)
        scaled[:, 2] = np.maximum(scaled[:, 2], -6)
        raw_by_arm[name], overall[name] = raw, scaled.mean(1)
        sums_by_arm[name] = stats["sums"].astype(np.int64)
        print(f"audited {name}", flush=True)
    contrast_max = {}
    for key, expected in result["contrasts"].items():
        a, b = key.split("-minus-")
        delta = overall[a] - overall[b]
        vals = [abs(float(delta[0]) - expected["overall"]),
                *np.abs(np.quantile(delta[1:], [.025, .975])
                        - [expected["overall_ci"]["lo"], expected["overall_ci"]["hi"]])]
        for j, m in enumerate(members):
            d = raw_by_arm[a][:, j] - raw_by_arm[b][:, j]
            e = expected["raw"][m]
            vals.extend([abs(d[0] - e["point"]),
                         *np.abs(np.quantile(d[1:], [.025, .975]) - [e["lo"], e["hi"]])])
        contrast_max[key] = float(max(vals))
        if contrast_max[key] > 1e-9:
            raise AssertionError(f"contrast replay drift: {key} {contrast_max[key]}")
    totals = {}
    for name in ("hamilton", "random_up"):
        delta = sums_by_arm[name] - sums_by_arm["floor"]
        rounding = json.loads((WORK / f"rounding-{name}.json").read_text())
        k = np.array([t["residual_k"] for t in rounding["targets"]])
        totals[name] = {
            "zero_one_differences": bool(np.all((delta == 0) | (delta == 1))),
            "per_target_mass_equals_residual": bool(np.array_equal(delta.sum(1), k)),
            "targets_with_any_move": int((delta != 0).any(1).sum())}
        if not (totals[name]["zero_one_differences"]
                and totals[name]["per_target_mass_equals_residual"]):
            raise AssertionError(f"saved totals violate the rounding contract: {name}")
    R.write_json(out, {
        "status": "PASS", "replay": replay, "independent_bootstrap_max_abs": contrast_max,
        "totals_contract_from_saved_statistics": totals,
        "result_sha256": R.sha(WORK / "results.json"),
        "scope": "sufficient-statistic scorer replay, dense weighted MSE bootstrap and the "
                 "0/1 totals contract from saved per-gene totals",
        "not_reproduced": "raw count regeneration"})
    print(json.dumps({"status": "PASS", "independent_bootstrap_max_abs": contrast_max},
                     indent=2))


if __name__ == "__main__":
    main()
