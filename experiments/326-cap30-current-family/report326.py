"""Frozen paired cap-dose/interaction report; reuses the tested 323 MSE formula."""

import json

import numpy as np
import pandas as pd
import polars as pl

import run326 as A
import report323 as B
from cell_eval2.metrics.delta import _across_pert_budget
from cell_eval2.prep import bulk_lognorm_means

CONTRASTS = {
    "cap30-minus-cap15": {"cap30": 1, "cap15": -1},
    "cap-dose-minus-plain-dose": {"cap30": 1, "cap15": -1, "plain30": -1, "plain15": 1},
    "cap30-minus-shuffled_cap30": {"cap30": 1, "shuffled_cap30": -1},
}


def linear_contrast(values, weights):
    return sum(values[k] * v for k, v in weights.items())


def main():
    A.verify()
    if (A.WORK / "results.json").exists():
        raise FileExistsError("Do not overwrite a terminal report")
    if not json.loads((A.WORK / "known-answer.json").read_text())["pass"]:
        raise AssertionError("cap15 known-answer has not passed")
    r = np.load(A.BASE / "work/reference.npz", allow_pickle=False)
    targets, genes = r["targets"], r["genes"]
    t_all = bulk_lognorm_means(r["sums"].astype(float), A.E.TS)
    truth, control = t_all[:-1], t_all[-1]
    positions = {g: i for i, g in enumerate(genes)}
    own = np.array([positions[t] for t in targets])
    n = len(targets)
    rng = np.random.default_rng(20260909)
    draws = rng.integers(0, n, (10000, n))
    counts = np.zeros((10000, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(10000), n), draws.ravel()), 1)
    counts = np.vstack((np.ones(n, dtype=np.int64), counts))
    members = list(A.E.MEMBERS)
    if members != ["pds", "mse", "nmae", "fid", "reach", "jac"]:
        raise ValueError("scaling function member order changed")
    anchors = json.loads((A.BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])
    raw_by_arm, overall_by_arm, records, gates = {}, {}, [], {}
    for name in ("cap15", "cap30", "plain15", "plain30", "shuffled_cap30"):
        folder = A.BASE / "work" if name.startswith("plain") else A.WORK
        metrics = json.loads((folder / f"metrics-{name}.json").read_text())
        stats_path = folder / f"stats-{name}.npz"
        A.require_same(f"{name} statistics SHA", A.R.sha(stats_path), metrics["statistics_sha256"])
        with np.load(stats_path, allow_pickle=False) as stats:
            if not np.array_equal(stats["targets"], targets):
                raise ValueError(f"{name}: target order changed")
            pm = bulk_lognorm_means(stats["sums"].astype(float), A.E.TS)
            jk = stats["jk"]
        mse, budget, rho, denominator = B.mse_bootstrap(
            counts, pm, truth, control, jk, r["jk"][:-1], r["jk"][-1], own)
        expected_budget = _across_pert_budget(
            pm, dict(enumerate(map(int, own))), np.full(n, 1/(len(genes)-1)))
        mse_error = abs(mse[0]-metrics["raw"]["mse"])
        budget_error = abs(budget[0]-expected_budget)
        if mse_error > 1e-10 or budget_error > 1e-10 or np.any(denominator <= 0):
            raise AssertionError(f"{name}: MSE identity or denominator failed")
        raw = np.empty((len(counts), len(members)))
        table = pl.read_csv(folder / f"per-target-{name}.csv")
        for j, member in enumerate(members):
            if member == "mse":
                raw[:, j] = mse
            else:
                rows = table.filter(pl.col("metric") == A.E.MEMBERS[member])
                mapping = dict(zip(rows["perturbation"], rows["value"]))
                values = np.array([mapping.get(str(t), np.nan) for t in targets])
                mask = np.isfinite(values)
                denominator_de = counts[:, mask].sum(1)
                if np.any(denominator_de <= 0):
                    raise ValueError(f"{name}/{member}: empty resampled cohort")
                raw[:, j] = (counts[:, mask] @ values[mask]) / denominator_de
            if abs(raw[0, j] - metrics["raw"][member]) > 1e-12:
                raise AssertionError(f"{name}/{member}: aggregate identity failed")
        if not np.isfinite(raw).all():
            raise ValueError(f"{name}: nonfinite resample")
        scaled = B.scaled(raw, base, replicate)
        raw_by_arm[name] = raw
        overall_by_arm[name] = scaled.mean(1)
        records.append({"arm": name, **metrics["raw"],
                        "local_h1_overall": float(scaled[0].mean())})
        gates[name] = {"mse_identity_abs_diff": float(mse_error),
                       "budget_identity_abs_diff": float(budget_error),
                       "rho_bootstrap_range": [float(rho[1:].min()), float(rho[1:].max())]}
    comparisons = {}
    for name, weights in CONTRASTS.items():
        overall = linear_contrast(overall_by_arm, weights)
        raw = linear_contrast(raw_by_arm, weights)
        interval = B.ci(overall[1:])
        comparisons[name] = {
            "local_h1_overall_delta": float(overall[0]), "overall_ci": interval,
            "overall_interval_excludes_zero": bool(interval["lo"] > 0 or interval["hi"] < 0),
            "raw_deltas": {m: {"point": float(raw[0, j]), **B.ci(raw[1:, j])}
                           for j, m in enumerate(members)},
        }
    negative_pass = comparisons["cap30-minus-shuffled_cap30"]["raw_deltas"]["pds"]["lo"] > 0
    result = {
        "status": "COMPLETE_DESCRIPTIVE_ONLY" if negative_pass else "NEGATIVE_CONTROL_FAILED",
        "scope": "one reused H1 context, union OFF; not an A/B/C prediction or submission clearance",
        "negative_control_pass": negative_pass,
        "bootstrap": {"draws": 10000, "seed": 20260909, "unit": "target",
                      "fixed": ["reference gallery", "local H1 anchors", "constructed models"],
                      "mse": "panel ratio and global claim cap recomputed"},
        "arms": records, "comparisons": comparisons, "identity_gates": gates,
        "input_manifest_sha256": A.R.sha(A.HERE / "inputs.json"),
        "anchors_sha256": A.R.sha(A.BASE / "work/h1-anchors.json"),
    }
    A.R.write_json(A.WORK / "results.json", result)
    pd.DataFrame(records).to_csv(A.WORK / "scores.csv", index=False)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
