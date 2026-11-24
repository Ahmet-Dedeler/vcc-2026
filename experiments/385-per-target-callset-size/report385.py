#!/usr/bin/env python3
"""385 step 3c: paired-target bootstrap contrasts vs constant-K cap15.

Mirrors 380's report380.py (paired 10,000-target resamples, seed 20260909,
323's frozen H1-local anchors, whole-panel MSE recomputed per resample via
report323 helpers, imported read-only). Adds mean-n_pred per arm (q<0.05 excl.
own target, from stats-*.npz) as a diagnostic. Writes work/results385.json.

Usage: uv run python report385.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
import report323 as B  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402
import run385  # noqa: E402

DRAWS, SEED = 10000, 20260909
BASELINE = "cap15"
ARMS = ["cap15", "Kc1", "Kc2", "Kshuf"]
LOWER_BETTER = ("mse", "nmae")


def resample_counts(n):
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, n, (DRAWS, n))
    counts = np.zeros((DRAWS, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(DRAWS), n), draws.ravel()), 1)
    return np.vstack([np.ones(n, dtype=np.int64), counts])


def arm_raw(name, r, truth, control, own, counts, members, src):
    metrics = json.loads((src[name] / f"metrics-{name}.json").read_text())
    path = src[name] / f"stats-{name}.json".replace(".json", ".npz")
    if R.sha(path) != metrics["statistics_sha256"]:
        raise ValueError(f"statistics changed since scoring: {name}")
    with np.load(path, allow_pickle=False) as z:
        if not np.array_equal(z["targets"], r["targets"]):
            raise ValueError(f"target order drift: {name}")
        pm = bulk_lognorm_means(z["sums"].astype(float), E.TS)
        jk = z["jk"]
        qq = z["q"]
    # own-target-excluded n_pred diagnostic (scorer convention)
    keep = np.load(BASE / "work/reference.npz", allow_pickle=False)["keep"]
    genes = r["genes"]
    own_pos = np.array([np.where(genes[keep] == t)[0][0] if t in genes[keep] else -1
                        for t in r["targets"]])
    n_pred = np.array([(qq[i] < 0.05).sum() - (1 if own_pos[i] >= 0 and qq[i][own_pos[i]] < 0.05 else 0)
                       for i in range(len(r["targets"]))], dtype=float)
    mse, budget, rho, den = B.mse_bootstrap(counts, pm, truth, control, jk,
                                            r["jk"][:-1], r["jk"][-1], own)
    if abs(mse[0] - metrics["raw"]["mse"]) > 1e-10:
        raise AssertionError(f"bootstrap MSE identity failed: {name}")
    table = pl.read_csv(src[name] / f"per-target-{name}.csv")
    raw = np.empty((len(counts), len(members)))
    for j, m in enumerate(members):
        if m == "mse":
            raw[:, j] = mse
            continue
        rows = table.filter(pl.col("metric") == E.MEMBERS[m])
        mapping = dict(zip(rows["perturbation"], rows["value"]))
        values = np.array([mapping.get(str(t), np.nan) for t in r["targets"]], dtype=float)
        valid = np.isfinite(values)
        denom = counts[:, valid].sum(1)
        raw[:, j] = (counts[:, valid] @ values[valid]) / denom
        if abs(raw[0, j] - metrics["raw"][m]) > 1e-12:
            raise AssertionError(f"per-target mean identity failed: {name}/{m}")
    return metrics, raw, n_pred


def main():
    run385.verify_with_known_drift()
    ka = json.loads((WORK / "known-answer.json").read_text())
    if not ka["pass"]:
        raise ValueError("known-answer gate has not passed")
    with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
        r = {k: z[k] for k in ("targets", "genes", "sums", "jk")}
    logmeans = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = logmeans[:-1], logmeans[-1]
    lookup = {g: i for i, g in enumerate(r["genes"])}
    own = np.array([lookup[t] for t in r["targets"]])
    n = len(r["targets"])
    counts = resample_counts(n)
    members = list(E.MEMBERS)
    anchors = json.loads((BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])
    src = {"cap15": BASE / "work", "Kc1": WORK, "Kc2": WORK, "Kshuf": WORK}
    raw, overall, rows, npreds = {}, {}, [], {}
    for name in ARMS:
        metrics, raw[name], npreds[name] = arm_raw(name, r, truth, control, own,
                                                   counts, members, src)
        scaled = B.scaled(raw[name], base, replicate)
        overall[name] = scaled.mean(1)
        rows.append({"arm": name, **metrics["raw"],
                     "mean_n_pred": float(npreds[name].mean()),
                     "local_h1_overall": float(overall[name][0]),
                     "local_scaled": {m: float(scaled[0, j]) for j, m in enumerate(members)}})
    contrasts = {}
    for name in ARMS:
        if name == BASELINE:
            continue
        delta = overall[name] - overall[BASELINE]
        entry = {"overall": float(delta[0]), "overall_ci": B.ci(delta[1:]), "raw": {}}
        for j, m in enumerate(members):
            d = raw[name][:, j] - raw[BASELINE][:, j]
            e = {"point": float(d[0]), **B.ci(d[1:])}
            lo, hi = e["lo"], e["hi"]
            if m in LOWER_BETTER:
                e["beneficial"] = bool(hi < 0)
            else:
                e["beneficial"] = bool(lo > 0)
            entry["raw"][m] = e
        # H1-scaled deltas per member
        b, rep = base, replicate
        entry["scaled"] = {}
        for j, m in enumerate(members):
            d = (raw[name][:, j] - raw[BASELINE][:, j]) / abs(rep[j] - b[j])
            entry["scaled"][m] = {"point": float(d[0]), **B.ci(d[1:])}
        contrasts[name] = entry
    result = {"baseline": BASELINE, "arms": rows, "contrasts": contrasts,
              "bootstrap": {"draws": DRAWS, "seed": SEED, "unit": "target", "paired": True},
              "anchors_sha256": R.sha(BASE / "work/h1-anchors.json"),
              "scope": "85-target reused H1 development panel; not A/B/C validation"}
    R.write_json(WORK / "results385.json", result)
    for name, c in contrasts.items():
        print(f"== {name} vs cap15: overall {c['overall']:+.5f} {c['overall_ci']}")
        for m in members:
            print(f"   {m:>6}: raw {c['raw'][m]['point']:+.6f} [{c['raw'][m]['lo']:+.6f},"
                  f"{c['raw'][m]['hi']:+.6f}]  scaled {c['scaled'][m]['point']:+.5f} "
                  f"[{c['scaled'][m]['lo']:+.5f},{c['scaled'][m]['hi']:+.5f}]")


if __name__ == "__main__":
    main()
