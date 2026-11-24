#!/usr/bin/env python3
"""346 step 4b — two additions the registration's controls imply but do not name.

1. `mean_response_b1`: a TARGET-CONSTANT arm. Every evaluation target gets the same vector,
   the mean of the Pan H1 training-target directions, integrated at beta=1 by the same
   norm-preserving rule. It uses the recipient's labels but no per-target donor information.
   If this reproduces the same-line ridge's sign gain, the regime buys a line-level program
   offset (cheap: a handful of perturbations), not per-target direction information
   (expensive: a full non-panel bank). That distinction is the whole point for context A.
2. A paired 10,000-target bootstrap of the design-stage weighted sign accuracy deltas, so the
   design-stage gain carries an interval and not just a point.

Writes only into this experiment's work/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/334-text-conditioned-transfer"))
import core334 as P      # noqa: E402
import fit346 as F       # noqa: E402
import effects346 as X   # noqa: E402

DRAWS, SEED = 10000, 20260909
ORDER = ("cap15", "source_mean_b1", "diag_cross_b1", "same_diag_b1", "crossline_full_b1",
         "same_ridge_b05", "same_ridge_b1", "shuffled_b1", "mean_response_b1")


def constant_prediction(axes):
    """Mean unit-RMS direction of the Pan H1 training-target rows, on the covered axis."""
    rows, _ = F.bank(F.PAN, axes["training"], axes["covered"])
    direction, _ = P.rms_direction(rows)
    return direction.mean(0)


def build_extra():
    axes = F.build_axes()
    with np.load(WORK / "predictions.npz", allow_pickle=False) as z:
        covered = [str(g) for g in z["genes"]]
        test = [str(t) for t in z["test_targets"]]
    if covered != axes["covered"] or test != axes["test"]:
        raise ValueError("axis drift against the frozen fit")
    vector = constant_prediction(axes)
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        genes = [str(g) for g in z["genes"]]
        targets = [str(t) for t in z["targets"]]
        cap15 = np.asarray(z["cap15"], dtype=np.float64)
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi[g] for g in covered], dtype=np.int64)
    ti = {t: i for i, t in enumerate(targets)}
    rows = np.array([ti[t] for t in test], dtype=np.int64)
    new = cap15.copy()
    for i in rows:
        new[i, cols] = X.blend(cap15[i, cols], vector, 1.0)
    outside = np.ones(len(genes), dtype=bool)
    outside[cols] = False
    if not np.array_equal(new[:, outside], cap15[:, outside]):
        raise AssertionError("mean_response_b1 touched a gene outside the covered axis")
    if float(np.abs(np.linalg.norm(new[:, cols], axis=1)
                    - np.linalg.norm(cap15[:, cols], axis=1)).max()) > 1e-10:
        raise AssertionError("mean_response_b1 changed a covered-axis row norm")
    np.savez_compressed(WORK / "effects-extra.npz", genes=np.array(genes),
                        targets=np.array(targets), covered_columns=cols,
                        modified_rows=rows, mean_response_b1=new,
                        constant_direction=vector)
    return new


def per_target_signs(arms):
    with np.load(BASE / "work/reference.npz", allow_pickle=False, mmap_mode="r") as z:
        genes_all = np.asarray(z["genes"]).astype(str)
        keep = np.asarray(z["keep"], dtype=bool)
        targets = np.asarray(z["targets"]).astype(str)
    genes = genes_all[keep]
    gi = {g: i for i, g in enumerate(genes)}
    with np.load(WORK / "effects346.npz", allow_pickle=False) as z:
        egenes = np.asarray(z["genes"]).astype(str)
    eg = {g: i for i, g in enumerate(egenes)}
    col = np.array([eg.get(g, -1) for g in genes])
    miss = col < 0
    mapped = {}
    for a, e in arms.items():
        m = np.asarray(e, dtype=np.float64)[:, np.where(miss, 0, col)]
        m[:, miss] = 0.
        mapped[a] = m
    ref = pl.read_parquet(BASE / "work/reference_de.parquet")
    values = {a: [] for a in arms}
    eligible = []
    for ti, t in enumerate(targets):
        sub = ref.filter(pl.col("target") == t)
        g = sub["feature"].to_numpy().astype(str)
        q = sub["p_adj"].to_numpy()
        l = sub["log2_fold_change"].to_numpy().astype(float)
        gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
        if gate.sum() < 10:
            continue
        idx = np.array([gi[x] for x in g[gate]])
        l = l[gate]
        w = np.abs(l)
        D = w.sum()
        eligible.append(str(t))
        for a in arms:
            se = np.sign(mapped[a][ti, idx]) == np.sign(l)
            values[a].append(float((w * se).sum() / D))
    return eligible, {a: np.asarray(v) for a, v in values.items()}


def main() -> None:
    extra = build_extra()
    with np.load(WORK / "effects346.npz", allow_pickle=False) as z:
        arms = {a: np.asarray(z[a], dtype=np.float64) for a in ORDER if a in z.files}
    arms["mean_response_b1"] = extra
    eligible, values = per_target_signs(arms)
    rng = np.random.default_rng(SEED)
    n = len(eligible)
    draws = rng.integers(0, n, (DRAWS, n))
    counts = np.zeros((DRAWS, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(DRAWS), n), draws.ravel()), 1)
    out = {"n_eligible_targets": n, "means": {a: float(v.mean()) for a, v in values.items()},
           "paired_intervals_vs_cap15": {}}
    for a, v in values.items():
        if a == "cap15":
            continue
        d = v - values["cap15"]
        boot = (counts @ d) / n
        out["paired_intervals_vs_cap15"][a] = {
            "point": float(d.mean()),
            "lo": float(np.percentile(boot, 2.5)),
            "hi": float(np.percentile(boot, 97.5)),
            "targets_improved": int((d > 0).sum()), "targets_worsened": int((d < 0).sum())}
    d = values["same_ridge_b1"] - values["mean_response_b1"]
    boot = (counts @ d) / n
    out["same_ridge_b1_minus_mean_response_b1"] = {
        "point": float(d.mean()), "lo": float(np.percentile(boot, 2.5)),
        "hi": float(np.percentile(boot, 97.5)),
        "targets_improved": int((d > 0).sum()), "targets_worsened": int((d < 0).sum())}
    out["bootstrap"] = {"draws": DRAWS, "seed": SEED, "unit": "target", "paired": True}
    out["note"] = ("mean_response_b1 is target-constant: it uses the recipient's own training "
                   "labels but no per-target donor information. Its share of same_ridge_b1's "
                   "sign gain says how much of the regime is a line-level program offset.")
    (WORK / "signs-extra.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
