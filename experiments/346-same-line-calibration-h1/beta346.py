#!/usr/bin/env python3
"""346 step 9 — the blend-weight trade-off curve, at the design stage only.

Step 8 established that design-stage cosine rank retrieval tracks the delivered `pds` member
across the five emitted arms at r = 0.999 and slope 1.01. This sweeps the blend weight beta
and reports both design-stage endpoints — |lfc|-weighted sign accuracy (which 343 ties to
nmae/fid/reach) and retrieval (which step 8 ties to pds) — so the shape of the trade-off is
visible without emitting nine more arms. These are proxies, not scored members: only the five
emitted arms carry six-member numbers. Writes only into work/.
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
for _p in (str(BASE), str(ROOT / "src"), str(ROOT / "experiments/334-text-conditioned-transfer")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import run323 as R           # noqa: E402
import evaluate323 as E      # noqa: E402
import core334 as P          # noqa: E402
import effects346 as X       # noqa: E402
from cell_eval2.prep import bulk_lognorm_means   # noqa: E402

BETAS = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)


def main() -> None:
    with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
        targets = np.asarray(z["targets"]).astype(str)
        genes = np.asarray(z["genes"]).astype(str)
        keep = np.asarray(z["keep"], dtype=bool)
        sums = z["sums"]
    logmeans = bulk_lognorm_means(sums.astype(float), E.TS)
    truth = logmeans[:-1] - logmeans[-1]
    ranked = np.array([g not in set(targets) for g in genes])
    dgenes = genes[keep]
    dgi = {g: i for i, g in enumerate(dgenes)}

    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        egenes = [str(g) for g in z["genes"]]
        cap15 = np.asarray(z["cap15"], dtype=np.float64)
    with np.load(WORK / "predictions.npz", allow_pickle=False) as z:
        covered = [str(g) for g in z["genes"]]
        test = [str(t) for t in z["test_targets"]]
        prediction = np.asarray(z["same_line_full"], dtype=np.float64)
    gi = {g: i for i, g in enumerate(egenes)}
    cols = np.array([gi[g] for g in covered], dtype=np.int64)
    ti = {t: i for i, t in enumerate(targets)}
    rows = np.array([ti[t] for t in test], dtype=np.int64)
    col = np.array([gi.get(g, -1) for g in dgenes])
    miss = col < 0

    ref = pl.read_parquet(BASE / "work/reference_de.parquet")
    gates = []
    for tix, t in enumerate(targets):
        sub = ref.filter(pl.col("target") == t)
        g = sub["feature"].to_numpy().astype(str)
        q = sub["p_adj"].to_numpy()
        l = sub["log2_fold_change"].to_numpy().astype(float)
        gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
        if gate.sum() < 10:
            continue
        gates.append((tix, np.array([dgi[x] for x in g[gate]]), l[gate]))

    curve = []
    for beta in BETAS:
        effect = cap15.copy()
        for j, i in enumerate(rows):
            effect[i, cols] = X.blend(cap15[i, cols], prediction[j], beta)
        mapped = R.map_axis(effect, np.array(egenes), genes)
        retrieval = float(P.endpoints(mapped[:, ranked], truth[:, ranked])["retrieval"].mean())
        design = effect[:, np.where(miss, 0, col)]
        design[:, miss] = 0.
        signs = []
        for tix, idx, l in gates:
            w = np.abs(l)
            se = np.sign(design[tix, idx]) == np.sign(l)
            signs.append(float((w * se).sum() / w.sum()))
        curve.append({"beta": beta, "design_weighted_sign": float(np.mean(signs)),
                      "design_retrieval": retrieval})
        print(json.dumps(curve[-1]), flush=True)

    base = curve[0]
    for row in curve:
        row["sign_gain"] = row["design_weighted_sign"] - base["design_weighted_sign"]
        row["retrieval_loss"] = base["design_retrieval"] - row["design_retrieval"]
    out = {"curve": curve,
           "beta_zero_reproduces_cap15": {
               "design_weighted_sign": base["design_weighted_sign"],
               "expected_343": 0.600185097698462,
               "absolute_difference": abs(base["design_weighted_sign"] - 0.600185097698462)},
           "note": ("proxies, not scored members. Step 8 measured design_retrieval -> delivered "
                    "pds at slope 1.01 (r=0.999) across the five emitted arms; 343 ties "
                    "design_weighted_sign to nmae/fid/reach. Only the five emitted arms in "
                    "results.json carry actual six-member numbers.")}
    (WORK / "beta-curve.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out["beta_zero_reproduces_cap15"], indent=2))


if __name__ == "__main__":
    main()
