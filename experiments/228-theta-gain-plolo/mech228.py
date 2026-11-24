"""228 — KA-3 and the mechanism check: is the gain fitted HERE 189's `theta`, and is it the
same object across folds?

224 measured `r(gain, 189 theta) = +0.81 / +0.90` and `r` between gains fitted for different
held lines `= +0.95`, on its own surface.  If this experiment's gain -- fitted on the official
axis, on production's own consensus, through the fold's own operator -- is a different object,
the two experiments are not measuring the same thing and 228's number does not price 224's.

    uv run python mech228.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import FOLDS, OFFICIAL, build_fold  # noqa: E402


def spearman(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


def main() -> None:
    keys = ("gain", "gain_reg")
    Z = {f: np.load(WORK / f"gain-{f}.npz", allow_pickle=True) for f in FOLDS}
    out: dict = {"folds": list(FOLDS)}

    print("1. STABILITY -- is the gain a property of the GENE or of the fold?")
    for key in keys:
        print(f"   [{key}]")
        for a in FOLDS:
            for b in FOLDS:
                if a >= b:
                    continue
                m = np.asarray(Z[a]["live"], bool) & np.asarray(Z[b]["live"], bool)
                r = float(np.corrcoef(Z[a][key][m], Z[b][key][m])[0, 1])
                s = spearman(Z[a][key][m], Z[b][key][m])
                out[f"stability_{key}_{a}_{b}"] = r
                print(f"     r({a:14s}, {b:14s}) = {r:+.4f}   spearman {s:+.4f}   "
                      f"n = {int(m.sum())}")

    print("\n2. KA-3 -- is it 189's `theta`?")
    p189 = ROOT / "experiments/189-per-gene-reliability/work/prior189.npz"
    if p189.exists():
        z = np.load(p189, allow_pickle=True)
        idx = {str(g): i for i, g in enumerate(z["genes"])}
        pos = [(i, idx[g]) for i, g in enumerate(OFFICIAL) if g in idx]
        rows = np.asarray([i for i, _ in pos]); cols = np.asarray([j for _, j in pos])
        for key in keys:
            gbar = np.mean([np.asarray(Z[f][key]) for f in FOLDS], axis=0)
            livebar = np.all([np.asarray(Z[f]["live"], bool) for f in FOLDS], axis=0)
            fullbar = np.all([np.asarray(Z[f]["full_cov"], bool) for f in FOLDS], axis=0)
            for restrict, mask in (("all-live", livebar), ("FULL-COVERAGE", livebar & fullbar)):
                for stat in ("theta", "eb", "raw"):
                    a = gbar[rows]; b = np.asarray(z[stat])[cols]
                    ok = np.isfinite(a) & np.isfinite(b) & mask[rows]
                    r = float(np.corrcoef(a[ok], b[ok])[0, 1])
                    out[f"corr_189_{stat}_{key}_{restrict}"] = r
                    print(f"   [{restrict:13s}] r({key:9s}, 189 {stat:6s}) = {r:+.4f}   "
                          f"spearman {spearman(a[ok], b[ok]):+.4f}   n = {int(ok.sum())}")
    else:
        print("   189's prior189.npz not on disk")

    print("\n3. What predicts the gain here?")
    fold = build_fold("orion_hct116")
    b_ctrl = np.asarray(fold["b_ctrl"], dtype=np.float64)
    for key in keys:
        gbar = np.mean([np.asarray(Z[f][key]) for f in FOLDS], axis=0)
        livebar = np.all([np.asarray(Z[f]["live"], bool) for f in FOLDS], axis=0)
        cc = np.mean([np.asarray(Z[f]["cov_count"]) for f in FOLDS], axis=0)
        for nm, v in (("basal b_ctrl", b_ctrl), ("fit-pool coverage", cc),
                      ("log den", np.log10(np.maximum(
                          np.mean([np.asarray(Z[f]["den"]) for f in FOLDS], axis=0), 1e-300)))):
            ok = livebar & np.isfinite(v)
            r = float(np.corrcoef(gbar[ok], v[ok])[0, 1])
            out[f"pred_{key}_{nm}"] = r
            print(f"   r({key:9s}, {nm:20s}) = {r:+.4f}")

    print("\n4. Distribution")
    for key in keys:
        for f in FOLDS:
            g = np.asarray(Z[f][key]); m = np.asarray(Z[f]["live"], bool)
            q = np.quantile(g[m], [0.1, 0.5, 0.9])
            print(f"   {key:9s} {f:14s} p10 {q[0]:+.4f}  p50 {q[1]:+.4f}  p90 {q[2]:+.4f}"
                  f"   frac<0 {float((g[m] < 0).mean()):.3f}")
            out[f"dist_{key}_{f}"] = [float(x) for x in q]

    (WORK / "mech228.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'mech228.json'}")


if __name__ == "__main__":
    main()
