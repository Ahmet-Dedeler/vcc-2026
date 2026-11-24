#!/usr/bin/env python3
"""454b — per-target optimal blend weight vs target strength (H1, analytic): cos(unit(ours) + w unit(Pan), truth) per
covered target over a w grid, and a strength-dependent w rule w_t = w0 * (s_t / median s)^g with s = OUR key-scaled norm
(non-oracle). Reads the same inputs as mag454.py.   usage: RUN426_EMITTER=prod uv run python .../wstrength454.py"""
from __future__ import annotations
import runpy
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

HERE = Path(__file__).resolve().parent
g = runpy.run_path(str(HERE / "mag454.py"), run_name="mag454_lib")
targets, eff_bulk, cols_ours, cols_ref = g["targets"], g["eff_bulk"], g["cols_ours"], g["cols_ref"]
perts, sums, ci, logpb, X, pt, pcol = g["perts"], g["sums"], g["ci"], g["logpb"], g["X"], g["pt"], g["pcol"]
W = [0, 0.25, 0.5, 1, 2, 4, 8, 32]
rec = []
for i, t in enumerate(targets):
    k = np.flatnonzero(pt == t)
    if not len(k):
        continue
    pi = int(np.flatnonzero(perts == t)[0])
    d = (logpb(sums[pi]) - logpb(sums[ci]))[cols_ref]
    o = eff_bulk[i, cols_ours]
    pv = np.where(pcol >= 0, X[k[0]][np.maximum(pcol, 0)], np.nan)
    p = np.where(np.isfinite(pv), pv, 0.0)
    cs = []
    for w in W:
        m = o / np.linalg.norm(o) + w * p / np.linalg.norm(p)
        cs.append(float(m @ d / (np.linalg.norm(m) * np.linalg.norm(d))))
    rec.append(dict(t=t, s=float(np.linalg.norm(o)), e=float(d @ d), cs=cs))
S = np.array([r["s"] for r in rec]); E = np.array([r["e"] for r in rec]); C = np.array([r["cs"] for r in rec])
print(f"{len(rec)} covered targets; mean cos by global w:")
for j, w in enumerate(W):
    print(f"  w {w:5g}: mean cos {C[:, j].mean():+.4f}  median {np.median(C[:, j]):+.4f}")
best = np.array(W)[C.argmax(1)]
print("per-target argmax w vs our strength s: spearman", round(spearmanr(S, best)[0], 3), "| vs true energy", round(spearmanr(E, best)[0], 3))
lo, hi = S <= np.median(S), S > np.median(S)
for lab, m in (("weaker half (our key)", lo), ("stronger half", hi)):
    print(f"  {lab}: mean cos by w " + " ".join(f"{w:g}:{C[m, j].mean():+.3f}" for j, w in enumerate(W)))
print("strength rule w_t = w0 (s_t/med s)^g (mean cos; interpolated on log-w grid):")
lw = np.log(np.array(W[1:]))
for w0 in (1, 2, 4):
    for gg in (0, 0.5, 1, 2):
        wt = w0 * (S / np.median(S)) ** gg
        vals = []
        for r, w in zip(rec, wt):
            m = None
            vals.append(np.interp(np.log(np.clip(w, W[1], W[-1])), lw, np.array(r["cs"][1:])))
        print(f"  w0 {w0} g {gg}: mean cos {np.mean(vals):+.4f}")
