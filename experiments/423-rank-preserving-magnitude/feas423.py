#!/usr/bin/env python3
"""Feasible fraction of the target LFC change (achieved dA / target dA, clipped to [-1, 2]) per touched gene, from
work/h1/genes-<variant>_s<seed>.npy (cols: target, gene, frac, lfc, dA). Split by direction and by |lfc| quantile
within each target, both gene-count-weighted and |dA|-weighted.
usage: uv run python feas423.py VARIANT[,VARIANT...] [SEEDTAG=s0]"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
tag = sys.argv[2] if len(sys.argv) > 2 else "s0"
for v in sys.argv[1].split(","):
    f = HERE / f"work/h1/genes-{v}_{tag}.npy"
    if not f.exists():
        print(v, "missing"); continue
    G = np.load(f)
    t, fr, l, dA = G[:, 0], G[:, 2], G[:, 3], G[:, 4]
    # |lfc| quantile within target
    qn = np.zeros(len(G))
    for ti in np.unique(t):
        m = t == ti
        a = np.abs(l[m])
        qn[m] = (np.argsort(np.argsort(a)) + 0.5) / a.size
    w = np.abs(dA)
    def row(m, name):
        if not m.any():
            return f"{name}: -"
        return (f"{name}: n/target {m.sum() / len(np.unique(t)):.0f}, any>0 {np.mean(fr[m] > 0.005):.3f}, "
                f"mean {fr[m].mean():.3f}, |dA|-w {np.sum(fr[m] * w[m]) / w[m].sum():.3f}, full(>=.95) {np.mean(fr[m] >= .95):.3f}")
    print(f"== {v} ({tag})")
    print("  " + row(np.ones(len(G), bool), "all"))
    print("  " + row(l > 0, "up"))
    print("  " + row(l < 0, "down"))
    for lo, hi in ((0, .5), (.5, .8), (.8, .95), (.95, 1.01)):
        print("  " + row((qn >= lo) & (qn < hi), f"|lfc| q{lo:.2f}-{min(hi, 1):.2f}"))
