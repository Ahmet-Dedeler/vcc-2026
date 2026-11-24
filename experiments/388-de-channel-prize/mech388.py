#!/usr/bin/env python3
"""388 mechanism gate — cell_eval2's own functions on synthetic blocks.

Two 400-cell x 2000-gene integer count blocks, IDENTICAL per-gene totals, different
per-cell library layouts:
  A (welded): every cell total = 20_000.
  B (split):  three library classes 14k/20k/26k x133/134/133 (grand total 8.00M,
              brief's approx 14k/20k/27k with the total pinned so bulks can match).
Background genes: counts proportional to library size (layout-neutral CPM).
Probe gene g* (index 0): fixed total S=4000; fraction f in the low class, rest even.
Reference group: 400 welded cells, same per-gene totals (no DE except via layout).

Shows: prep.bulk_lognorm_means bit-identical; compute_lfc_table (arithmetic,
target_sum 1e6, eps 1e-9) log2FC differs. Sweep f -> |dlog2FC| curve.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import json
import numpy as np
import anndata as ad
import pandas as pd
from pathlib import Path

from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.de_compute import _to_linear, compute_lfc_table

OUT = Path("~/Code/virtual-cell-challenge-2026/"
           "experiments/388-de-channel-prize/work")

rng = np.random.default_rng(388)
N, G = 400, 2000
LIBS = np.array([14_000, 20_000, 26_000])
CLS = np.array([133, 134, 133])
WELD = 20_000
S_PROBE = 4000
GPROBE = 0


def lib_of(layout):
    return np.full(N, WELD) if layout == "weld" else np.repeat(LIBS, CLS)


def spread(total, weights):
    """Integer vector summing to `total`, proportional to `weights` (largest remainder)."""
    w = np.asarray(weights, dtype=float)
    x = np.floor(total * w / w.sum()).astype(np.int64)
    d = int(total - x.sum())
    if d > 0:
        frac = total * w / w.sum() - x
        order = np.argsort(-frac + rng.random(w.size) * 1e-9, kind="stable")
        x[order[:d]] += 1
    return x


# background totals: sum to grand - probe (same grand 8M for both layouts)
GRAND = N * WELD
assert int(np.repeat(LIBS, CLS).sum()) == GRAND
raw = rng.integers(200, 3000, size=G - 1).astype(float)
base = spread(GRAND - S_PROBE, raw)  # exact-sum background totals


def build_block(layout, f_low):
    lib = lib_of(layout)
    X = np.zeros((N, G), dtype=np.int64)
    for j in range(1, G):
        X[:, j] = spread(int(base[j - 1]), lib)  # layout-neutral: proportional to L_i
    if layout == "weld":
        X[:, GPROBE] = spread(S_PROBE, np.ones(N))
    else:
        n_low = CLS[0]
        s_low = int(round(f_low * S_PROBE))
        X[:n_low, GPROBE] = spread(s_low, np.ones(n_low))
        rest = np.arange(n_low, N)
        X[rest, GPROBE] = spread(S_PROBE - s_low, np.ones(N - n_low))
    assert (X[:, GPROBE].sum() == S_PROBE)
    assert (X[:, 1:].sum(axis=0) == base).all()
    # libraries approximate (rounding noise << class separation); classes must stay distinct
    dev = np.abs(X.sum(axis=1) - lib).max()
    assert dev <= 150, dev
    return X


def lfc_of(Xt, Xref):
    X = np.vstack([Xt, Xref]).astype(np.float64)
    obs = pd.DataFrame({"grp": ["pert"] * N + ["ctrl"] * N})
    a = ad.AnnData(X=X, obs=obs)
    a.var_names = [f"g{j}" for j in range(G)]
    lin = _to_linear(a, "counts", target_sum=1e6)
    tab = compute_lfc_table(lin, groupby="grp", reference="ctrl",
                            mean_calc="arithmetic", epsilon=1e-9)
    return tab["log2_fold_change"].to_numpy()


Xref = build_block("weld", 0.0)  # reference: welded, probe uniform
results = {}
lfc_bg = {}
for f_low in [0.0, 1 / 3, 0.5, 0.75, 1.0]:
    Xa = build_block("weld", f_low)
    Xb = build_block("split", f_low)
    assert (Xa.sum(axis=0) == Xb.sum(axis=0)).all(), "per-gene totals identical"
    sa = np.vstack([Xa.sum(axis=0, keepdims=True), Xref.sum(axis=0, keepdims=True)]).astype(float)
    sb = np.vstack([Xb.sum(axis=0, keepdims=True), Xref.sum(axis=0, keepdims=True)]).astype(float)
    ba = bulk_lognorm_means(sa, 5e4)
    bb = bulk_lognorm_means(sb, 5e4)
    la, lb = lfc_of(Xa, Xref), lfc_of(Xb, Xref)
    results[str(f_low)] = dict(bulk_maxdiff=float(np.max(np.abs(ba - bb))),
                               lfc_weld=float(la[GPROBE]),
                               lfc_split=float(lb[GPROBE]),
                               dlfc=float(lb[GPROBE] - la[GPROBE]))
    if f_low == 0.5:
        d = np.abs(lb - la)
        lfc_bg = dict(max=float(d[1:].max()), med=float(np.median(d[1:])))

for f, v in results.items():
    print(f"f_low={f:>5}: bulk_max|d|={v['bulk_maxdiff']:.3e} "
          f"LFC_weld={v['lfc_weld']:+.4f} LFC_split={v['lfc_split']:+.4f} "
          f"Delta={v['dlfc']:+.4f}")
print(f"background genes (f=0.5): max|Delta|={lfc_bg['max']:.4f} median={lfc_bg['med']:.2e}")
(OUT / "mechanism388.json").write_text(json.dumps(
    {"sweep": results, "background": lfc_bg,
     "libs": [14000, 20000, 26000], "cells": [133, 134, 133],
     "probe_total": S_PROBE, "note": "bulk fns: prep.bulk_lognorm_means(prep.py:329); "
     "LFC fns: de_compute._to_linear(:194)+compute_lfc_table(:334), arithmetic, eps 1e-9"},
    indent=1))
print("wrote work/mechanism388.json")
