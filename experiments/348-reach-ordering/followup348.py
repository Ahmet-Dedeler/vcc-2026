#!/usr/bin/env python3
"""348 follow-up: paired-target bootstrap CIs for ordering deltas; per-target
O4-O1 anatomy; two post-hoc combos (labelled exploratory)."""
from __future__ import annotations
import json
import numpy as np
import polars as pl
from pathlib import Path
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

df = pl.read_csv(WORK / "per-target-orderings.csv").sort("target")
targets = df["target"].to_list()
ordr = json.load(open(WORK / "orderings.json"))
schemes = [c for c in df.columns if c != "target"]
Y = {c: df[c].to_numpy() for c in schemes}  # NaN where pool empty
scored = ~np.isnan(Y["O1_current"])
print(f"scored targets: {scored.sum()}", flush=True)
names = [t for t, s in zip(targets, scored) if s]

rng = np.random.default_rng(34809)
B = 10000
idx = np.flatnonzero(scored)
base = Y["O1_current"][idx]
out = {}
for c in schemes:
    if c == "O1_current":
        continue
    d = Y[c][idx] - base
    boots = rng.choice(len(idx), size=(B, len(idx)), replace=True)
    means = d[boots].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    out[c] = {"delta": float(d.mean()), "ci95": [float(lo), float(hi)],
              "frac_targets_pos": float(np.mean(d > 0)),
              "frac_targets_neg": float(np.mean(d < 0))}
    print(f"{c}: delta={d.mean():+.6f} 95%CI=[{lo:+.6f},{hi:+.6f}] "
          f"pos={out[c]['frac_targets_pos']:.2f} neg={out[c]['frac_targets_neg']:.2f}",
          flush=True)

# O4-O1 anatomy: pool size vs delta
ROOT = Path(__file__).resolve().parents[2]
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz")
rq = ref["q"]
rt = ref["targets"].astype(str)
rg = ref["genes"].astype(str)[ref["keep"]]
pos = {t: i for i, t in enumerate(rt)}
nconf = {}
for t in names:
    i = pos[t]
    nconf[t] = int(((rq[i] < 0.05) & (rg != t)).sum())
d4 = Y["O4_abs_design"][idx] - base
for t, dd in sorted(zip(names, d4), key=lambda x: -x[1])[:8]:
    print(f"  top O4 gain: {t} n_conf={nconf[t]} d={dd:+.5f}", flush=True)
for t, dd in sorted(zip(names, d4), key=lambda x: x[1])[:5]:
    print(f"  top O4 loss: {t} n_conf={nconf[t]} d={dd:+.5f}", flush=True)
big = np.array([nconf[t] for t in names])
print(f"corr(n_conf, O4delta)={np.corrcoef(big, d4)[0,1]:+.3f} "
      f"mean d | n_conf<=20: {d4[big<=20].mean():+.5f} (n={(big<=20).sum()}) "
      f"n_conf>20: {d4[big>20].mean():+.5f} (n={(big>20).sum()})", flush=True)
json.dump(out, open(WORK / "delta-ci.json", "w"), indent=2)
