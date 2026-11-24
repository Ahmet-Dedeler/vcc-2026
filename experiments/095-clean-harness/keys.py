"""Which computable confidence key orders the reference's DE genes best?

`reach` walks the reference's DE pool in OUR order. The pool and each gene's correctness
are fixed by the model; the order is a free parameter realised through the per-gene
emission spread. Production uses the emitted p-value order; `--keep-top` uses |effect|.
The oracle -- correct signs first -- is worth +0.103 overall, so the question is which
BUILDABLE key gets closest.

Experiment 013 nominated `sqrt(|A|.|B|) x sign-agreement` across donors on a LODO proxy
(+24% over magnitude ordering) and it was never tested against the real metric. This does
that, on the depth-repaired harness, with the two donors the prediction actually used.
"""
from __future__ import annotations
import numpy as np, os, sys

ROOT = "~/Code/virtual-cell-challenge-2026"
d = np.load(f"{ROOT}/experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
genes, perts = d["genes"].astype(str), d["perts"].astype(str)
rq, rl, pp, pq, pl = d["real_padj"], d["real_lfc"], d["pred_p"], d["pred_padj"], d["pred_lfc"]
G, ALPHA = genes.size, 0.05
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)
B, R = 0.079022, 0.964145

# ---- donor deltas on the DE-table gene axis, for the two donors the LODO run used ----
DON = {}
for name in ("k562_gwps", "orion_hek293t"):
    root = f"{ROOT}/data/donors/{name}"
    dg = np.array(open(f"{root}/genes.txt").read().split())
    dp = np.array(open(f"{root}/perts.txt").read().split())
    delta = np.load(f"{root}/delta.npy", mmap_mode="r")
    gi = {g: i for i, g in enumerate(dg)}
    cols = np.array([gi.get(g, -1) for g in genes])
    pi = {p: i for i, p in enumerate(dp)}
    DON[name] = (delta, cols, pi)
    print(f"{name}: {delta.shape}, {(cols>=0).sum()}/{G} of the DE axis covered, "
          f"{sum(1 for t in perts if t in pi)}/{perts.size} targets carried")

def donor_vec(name, target):
    delta, cols, pi = DON[name]
    if target not in pi: return None
    row = np.asarray(delta[pi[target]], dtype=np.float64)
    out = np.zeros(G); ok = cols >= 0
    out[ok] = row[cols[ok]]
    return out

def reach_of(order, idn, mt, nc):
    nd, nmt = np.cumsum(idn[order]), np.cumsum(mt[order])
    with np.errstate(invalid="ignore", divide="ignore"):
        pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
    good = np.flatnonzero((nd > 0) & (pur >= 0.9))
    return (nd[good[-1]] / nc) if good.size else 0.0

KEYS = ["production", "|lfc|", "agree_geo", "agree_only", "geo_only", "min_ab", "oracle"]
out = {k: [] for k in KEYS}
n_both = 0
for i, t in enumerate(perts):
    ok = np.ones(G, bool); j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA); ps = ok & (pq[i] < ALPHA)
    gs = np.flatnonzero(rs); nc = gs.size
    if nc == 0: continue
    rlv, plv = rl[i][gs], pl[i][gs]
    idn = (~np.isnan(rlv)) & (rlv != 0)
    mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
    absl = np.abs(np.nan_to_num(plv, nan=-np.inf)); sig = ps[gs]
    A = donor_vec("k562_gwps", t); Bv = donor_vec("orion_hek293t", t)
    have = A is not None and Bv is not None
    n_both += have
    out["production"].append(reach_of(np.lexsort((fr[gs], -absl, pp[i][gs], pq[i][gs], ~sig)), idn, mt, nc))
    out["|lfc|"].append(reach_of(np.lexsort((fr[gs], -absl)), idn, mt, nc))
    out["oracle"].append(reach_of(np.lexsort((fr[gs], ~mt)), idn, mt, nc))
    if have:
        a, b = A[gs], Bv[gs]
        geo = np.sqrt(np.abs(a) * np.abs(b))
        agree = (np.sign(a) == np.sign(b)) & (a != 0) & (b != 0)
        out["agree_geo"].append(reach_of(np.lexsort((fr[gs], -(geo * agree))), idn, mt, nc))
        out["agree_only"].append(reach_of(np.lexsort((fr[gs], -absl, ~agree)), idn, mt, nc))
        out["geo_only"].append(reach_of(np.lexsort((fr[gs], -geo)), idn, mt, nc))
        out["min_ab"].append(reach_of(np.lexsort((fr[gs], -np.minimum(np.abs(a), np.abs(b)))), idn, mt, nc))

print(f"\nperturbations carried by BOTH donors: {n_both}/{perts.size}\n")
base = np.mean(out["production"])
print(f"{'ordering key':<34}{'raw reach':>11}{'scaled':>10}{'Δoverall':>12}")
for k in KEYS:
    v = np.mean(out[k])
    note = "" if len(out[k]) == len(out["production"]) else f"  (n={len(out[k])})"
    print(f"{k:<34}{v:>11.5f}{(v-B)/(R-B):>10.4f}{(v-base)/(R-B)/6:>+12.5f}{note}")
