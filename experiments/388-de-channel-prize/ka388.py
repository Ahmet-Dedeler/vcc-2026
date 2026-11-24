#!/usr/bin/env python3
"""388 KA — reproduce 343's numbers + 348 O1 with own code, read-only on 323/work.

Conventions (from decompose343.py / signs343.py / reach348.py, reimplemented here):
  gate: reference q<0.05, finite nonzero lfc, own target removed, >=10 genes/target.
  weights: |l_ref|. KA targets: nmae 0.950836, pw delivered 0.588, designed 0.600,
  cf all-signs-right 0.803, best-scalar 0.928, both 0.658, O1 reach 0.25434726.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import json
import numpy as np
import polars as pl
from pathlib import Path

REPO = Path("~/Code/virtual-cell-challenge-2026")
W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/388-de-channel-prize/work"
OUT.mkdir(parents=True, exist_ok=True)

s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz")
e = np.load(W / "effects.npz")
genes_all = r["genes"].astype(str)
keep = r["keep"].astype(bool)
genes = genes_all[keep]
targets = s["targets"].astype(str)
assert (targets == r["targets"].astype(str)).all()
assert (e["targets"].astype(str) == targets).all()
ref = pl.read_parquet(W / "reference_de.parquet")

# reference.npz q/lfc vs parquet agreement (reach348 build_long asserts exact)
qrm = np.asarray(r["q"], dtype=float)
lrm = np.asarray(r["lfc"], dtype=float)
gi = {g: i for i, g in enumerate(genes)}
ti = {t: i for i, t in enumerate(targets)}
tcol = ref["target"].to_numpy().astype(str)
gcol = ref["feature"].to_numpy().astype(str)
ii = np.array([ti[t] for t in tcol])
jj = np.array([gi[g] for g in gcol])
qq = np.empty((len(targets), len(genes)))
ll = np.empty((len(targets), len(genes)))
qq[ii, jj] = ref["p_adj"].to_numpy().astype(float)
ll[ii, jj] = ref["log2_fold_change"].to_numpy().astype(float)
print(f"parquet-vs-npz: max|dq|={np.max(np.abs(qq-qrm)):.3e} "
      f"max|dlfc|={np.nanmax(np.abs(ll-lrm)):.3e}")

lhat_all = np.asarray(s["lfc"], dtype=float)
qhat_all = np.asarray(s["q"], dtype=float)
phat_all = np.asarray(s["p"], dtype=float)

# E1: designed cap15 mapped by gene name onto universe (signs343 convention)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes])
miss = col < 0
print("DE-universe genes absent from effects axis:", int(miss.sum()))
eff = np.asarray(e["cap15"], dtype=float)[:, np.where(miss, 0, col)]
eff[:, miss] = 0.0

grid = np.concatenate([np.linspace(0.25, 1, 4), np.linspace(1.25, 6, 20)])
rows = []
pw_des = []
for k, t in enumerate(targets):
    m = (tcol == t)
    g = gcol[m]
    q = ref["p_adj"].to_numpy().astype(float)[m]
    l = ref["log2_fold_change"].to_numpy().astype(float)[m]
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    idx = np.array([gi[x] for x in g[gate]])
    l = l[gate]
    lh = np.where(np.isfinite(lhat_all[k, idx]), lhat_all[k, idx], 0.)
    ef = eff[k, idx]
    D = np.abs(l).sum()
    nmae = np.abs(lh - l).sum() / D
    pw = (np.abs(l) * (np.sign(lh) == np.sign(l))).sum() / D
    pw_des.append((np.abs(l) * (np.sign(ef) == np.sign(l))).sum() / D)
    cf_sign = np.abs(np.abs(lh) - np.abs(l)).sum() / D
    vals = [np.abs(a * lh - l).sum() / D for a in grid]
    j = int(np.argmin(vals))
    cf_scale = vals[j]
    # "both" (343's exact cf_signs_and_scalar): LS scalar + correct signs
    a_ls = (np.abs(lh) * np.abs(l)).sum() / max((lh ** 2).sum(), 1e-12)
    cf_both = np.abs(a_ls * np.abs(lh) - np.abs(l)).sum() / D
    # extra: scalar re-optimized GIVEN correct signs (grid search on |lh| vs |l|)
    cf_both_reopt = min(np.abs(a * np.abs(lh) - np.abs(l)).sum() / D for a in grid)
    rows.append(dict(nmae=nmae, pw=pw, cf_sign=cf_sign, cf_scale=cf_scale,
                     cf_both=cf_both, cf_both_reopt=cf_both_reopt))

nmae = float(np.mean([x["nmae"] for x in rows]))
pw_m = float(np.mean([x["pw"] for x in rows]))
pw_d = float(np.mean(pw_des))
cf_s = float(np.mean([x["cf_sign"] for x in rows]))
cf_sc = float(np.mean([x["cf_scale"] for x in rows]))
cf_b = float(np.mean([x["cf_both"] for x in rows]))
cf_br = float(np.mean([x["cf_both_reopt"] for x in rows]))
print(f"n_targets={len(rows)}")
print(f"nmae        {nmae:.6f}  (KA 0.950836)")
print(f"pw delivered{pw_m:.6f}  (KA 0.588)")
print(f"pw designed {pw_d:.6f}  (KA 0.600)")
print(f"cf signs    {cf_s:.6f}  (KA 0.803)")
print(f"cf scalar   {cf_sc:.6f}  (KA 0.928)")
print(f"cf both     {cf_b:.6f}  (KA 0.658)")
print(f"cf both reopt {cf_br:.6f}  (extra, grid re-optimized)")

# reach O1: verbatim walk_reach + current order key (reach348.py:125-145,178-187)
def walk_reach(order_idx, match, adjud, n_conf):
    """order_idx: per-target ascending rank positions over pool genes."""
    n = match.shape[0]
    R = np.full(n, np.nan)
    Ks = np.zeros(n, dtype=int)
    for i in range(n):
        if n_conf[i] == 0:
            continue
        m = match[i][order_idx[i]]
        a = adjud[i][order_idx[i]]
        ma = m[a]
        denom = np.arange(1, ma.size + 1)
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / denom
        ok = np.flatnonzero(pur >= 0.9)
        k = int(ok.max() + 1) if ok.size else 0
        Ks[i] = k
        R[i] = k / n_conf[i]
    return R, Ks

n, gg = len(targets), len(genes)
ref_sig = qrm < 0.05
own = np.equal.outer(targets, genes)
pool = ref_sig & ~own
n_conf = pool.sum(1)
lr = lrm
adjud = np.isfinite(lr) & (lr != 0)
lp = lhat_all
commit = np.isfinite(lp) & (lp != 0)
match = adjud & commit & (np.sign(lp) == np.sign(lr))
feat = np.tile(genes, (n, 1))
qm = np.where(np.isfinite(qhat_all), qhat_all, np.inf)
pm = np.where(np.isfinite(phat_all), phat_all, np.inf)
ala = np.abs(lp)
ala_m = np.where(np.isfinite(ala), -ala, np.inf)
sig = np.isfinite(qhat_all) & (qhat_all < 0.05)
cur = np.empty((n, gg), dtype=int)
for i in range(n):
    cur[i] = np.lexsort((feat[i], ala_m[i], pm[i], qm[i], ~sig[i]))
Rmine, Ks = walk_reach(cur, match & pool, adjud & pool, n_conf)
print(f"O1 walk mean {float(np.nanmean(Rmine)):.8f}  (KA 0.25434726)")
print(f"pool-empty targets {(n_conf==0).sum()}")

ka = dict(n_targets=len(rows), nmae=nmae, pw_delivered=pw_m, pw_designed=pw_d,
          cf_sign=cf_s, cf_scalar=cf_sc, cf_both=cf_b, cf_both_reopt=cf_br,
          o1_reach=float(np.nanmean(Rmine)))
(OUT / "ka388.json").write_text(json.dumps(ka, indent=1))
checks = dict(nmae=abs(nmae - 0.950836) <= 1e-4,
              pw_d=abs(pw_m - 0.588) <= 1e-3 + 1e-4,
              pw_des=abs(pw_d - 0.600) <= 1e-3 + 1e-4,
              cf_sign=abs(cf_s - 0.803) <= 1e-3,
              cf_scale=abs(cf_sc - 0.928) <= 1e-3,
              cf_both=abs(cf_b - 0.658) <= 1e-3,
              o1=abs(float(np.nanmean(Rmine)) - 0.25434726) <= 1e-4)
print(json.dumps(checks, indent=1))
