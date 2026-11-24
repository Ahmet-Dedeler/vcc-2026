"""280 KA-2: the pooled four-donor Rest response from the pseudobulk vs the on-disk cd4_compact Rest bank
(delta = ln(2) x DESeq2 log2FC).  Different estimators, so the preregistered pass rule is an AGREEMENT rule:
median per-target cosine >= 0.50 on targets with >= 400 pooled cells, and > 0 for >= 95 % of them.
Also the descriptive tail (donor-variable targets) is written here."""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from lib280 import DORDER, Bank, cosines, residualise  # noqa: E402

ROOT = "~/Code/virtual-cell-challenge-2026"
BK = f"{ROOT}/data/external/cd4t/cloud-compact/cd4_compact"
P = lambda *a: print(*a, flush=True)  # noqa: E731
bank = Bank(".")
SC, CC, CO = bank.scored, bank.cols_cd4, bank.cols_off
o = bank.obs
res = {}

for cond, sub in (("Rest", "cd4_rest"), ("Stim8hr", "cd4_stim8hr"), ("Stim48hr", "cd4_stim48hr")):
    perts = open(f"{BK}/{sub}/perts.txt").read().split("\n")
    genes = open(f"{BK}/{sub}/genes.txt").read().split("\n")
    delta = np.load(f"{BK}/{sub}/delta.npy", mmap_mode="r")
    assert delta.shape == (len(perts), len(genes)), (delta.shape, len(perts), len(genes))
    # bank genes -> official positions; keep those also mapped into the CD4 pseudobulk axis
    official = {g: i for i, g in enumerate([l.strip() for l in
                open(f"{ROOT}/data/challenge/gene_names.csv").read().strip().split("\n")[1:]])}
    bpos = np.asarray([official.get(g, -1) for g in genes])
    off2map = -np.ones(len(official), dtype=np.int64)
    off2map[CO] = np.arange(len(CO))
    take = (bpos >= 0) & (off2map[np.clip(bpos, 0, None)] >= 0)
    mapcols = off2map[bpos[take]]
    keep = SC[mapcols]                         # scored (non-panel) columns only
    mapcols, take_idx = mapcols[keep], np.where(take)[0][keep]

    tg = o[(o.gtype == "targeting") & o.pert.isin(perts)]
    tt = sorted(set(tg.pert))
    rows_by_t = {t: tg[tg.pert == t].i.values for t in tt}
    cellsum = {t: float(tg[tg.pert == t].n_cells.sum()) for t in tt}
    b = bank.ctrl_profile(cond, DORDER)[CC]
    M = np.zeros((len(tt), len(CC)))
    for k, t in enumerate(tt):
        M[k] = bank.profile(bank.counts(rows_by_t[t]))[CC] - b
    R = residualise(M)
    R[:, ~SC] = 0.0
    mine = R[:, mapcols]
    theirs_raw = np.asarray(delta[[perts.index(t) for t in tt]][:, take_idx], dtype=np.float64)
    theirs = theirs_raw - theirs_raw.mean(0, keepdims=True)          # same residual centring
    cs = cosines(mine, theirs)
    cs_unc = cosines(R[:, mapcols], theirs_raw)
    cells = np.array([cellsum[t] for t in tt])
    deep = cells >= 400
    md = float(np.nanmedian(cs[deep]))
    frac_pos = float(np.nanmean(cs[deep] > 0))
    scale = float(np.nanmedian(np.linalg.norm(theirs, axis=1) / np.maximum(np.linalg.norm(mine, axis=1), 1e-12)))
    P(f"KA-2 {cond}: {len(tt)} shared targets, {int(deep.sum())} with >=400 pooled cells, "
      f"{len(mapcols)} shared scored genes")
    P(f"   residual-centred cosine: median(all) {np.nanmedian(cs):+.4f}  median(>=400 cells) {md:+.4f}  "
      f"frac>0 {frac_pos:.3f}  q10 {np.nanquantile(cs[deep],.1):+.4f}  q90 {np.nanquantile(cs[deep],.9):+.4f}")
    P(f"   un-centred cosine:       median(>=400) {np.nanmedian(cs_unc[deep]):+.4f}"
      f"   |theirs|/|mine| median {scale:.2f}   max|Δ| elementwise {np.abs(mine-theirs).max():.4f}")
    P(f"   VERDICT: {'PASS' if (md >= 0.50 and frac_pos >= 0.95) else 'FAIL'} "
      f"(rule: median >= 0.50 and >0 for >= 95 %)")
    res[cond] = dict(n=len(tt), n_deep=int(deep.sum()), n_genes=int(len(mapcols)),
                     median_all=float(np.nanmedian(cs)), median_deep=md, frac_pos=frac_pos,
                     median_uncentred=float(np.nanmedian(cs_unc[deep])), norm_ratio=scale,
                     max_abs_elementwise=float(np.abs(mine - theirs).max()),
                     passed=bool(md >= 0.50 and frac_pos >= 0.95))
    np.savez(f"ka2_{cond}.npz", targets=np.asarray(tt), cos=cs, cos_uncentred=cs_unc, cells=cells)

json.dump(res, open("ka2_280.json", "w"), indent=1)
