#!/usr/bin/env python3
"""390c -- compositional collateral of steering the CPM block's arithmetic means down (H1, analytic).
Per cell CPM sums to 1e6, so sum_g A_g == 1e6 for any placement: whatever the steered block loses in arithmetic-mean
CPM, the unsteered genes gain. Head = top-400 |design| in the universe, block = top-2200 ref_lin among non-head (361).
A0 = delivered arithmetic mean under welding = prediction bulk CPM. Targets: every block gene at ref_lin * 2^-eps.
Reports the implied uniform log2 shift on all other genes and how many head/other calls it would flip."""
import numpy as np, json
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
keep = r["keep"].astype(bool); genes = r["genes"].astype(str)[keep]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); sig = np.isfinite(s["q"]) & (s["q"] < .05)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
ref_lin = r["ref_lin"].astype(float)
sums = s["sums"].astype(float); A0full = sums / sums.sum(1, keepdims=True) * 1e6; A0 = A0full[:, keep]
nT = len(targets); rows = []
for eps in (0.002, 0.01, 0.05):
    shifts, bl_mean, flips_head, flips_other, frac_up = [], [], [], [], []
    for i in range(nT):
        head = np.zeros(len(genes), bool); head[np.argsort(-np.abs(E1[i]), kind="stable")[:400]] = True
        cpm = np.where(head, -np.inf, ref_lin); block = np.zeros(len(genes), bool); block[np.argsort(-cpm, kind="stable")[:2200]] = True
        SA0 = A0[i, block].sum(); SR = ref_lin[block].sum()
        R = (ref_lin[block] * 2 ** -eps).sum() - SA0            # change in the block's arithmetic mass
        O = 1e6 - SA0                                            # everything else carries -R
        sh = np.log2((O - R) / O); shifts.append(sh); bl_mean.append(np.log2(SA0 / SR))
        frac_up.append(float((E0[i, block] > 0).mean()))
        oth = sig[i] & ~block
        new = E0[i] + sh
        flips = (np.sign(new) != np.sign(E0[i])) & oth
        flips_head.append(float((flips & head).sum() / max((oth & head).sum(), 1)))
        flips_other.append(float((flips & ~head).sum() / max((oth & ~head).sum(), 1)))
    q = lambda v: [round(float(x), 5) for x in np.quantile(v, [.1, .5, .9])]
    row = {"eps": eps, "other_gene_log2_shift_q10_50_90": q(shifts), "block_massweighted_delivered_lfc_q10_50_90": q(bl_mean),
           "block_frac_delivered_up_median": round(float(np.median(frac_up)), 4),
           "head_calls_flipped_q10_50_90": q(flips_head), "other_calls_flipped_q10_50_90": q(flips_other)}
    rows.append(row); print(json.dumps(row))
(OUT / "collateral390.json").write_text(json.dumps(rows, indent=1))
