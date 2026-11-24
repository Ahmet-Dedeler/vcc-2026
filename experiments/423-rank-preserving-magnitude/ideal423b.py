"""Unconstrained ceiling per direction at larger k (H1 base_s0 stats)."""
import numpy as np
from pathlib import Path
import run423 as RN, rp423
z = np.load(RN.ROOT / "experiments/323-current-family-validation/work/reference.npz")
genes_k = z["genes"][z["keep"]]; targets = z["targets"]; RQ = z["q"]; RL = z["lfc"]
L = np.load(RN.ROOT / "experiments/416-de-carrier-menu/work/h1/stats-base_s0.npz")["lfc"]
def nm(Lx):
    out = []
    for i in range(len(targets)):
        S = (RQ[i] < .05) & np.isfinite(RL[i]) & (RL[i] != 0) & (genes_k != targets[i])
        if S.sum() < 10: continue
        lp = np.nan_to_num(Lx[i]); out.append(np.abs(lp[S] - RL[i][S]).sum() / np.abs(RL[i][S]).sum())
    return np.mean(out)
b = nm(L)
for sel in ("dn", "up", "q80dn", "q80up"):
    row = []
    for k in (1.4, 2.0, 2.5, 3.0, 4.0):
        L2 = L.copy()
        for i in range(L.shape[0]):
            l = L[i]; ok = np.isfinite(l) & (l != 0)
            t = rp423.select(sel, l, ok); L2[i] = np.where(t, k * l, l)
        row.append(f"k{k}: {nm(L2) - b:+.4f}")
    print(sel, " ".join(row))
