"""Analytic ceiling on the dumped targets: nmae if every selected gene hit k*lfc exactly (no rank constraint)."""
import numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
# use the full base_s0 stats from 416 for all 85 targets (identical construction)
st = np.load(RN.ROOT / "experiments/416-de-carrier-menu/work/h1/stats-base_s0.npz")
L = st["lfc"]
def nmae(lfc_all):
    out = []
    for i in range(len(r["targets"])):
        ql = r["q"][i]; ll = r["lfc"][i]
        S = (ql < .05) & np.isfinite(ll) & (ll != 0) & (genes_k != r["targets"][i])
        if S.sum() < 10: continue
        lp = np.where(np.isfinite(lfc_all[i]), lfc_all[i], 0.0)
        out.append(np.abs(lp[S] - ll[S]).sum() / np.abs(ll[S]).sum())
    return np.mean(out), len(out)
print("base", nmae(L))
for sel in ("all", "q50", "q80", "up", "dn"):
    row = []
    for k in (0.8, 1.2, 1.3, 1.4, 1.6, 2.0):
        L2 = L.copy()
        for i in range(L.shape[0]):
            l = L[i]; ok = np.isfinite(l) & (l != 0)
            t = rp423.select(sel, l, ok)
            L2[i] = np.where(t, k * l, l)
        row.append(f"k{k}: {nmae(L2)[0] - nmae(L)[0]:+.4f}")
    print(sel, " ".join(row))
