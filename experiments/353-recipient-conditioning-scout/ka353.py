"""353 known-answer check: independently recompute 343's delivered |lfc|-weighted
sign accuracy (recorded 0.588) and mean NMAE (recorded 0.9508; 323's
metrics-cap15.json raw nmae 0.9508364630325662) from 323's saved arrays.
Independent code path (pandas, own loop); reference.npz opened mmap'd so only
the small genes/keep arrays are touched."""
import numpy as np, pandas as pd
from pathlib import Path

W = Path("~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work")
s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz", mmap_mode="r")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
lhat = np.asarray(s["lfc"], dtype=np.float64)
ref = pd.read_parquet(W / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}

wacc_num = wacc_den = 0.0
nmaes, n_gate = [], []
for ti, t in enumerate(targets):
    sub = ref[ref["target"] == t]
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy()
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < 0.05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    idx = np.array([gi[x] for x in g[gate]])
    l = l[gate]
    w = np.abs(l)
    lh = np.nan_to_num(lhat[ti, idx])
    wacc_num += (w * (np.sign(lh) == np.sign(l))).sum()
    wacc_den += w.sum()
    nmaes.append(np.abs(lh - l).sum() / w.sum())
    n_gate.append(int(gate.sum()))

print("n_targets:", len(nmaes), "mean_gate:", float(np.mean(n_gate)))
print("weighted_sign_all:", wacc_num / wacc_den, "(recorded 0.588)")
print("mean_nmae:", float(np.mean(nmaes)), "(recorded 0.9508364630325662)")
