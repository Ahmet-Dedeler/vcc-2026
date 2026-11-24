"""KNOWN-ANSWER ROW for the DE pipeline: strong K562 knockdowns (erythroid TFs, Pol II, ribosome) must give
large Wilcoxon sets under the same recipe, else the small panel-target sets are a pipeline bug, not biology.
Same 18,400 control rows as extract_and_de.py (extract_rows.npz), same tie-corrected recipe, cap 400 cells,
plus a 150-cell subsample of each (the panel's typical count)."""
import time, h5py, numpy as np, pandas as pd, resource
from pathlib import Path
from scipy.stats import rankdata, norm
P = "~/Code/virtual-cell-challenge-2026/data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
W = Path("~/Code/virtual-cell-challenge-2026/experiments/260-raw-de-diagnostic-feasibility/work")
STRONG = ["GATA1", "TAL1", "LMO2", "POLR2A", "RPL3", "RPS3", "HSPA5", "MYC", "SUPT5H", "BRD4", "KLF1", "NFE2"]
CAP, BLOCK, ALPHA = 400, 3886, 0.05
rng = np.random.default_rng(0)
ctrl_rows = np.load(W / "extract_rows.npz")["ctrl_rows"]
with h5py.File(P, "r") as f:
    cats = f["obs/perturbation/categories"][:].astype(str); codes = f["obs/perturbation/codes"][:]
    genes = f["var/gene_name"][:].astype(str); idx = {c: i for i, c in enumerate(cats)}
    keep = {}
    for t in STRONG:
        if t in idx:
            r = np.flatnonzero(codes == idx[t])
            if len(r) >= 150: keep[t] = np.sort(rng.choice(r, min(CAP, len(r)), replace=False))
    print("strong targets present with >=150 cells:", {t: len(r) for t, r in keep.items()}, flush=True)
    sel = np.concatenate([ctrl_rows] + [keep[t] for t in keep]); order = np.argsort(sel); ss = sel[order]
    X = f["X"]; n = X.shape[0]; M = np.zeros((len(sel), X.shape[1]), np.float32); t0 = time.time()
    for s in range(0, n, BLOCK):
        e = min(s + BLOCK, n); lo, hi = np.searchsorted(ss, s), np.searchsorted(ss, e)
        if hi > lo:
            blk = X[s:e, :]; M[order[lo:hi]] = blk[ss[lo:hi] - s]
    print(f"pass {time.time() - t0:.0f}s", flush=True)
tot = M.sum(1, keepdims=True); N_CTRL = len(ctrl_rows)
ctrl_cpm = (M[:N_CTRL] / tot[:N_CTRL] * 1e6).astype(np.float32); ctrl_mean = ctrl_cpm.mean(0)
G = np.flatnonzero(ctrl_mean > 5.0); C = ctrl_cpm[:, G]
def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1.0)
def mwu_tc(A, B, chunk=800):
    n1, n2 = A.shape[0], B.shape[0]; N = n1 + n2; mu = n1 * (N + 1) / 2.0; p = np.empty(A.shape[1])
    for s in range(0, A.shape[1], chunk):
        e = min(s + chunk, A.shape[1]); Z = np.concatenate([A[:, s:e], B[:, s:e]], 0).astype(np.float64)
        R1 = rankdata(Z, axis=0)[:n1].sum(0); Zs = np.sort(Z, axis=0); T = np.zeros(e - s)
        for j in range(e - s):
            _, c = np.unique(Zs[:, j], return_counts=True); T[j] = (c ** 3 - c).sum()
        var = n1 * n2 / 12.0 * ((N + 1) - T / (N * (N - 1))); z = np.where(var > 0, (R1 - mu) / np.sqrt(np.maximum(var, 1e-300)), 0.0)
        p[s:e] = 2 * norm.sf(np.abs(z))
    return p
rows = []; pos = N_CTRL; gpos = {g: i for i, g in enumerate(genes)}
for t in keep:
    n_p = len(keep[t]); A_all = M[pos:pos + n_p] / tot[pos:pos + n_p] * 1e6; pos += n_p
    for nc in sorted({150, n_p}):
        sub = np.random.default_rng(1).choice(n_p, nc, replace=False) if nc < n_p else np.arange(n_p)
        Ag = A_all[sub][:, G]; q = bh(mwu_tc(Ag, C)); self_mask = np.array([genes[g] != t for g in G])
        lfc = np.log2((Ag.mean(0) + 1e-9) / (ctrl_mean[G] + 1e-9)); sig = (q < ALPHA) & self_mask
        kd = float(lfc[G == gpos[t]][0]) if (t in gpos and gpos[t] in set(G.tolist())) else np.nan
        rows.append(dict(target=t, cells=nc, n_real_tc=int(sig.sum()), frac_neg=float((lfc[sig] < 0).mean()) if sig.any() else np.nan, target_lfc=kd))
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(W / "ka_strong_targets.csv", index=False)
print(f"peak rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3:.2f} GiB")
