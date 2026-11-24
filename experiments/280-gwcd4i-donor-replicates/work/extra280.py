"""Extract a random sample of NON-PANEL targeting pseudobulks (Stim48hr only) so the residual centring can be
re-computed on a genome-wide sample instead of the 297 panel targets -- a robustness check on the centring."""
import numpy as np, h5py, time
from parts_file import PartsFile
z = np.load("obs280.npz", allow_pickle=True)
indptr = z["indptr"]; pert = z["pert"].astype(str); gtype = z["gtype"].astype(str); cond = z["cond"].astype(str)
ROOT = "~/Code/virtual-cell-challenge-2026"
panel = set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:])
cand = sorted(set(pert[(gtype == "targeting") & (cond == "Stim48hr")]) - panel)
rng = np.random.default_rng(20280906)
pick = set(rng.choice(cand, 500, replace=False))
sel = np.where((gtype == "targeting") & (cond == "Stim48hr") & np.isin(pert, list(pick)))[0]
print(f"{len(cand)} non-panel Stim48hr targets; sampled 500 -> {len(sel)} rows")
pf = PartsFile([f"parts/p{i}" for i in range(8)]); f = h5py.File(pf, "r")
D, I = f["X/data"], f["X/indices"]
out = np.lib.format.open_memmap("extra280_X.npy", mode="w+", dtype=np.float32, shape=(len(sel), len(z["genes"])))
runs = np.split(sel, np.where(np.diff(sel) != 1)[0] + 1); w = 0; t0 = time.time()
for run in runs:
    a, b = int(indptr[run[0]]), int(indptr[run[-1] + 1])
    dat = np.asarray(D[a:b]); idx = np.asarray(I[a:b]); off = indptr[run[0]]
    for r in run:
        s, e = int(indptr[r] - off), int(indptr[r + 1] - off)
        out[w, idx[s:e]] = dat[s:e]; w += 1
out.flush(); np.save("extra280_rows.npy", sel)
tot = z["total"][sel]
X = np.load("extra280_X.npy", mmap_mode="r")
mx = max(float(np.abs(np.asarray(X[lo:lo+2000], dtype=np.float64).sum(1) - tot[lo:lo+2000]).max())
         for lo in range(0, len(sel), 2000))
print(f"done {time.time()-t0:.0f}s, {len(sel)} rows, max|rowsum - obs.total| = {mx:g}")
