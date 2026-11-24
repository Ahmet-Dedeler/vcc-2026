"""How many official-axis genes pass the scorer's 5-CPM control gate on the CHALLENGE's own controls
(context A), streamed from h5py in row chunks. Comparator for K562-GWPS's 8,248/8,248."""
import h5py, numpy as np
p = "data/challenge/context_A.h5ad"
with h5py.File(p, "r") as f:
    X = f["X"]; print("X", dict(X.attrs)); shape = tuple(X.attrs["shape"])
    indptr = X["indptr"][:]
    n = shape[0]; G = shape[1]
    colsum = np.zeros(G); tots = np.zeros(n)
    for s in range(0, n, 2000):
        e = min(s + 2000, n); a, b = indptr[s], indptr[e]
        d = X["data"][a:b].astype(np.float64); ix = X["indices"][a:b]
        rp = indptr[s:e + 1] - a
        row = np.repeat(np.arange(e - s), np.diff(rp))
        t = np.bincount(row, weights=d, minlength=e - s); tots[s:e] = t
        colsum += np.bincount(ix, weights=d / t[row] * 1e6, minlength=G)
    mean_cpm = colsum / n
    print(f"context A: {n:,} cells x {G:,} genes; median UMI {np.median(tots):,.0f}; genes with mean control CPM > 5: {(mean_cpm > 5).sum():,} / {G:,}")
    obs = f["obs"]; print("obs keys", list(obs.keys()))
