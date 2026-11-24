"""Verify the extracted bank against the pseudobulk obs: row sums must equal obs.total_counts exactly."""
import numpy as np
z = np.load("obs280.npz", allow_pickle=True)
sel = np.load("bank280_rows.npy")
X = np.load("bank280_X.npy", mmap_mode="r")
tot = z["total"][sel]
rng = np.random.default_rng(1)
idx = np.sort(rng.choice(len(sel), 400, replace=False))
s = np.asarray(X[idx], dtype=np.float64).sum(1)
d = np.abs(s - tot[idx])
print(f"400 random rows: max|rowsum - obs.total_counts| = {d.max():.6g}  ({int((d==0).sum())}/400 exact)")
# full check in blocks
mx = 0.0; nex = 0
for lo in range(0, len(sel), 2000):
    b = np.asarray(X[lo:lo+2000], dtype=np.float64).sum(1)
    dd = np.abs(b - tot[lo:lo+2000]); mx = max(mx, float(dd.max())); nex += int((dd == 0).sum())
print(f"ALL {len(sel)} rows: max|d| = {mx:.6g}; exact on {nex}/{len(sel)}")
nz = (np.asarray(X[:2000]) != 0).sum(1)
print(f"nnz per row (first 2000): median {np.median(nz):.0f} of {X.shape[1]}")
