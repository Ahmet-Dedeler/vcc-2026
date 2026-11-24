"""Real-control group sums and jackknife per challenge context.

`data/challenge/context_{A,B,C}.h5ad` ARE the reference bundle's control cells (18,400 =
46 ntc_id x 400, manifest.json), and the metric uses the reference control as the origin on
both sides. So b_ctrl and C_ctrl are exactly computable for the scored panel without holding
any of the answers.
"""
import sys, time
import numpy as np, anndata as ad
from scipy.sparse import csr_matrix, issparse
from cell_eval2.moments import jackknife_correction
TS = 50_000.0
out = {}
for c in "ABC":
    t0 = time.time()
    A = ad.read_h5ad(f"~/Code/virtual-cell-challenge-2026/data/challenge/context_{c}.h5ad")
    X = A.X.tocsr() if issparse(A.X) else csr_matrix(A.X)
    n = X.shape[0]
    out[f"sum_{c}"] = np.asarray(X.sum(0)).ravel().astype(np.int64)
    out[f"n_{c}"] = np.float64(n)
    out[f"jk_{c}"] = jackknife_correction(X, np.zeros(n, np.intp), 1, TS, chunk=64)[0]
    out["genes"] = np.asarray([str(g) for g in A.var_names])
    print(f"{c}: n={n} jk={out[f'jk_{c}']:.6f}  {time.time()-t0:.0f}s", flush=True)
np.savez_compressed("~/Code/virtual-cell-challenge-2026/experiments/096-mse-floor/work/ctrl-stats.npz", **out)
print("done")
