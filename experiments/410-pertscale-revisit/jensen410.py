"""The Jensen bias of a 400-cell group's log-pseudobulk against the full control, from the challenge's own controls."""
import sys
from pathlib import Path
import numpy as np
import anndata as ad
ROOT = Path(__file__).resolve().parents[2]
TS = 5e4
rng = np.random.default_rng(0)
for ctx in ("A", "B", "C"):
    A = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
    X = A.X.tocsr()
    tot = np.asarray(X.sum(0)).ravel()
    b_full = np.log1p(TS * tot / tot.sum())
    D = []
    for _ in range(60):
        idx = rng.choice(X.shape[0], 400, replace=False)
        s = np.asarray(X[idx].sum(0)).ravel()
        D.append(np.log1p(TS * s / s.sum()) - b_full)
    D = np.array(D)
    J = D.mean(0); sd = D.std(0)
    exp400 = tot / tot.sum() * 400 * 20000
    bins = [0, 0.5, 2, 5, 20, 100, 1000, 1e12]
    print(f"context {ctx}: J (mean over 60 subsets) by expected count in a 400-cell group")
    for k in range(len(bins) - 1):
        m = (exp400 >= bins[k]) & (exp400 < bins[k + 1])
        print(f"   [{bins[k]:>6g},{bins[k+1]:>7g}) n {m.sum():5d}  J mean {J[m].mean():+.5f}  per-subset sd {sd[m].mean():.4f}")
    np.save(ROOT / f"experiments/410-pertscale-revisit/work/jensen_{ctx}.npy", J)
    del A, X
