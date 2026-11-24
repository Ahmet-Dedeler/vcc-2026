"""700: log1p(CP10k) control h5ad of the h1-test harness for a PIE controls-only dir (context ARC_H1, dataset h1test).

The dataset name differs from PIE's arc_vcc_25, so PIE serves no same-context evidence (keyed by (dataset, context)); the
context text is H1's (keyed by context name). h1-test's 85 targets are not in PIE's VCC25 train/val split (190 perts).
"""
import sys, numpy as np, anndata as ad, scipy.sparse as sp
src, out = sys.argv[1], sys.argv[2]
a = ad.read_h5ad(src)
m = (a.obs["target_gene"] == "non-targeting").to_numpy()
X = a.X[m]
X = X.tocsr().astype(np.float32) if sp.issparse(X) else sp.csr_matrix(np.asarray(X, dtype=np.float32))
lib = np.asarray(X.sum(1)).ravel()
X = sp.diags((1e4 / lib).astype(np.float32)) @ X
X.data = np.log1p(X.data)
obs = a.obs.loc[m, ["target_gene"]].copy()
obs["context"] = "ARC_H1"
ad.AnnData(X=X.tocsr(), obs=obs, var=a.var[[]].copy()).write_h5ad(out)
print(X.shape, float(np.median(lib)), flush=True)
