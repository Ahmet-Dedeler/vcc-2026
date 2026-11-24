"""700: write log1p(CP10k) control-cell h5ads for PIE's controls-only prep (one per context, CSR float32).

PIE's `pie process` normalizes raw counts to CP10k + natural log1p; prep then averages X over control cells to get
ctrl_means. Context labels follow PIE's context_text keys: A -> jurkat (Replogle's Jurkat text), B -> CVCL_0030 (HeLa),
C -> CVCL_1108 (CAL-33); both new texts were added by 700 (data/external/pie/local_sources/context_text).
usage: .venv/bin/python work/make_ctrl_h5ad700.py <bundle dir> <out dir> [A=jurkat B=CVCL_0030 C=CVCL_1108]
"""
import sys, numpy as np, anndata as ad, scipy.sparse as sp
bundle, out = sys.argv[1], sys.argv[2]
mapping = dict(a.split("=") for a in (sys.argv[3:] or ["A=jurkat", "B=CVCL_0030", "C=CVCL_1108"]))
for ctx, name in mapping.items():
    a = ad.read_h5ad(f"{bundle}/context_{ctx}.h5ad")
    X = a.X.tocsr().astype(np.float32) if sp.issparse(a.X) else sp.csr_matrix(a.X.astype(np.float32))
    lib = np.asarray(X.sum(1)).ravel()
    X = sp.diags((1e4 / lib).astype(np.float32)) @ X
    X.data = np.log1p(X.data)
    obs = a.obs[["target_gene"]].copy()
    obs["context"] = name
    assert (obs["target_gene"] == "non-targeting").all()
    b = ad.AnnData(X=X.tocsr(), obs=obs, var=a.var[[]].copy())
    b.write_h5ad(f"{out}/ctrl_{ctx}_{name}.h5ad")
    print(ctx, name, b.shape, X.nnz, float(np.median(lib)), flush=True)
