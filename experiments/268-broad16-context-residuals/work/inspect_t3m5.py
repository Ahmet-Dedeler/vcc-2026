import numpy as np, pandas as pd, h5py, scanpy as sc, time
t0 = time.time()
g = pd.read_csv("sc/T3M5_CRISPR.filt", index_col=0)
print("CRISPR.filt shape", g.shape, "dtype", g.dtypes.iloc[0], f"{time.time()-t0:.0f}s")
print("row names (first 12):", list(g.index[:12]))
print("row name patterns:", sorted(set(i.split("_")[0] if "_" in i else i.split("-")[0] for i in g.index))[:20], "n unique prefixes", len(set(i.rsplit("_",1)[0] for i in g.index)))
vals = g.to_numpy()
print("value stats: min", vals.min(), "max", vals.max(), "frac nonzero", (vals > 0).mean(), "integral", np.all(np.mod(vals, 1) == 0))
nz_per_cell = (vals > 0).sum(0)
print("guides>0 per cell:", np.bincount(nz_per_cell)[:8], "cells", vals.shape[1])
print("per-cell total guide UMI quantiles:", np.quantile(vals.sum(0), [0, .1, .5, .9, 1]))
with h5py.File("sc/T3M5_GE_filtered_feature_bc_matrix", "r") as h:
    def show(name, obj):
        if isinstance(obj, h5py.Dataset): print(" ", name, obj.shape, obj.dtype)
    h.visititems(show)
    ft = h["matrix/features/feature_type"][:]
    print("feature types:", {k: int(v) for k, v in zip(*np.unique(ft, return_counts=True))})
a = sc.read_10x_h5("sc/T3M5_GE_filtered_feature_bc_matrix", gex_only=False)
print("adata", a.shape, a.var["feature_types"].value_counts().to_dict(), f"{time.time()-t0:.0f}s")
print("obs names sample:", list(a.obs_names[:3]), "var names sample:", list(a.var_names[:3]), "var unique", a.var_names.is_unique)
bc_h5 = set(a.obs_names); bc_g = set(g.columns)
print("barcodes h5", len(bc_h5), "guide-file", len(bc_g), "overlap", len(bc_h5 & bc_g))
bc_h5s = set(b.split("-")[0] for b in a.obs_names); bc_gs = set(b.split("-")[0] for b in g.columns)
print("overlap after stripping suffix", len(bc_h5s & bc_gs))
