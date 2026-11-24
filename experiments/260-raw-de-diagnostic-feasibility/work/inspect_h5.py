"""Inspect an h5ad's structure with h5py only. Never touches X beyond metadata."""
import sys, h5py, numpy as np
p = sys.argv[1]
with h5py.File(p, "r") as f:
    print("top keys:", list(f.keys()))
    X = f["X"]
    if isinstance(X, h5py.Group):
        print("X group attrs:", dict(X.attrs)); 
        for k in X: print("  X/", k, X[k].shape, X[k].dtype, X[k].chunks)
    else:
        print("X dataset:", X.shape, X.dtype, X.chunks)
    for grp in ("obs", "var", "layers", "uns", "obsm", "raw"):
        if grp in f:
            g = f[grp]
            print(f"{grp} attrs:", {k: (v if not hasattr(v, 'tolist') else v.tolist()) for k, v in g.attrs.items()})
            for k in g:
                o = g[k]
                if isinstance(o, h5py.Group):
                    sub = {kk: (o[kk].shape, str(o[kk].dtype)) for kk in o}
                    print(f"  {grp}/{k}: group {dict(o.attrs).get('encoding-type')} {sub}")
                else:
                    print(f"  {grp}/{k}: {o.shape} {o.dtype}")
