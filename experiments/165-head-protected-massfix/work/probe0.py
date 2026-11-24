import h5py, numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
for p in (f"{ROOT}/data/harness/h1-test.h5ad",
          f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"):
    with h5py.File(p,"r") as f:
        X=f["X"]
        print(p.split("/")[-1], "encoding", dict(X.attrs).get("encoding-type"), "shape", dict(X.attrs).get("shape"))
        def idx(g):
            v=g["_index"]
            if isinstance(v,h5py.Group): v=v["values"]
            return v
        vn=idx(f["var"]); on=idx(f["obs"])
        print("   n_var",vn.shape,"n_obs",on.shape)
        print("   obs keys", list(f["obs"].keys()))
