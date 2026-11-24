import h5py, sys, numpy as np

def info(p):
    try:
        f = h5py.File(p, "r")
    except Exception as e:
        print(p, "ERR", e); return
    with f:
        X = f["X"]
        if isinstance(X, h5py.Group):
            shape = tuple(X.attrs.get("shape", ("?","?")))
            enc = X.attrs.get("encoding-type")
            dt = X["data"].dtype
        else:
            shape = X.shape; enc="dense"; dt=X.dtype
        obs = f["obs"]
        cols = list(obs.attrs.get("column-order", []))
        print(f"\n== {p}")
        print("  X:", enc, shape, dt)
        print("  obs cols:", cols)
        try:
            vidx = f["var"].attrs.get("_index")
            vn = f["var"][vidx][:5]
            print("  n_var:", f["var"][vidx].shape, "first:", [x.decode() if isinstance(x,bytes) else x for x in vn])
        except Exception as e:
            print("  var err", e)
        for c in cols:
            node = obs[c]
            if isinstance(node, h5py.Group) and "categories" in node:
                cats = node["categories"][:]
                cats = [x.decode() if isinstance(x, bytes) else str(x) for x in cats]
                print(f"  obs[{c}]: categorical n={len(cats)} sample={cats[:4]}{'...' if len(cats)>4 else ''}")
            else:
                print(f"  obs[{c}]: {node.dtype} shape={node.shape}")

for p in sys.argv[1:]:
    info(p)
