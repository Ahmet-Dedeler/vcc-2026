"""413: exact panel overlap for (a) head-to-head CRISPRko/CRISPRi A549 d7 library (Zenodo 21136747) and (b) DeSCOPE PC9 h5ad obs (Zenodo 21622013)."""
import h5py, json, numpy as np, pandas as pd
from overlap import report
out = {}
g = [l.strip().split("_")[0] for l in open("dl/h2h_A549_targets.txt") if l.strip()]
out["h2h_A549"] = report("head-to-head A549 d7 guides", g) | {"n_guides": len(g)}
f = h5py.File("dl/PC9.h5ad", "r")
print("PC9 obs keys", list(f["obs"].keys())[:30]); print("X", dict(f["X"].attrs) if isinstance(f["X"], h5py.Group) else f["X"].shape)
print("var n", len(f["var"][f["var"].attrs.get("_index", "_index")]))
def col(name):
    x = f["obs"][name]
    if isinstance(x, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in x["categories"][:]]); return pd.Series(cats[x["codes"][:]])
    return pd.Series([v.decode() if isinstance(v, bytes) else v for v in x[:]])
for k in f["obs"].keys():
    try:
        s = col(k)
        if s.dtype == object and 1 < s.nunique() < 50000: print(" ", k, s.nunique(), s.value_counts().head(6).to_dict())
    except Exception as e: pass
json.dump(out, open("small_sets_overlap.json", "w"), indent=1)
print(out["h2h_A549"]["n_unique"], out["h2h_A549"]["panel_exact"], out["h2h_A549"]["panel_hits"])
