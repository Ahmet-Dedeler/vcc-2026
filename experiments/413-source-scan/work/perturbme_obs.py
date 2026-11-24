"""413: Perturb-ME A375 MOI1 h5ad — per-target cell counts from obs only (h5py, no X read)."""
import h5py, json, re, numpy as np, pandas as pd
from overlap import report, to_current, PC
f = h5py.File("dl/PerturbME_A375_HLA_MOI1.h5ad", "r")
obs = f["obs"]
def col(name):
    g = obs[name]
    if isinstance(g, h5py.Group):  # categorical
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
        return pd.Series(cats[g["codes"][:]])
    v = g[:]
    return pd.Series([x.decode() if isinstance(x, bytes) else x for x in v])
print(list(obs.keys()))
sg, cond = col("sgRNA"), col("Condition")
print(cond.value_counts().to_dict()); print(sg.head().tolist())
tgt = sg.str.replace(r"[_-](sg)?\d+$", "", regex=True).str.split("_").str[0]
print(tgt.head().tolist())
out = {}
for c in ["Control", "all"]:
    t = tgt if c == "all" else tgt[cond == c]
    vc = t.value_counts()
    out[c] = {"cells": int(len(t)), "n_targets": int(len(vc)), "median_cells_per_target": float(vc.median()),
              "ge1": report(f"PerturbME {c} >=1", vc.index.tolist())["panel_exact"],
              "ge20": report(f"PerturbME {c} >=20", vc[vc >= 20].index.tolist())["panel_exact"],
              "ge50": report(f"PerturbME {c} >=50", vc[vc >= 50].index.tolist())["panel_exact"]}
    pan = vc[vc.index.map(to_current).isin(set(PC))]
    out[c]["panel_cells_median"] = float(pan.median()) if len(pan) else 0.0
print(json.dumps(out, indent=1))
json.dump(out, open("perturbme_obs_overlap.json", "w"), indent=1)
