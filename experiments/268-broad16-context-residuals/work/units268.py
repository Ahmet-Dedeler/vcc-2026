"""Unit table for RESULTS.md from the per-line build summaries (working lines) + metadata."""
import json, glob, numpy as np, pandas as pd
from load268 import *
rows = []
for l in WORKING:
    s = json.load(open(W / "halves" / f"{l}.json"))["summary"]; z = np.load(W / "halves" / f"{l}.npz", allow_pickle=True); nc = z["ncell"]; tl = list(z["targets"])
    real = np.array([nc[i, 0] for i, t in enumerate(tl) if t not in CONTROLS and nc[i, 0] > 0]); ctrl = np.array([nc[i, 0] for i, t in enumerate(tl) if t in CONTROLS])
    rows.append({"line": l, "batch": META.loc[l, "batch"], "lineage": META.loc[l, "OncotreeLineage"], "h5 cells": s["h5_cells"], "cells kept": s["pass_qc_cells"], "KOs": s["targets_with_cells"],
                 "cells/KO median": int(np.median(real)), "min": int(real.min()), "max": int(real.max()), "control cells (10 OR KOs)": int(ctrl.sum()), "UMI/cell": int(s["lib_mean"]),
                 "targets < 6 cells": int((real < 6).sum()), "KA exact": f"{s['ka_exact_targets']}/{s['ka_targets']}"})
df = pd.DataFrame(rows)
def md(df):
    cols = list(df.columns); out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows(): out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)
print(md(df))
print(); print("total kept cells", df["cells kept"].sum(), "| targets usable (>=6 cells on all 12 lines):", len(usable_targets(load_halves(WORKING)[0], WORKING, 6)), "of", len(COMPLETE), "complete")
