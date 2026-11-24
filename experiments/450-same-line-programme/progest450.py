"""450i: which estimator of the same-line programme predicts h1-test's block signs best? Pan rows (h1-test targets out),
compared on genes CLEAN4 also covers (gate-block and fid-like block, as in common450)."""
import runpy
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
P, lfc, gate, block, adj, truth, T, on_axis = (g[k] for k in ("P", "lfc", "gate", "block", "adj", "truth", "T", "on_axis"))
cat = lambda d: d["categories"][:].astype(str)[d["codes"][:]] if isinstance(d, h5py.Group) else d[:].astype(str)
with h5py.File(ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad", "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    rows = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, T))[0])
    X = np.where(f["layers/measure_mask"][rows, :] > 0, f["X"][rows, :], np.nan)
norm = np.sqrt(np.nansum(X ** 2, axis=1))
strong = norm >= np.median(norm)
Z = np.nan_to_num(X)
with np.errstate(invalid="ignore", divide="ignore"):
    est = {
        "P_down all (396 rule-like)": np.nanmean(np.where(np.isnan(X), np.nan, (X < 0).astype(float)), axis=0) - 0.5,
        "P_down, strong half": np.nanmean(np.where(np.isnan(X[strong]), np.nan, (X[strong] < 0).astype(float)), axis=0) - 0.5,
        "P_down, weak half": np.nanmean(np.where(np.isnan(X[~strong]), np.nan, (X[~strong] < 0).astype(float)), axis=0) - 0.5,
        "-mean delta": -np.nanmean(X, axis=0),
        "-median delta": -np.nanmedian(X, axis=0),
        "-mean of unit rows": -np.nanmean(Z / np.maximum(norm[:, None], 1e-9), axis=0),
        "-mean sign weighted by 1/norm": -np.nanmean(np.sign(Z) / np.maximum(norm[:, None], 1e-9), axis=0),
    }
both0 = np.isfinite(P["clean4"])
for name, v in est.items():
    s = on_axis(vg, v)
    sgn = np.where(np.isnan(s), np.nan, np.where(s > 0, -1.0, 1.0))   # positive = down-leaning
    both = both0 & np.isfinite(sgn)
    out = []
    for nm, sub in (("gate-block", gate & block), ("fid-like block", adj & block)):
        mm = sub & both[None, :]
        out.append(f"{nm} {((np.broadcast_to(sgn, lfc.shape) == truth) & mm).sum() / mm.sum():.4f}")
    print(f"  {name:34s} " + "   ".join(out))
mm = gate & block & both0[None, :]
print(f"  {'CLEAN4':34s} gate-block {((np.broadcast_to(P['clean4'], lfc.shape) == truth) & mm).sum() / mm.sum():.4f}")
