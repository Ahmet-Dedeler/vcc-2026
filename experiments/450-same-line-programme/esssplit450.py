"""450f: is the same-line programme carried by essential perturbations too? (Nadig Jurkat, the bank context A can use, is
an essential screen.) Pan's programme recomputed from (a) its K562-essential-list perturbations only, (b) the rest only,
scored on h1-test like sameline450 (gate-block / fid-like block), common coverage with CLEAN4."""
import runpy
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
P, lfc, gate, block, adj, truth, T, G, on_axis = (g[k] for k in ("P", "lfc", "gate", "block", "adj", "truth", "T", "G", "on_axis"))
ess = set(l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt"))
cat = lambda d: d["categories"][:].astype(str)[d["codes"][:]] if isinstance(d, h5py.Group) else d[:].astype(str)
with h5py.File(ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad", "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    rows = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, T))[0])
    X = np.where(f["layers/measure_mask"][rows, :] > 0, f["X"][rows, :], np.nan)
    names = pg[rows]
is_e = np.array([n in ess for n in names])
print(f"Pan rows {len(names)}: K562-essential-list {is_e.sum()}, other {(~is_e).sum()}")
for lab, m in (("essential only", is_e), ("non-essential only", ~is_e), ("all", np.ones_like(is_e))):
    with np.errstate(invalid="ignore"):
        pdn = np.nanmean(np.where(np.isnan(X[m]), np.nan, (X[m] < 0).astype(float)), axis=0)
    pr = on_axis(vg, pdn)
    s = np.where(np.isnan(pr), np.nan, np.where(pr > 0.5, -1.0, 1.0))
    both = np.isfinite(s) & np.isfinite(P["clean4"])
    out = []
    for nm, sub in (("gate-block", gate & block), ("fid-like block", adj & block)):
        mm = sub & both[None, :]
        a = ((np.broadcast_to(s, lfc.shape) == truth) & mm).sum() / mm.sum()
        c = ((np.broadcast_to(P["clean4"], lfc.shape) == truth) & mm).sum() / mm.sum()
        out.append(f"{nm} {a:.4f} (CLEAN4 {c:.4f})")
    print(f"  {lab:20s} n={int(m.sum()):5d}  " + "   ".join(out))
