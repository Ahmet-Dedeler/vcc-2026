"""450d: per-target DE signs from a same-line bank's own rows (the 402 question on H1, cross-study).
For h1-test targets Pan covers: sign accuracy on each target's gate genes (q < 0.05) of
  ours    = sign of the cross-line composition (323 cap15),
  pan_row = sign of Pan's row for that target,
  pan_prog= Pan's programme (P_down over its other rows, as in sameline450),
  clean4  = the shipped per-gene prior.
"""
import runpy
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
P, lfc, gate, block, truth, T, G, gidx = (g[k] for k in ("P", "lfc", "gate", "block", "truth", "T", "G", "gidx"))
eff = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz")
cg = {x: i for i, x in enumerate(eff["genes"].astype(str))}
ours = np.full(lfc.shape, np.nan)
col = np.array([cg.get(x, -1) for x in G])
ti = {t: i for i, t in enumerate(eff["targets"].astype(str))}
for i, t in enumerate(T):
    r = eff["cap15"][ti[t]]
    ours[i] = np.where(col >= 0, np.sign(r[np.maximum(col, 0)]), np.nan)
ours[ours == 0] = np.nan


def cat(d):
    return d["categories"][:].astype(str)[d["codes"][:]] if isinstance(d, h5py.Group) else d[:].astype(str)


with h5py.File(ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad", "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    sel = np.sort(np.where((ds == "Pan2025_hESC") & np.isin(pg, T))[0])
    Xp = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
    pt = pg[sel]
vi = {x: i for i, x in enumerate(vg)}
pcol = np.array([vi.get(x, -1) for x in G])
panrow = np.full(lfc.shape, np.nan)
for i, t in enumerate(T):
    rows = np.where(pt == t)[0]
    if len(rows):
        r = np.nanmean(Xp[rows], axis=0)
        panrow[i] = np.where(pcol >= 0, np.sign(r[np.maximum(pcol, 0)]), np.nan)
panrow[panrow == 0] = np.nan
cov_t = np.array([np.isfinite(panrow[i]).any() for i in range(len(T))])
print(f"targets with a Pan row: {cov_t.sum()}")
for name, sub in (("gate (covered targets)", gate & cov_t[:, None]), ("gate-block (covered)", gate & block & cov_t[:, None])):
    both = sub & np.isfinite(ours) & np.isfinite(panrow)
    w = np.abs(np.where(both, lfc, 0))
    out = []
    for nm, pr in (("ours", ours), ("pan_row", panrow), ("pan_prog", np.broadcast_to(P["pan_pdown"], lfc.shape)),
                   ("clean4", np.broadcast_to(P["clean4"], lfc.shape))):
        ok = both & np.isfinite(pr)
        hit = (pr == truth) & ok
        out.append(f"{nm} {hit.sum() / ok.sum():.4f}/{(w * hit).sum() / (w * ok).sum():.4f}")
    print(f"{name:24s} n {int(both.sum()):6d}  " + "  ".join(out) + "   (unweighted/|lfc|-weighted)")
# per-target, ours vs pan_row on gate
d = []
for i in np.where(cov_t)[0]:
    m = gate[i] & np.isfinite(ours[i]) & np.isfinite(panrow[i])
    if m.sum() >= 10:
        d.append((np.mean(panrow[i][m] == truth[i][m]), np.mean(ours[i][m] == truth[i][m])))
d = np.array(d)
print(f"per target (gate >= 10): n {len(d)}; pan_row {d[:, 0].mean():.4f} ours {d[:, 1].mean():.4f}; pan better on "
      f"{int((d[:, 0] > d[:, 1]).sum())}")
