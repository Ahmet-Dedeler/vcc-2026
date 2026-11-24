"""450h: the same-line programme in the BULK. h1-test truth rows (log pseudobulk CPM vs control, genes >= 10 CPM, panel
out) carry a shared component (their mean over the 85 targets). Compare that truth mean with our cross-line prediction's
mean (323 cap15) and Pan's mean delta (h1-test targets excluded), and measure what swapping our shared component for
Pan's does to each target's raw cosine and to a scale-optimal squared error (what mse reads)."""
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
eff = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz")
genes_h1 = ref["genes"].astype(str); perts = ref["perts"].astype(str); T = ref["targets"].astype(str)
sums = ref["sums"].astype(float); b_ctrl = ref["b_ctrl"].astype(float); ctrl_cpm = b_ctrl / b_ctrl.sum() * 1e6
cat = lambda d: d["categories"][:].astype(str)[d["codes"][:]] if isinstance(d, h5py.Group) else d[:].astype(str)
with h5py.File(ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad", "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    rows = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, T))[0])
    X = np.where(f["layers/measure_mask"][rows, :] > 0, f["X"][rows, :], np.nan)
pan_mean_all = np.nanmean(X, axis=0)
cg = eff["genes"].astype(str)
panel = set(T)
common = [g for g in genes_h1 if g in set(cg) and g in set(vg) and g not in panel]
ih = {g: i for i, g in enumerate(genes_h1)}; ic = {g: i for i, g in enumerate(cg)}; ip = {g: i for i, g in enumerate(vg)}
gh = np.array([ih[g] for g in common]); gc = np.array([ic[g] for g in common]); gp = np.array([ip[g] for g in common])
keep = (ctrl_cpm[gh] >= 10) & np.isfinite(pan_mean_all[gp])
gh, gc, gp = gh[keep], gc[keep], gp[keep]
pi = {t: i for i, t in enumerate(perts)}; ti = {t: i for i, t in enumerate(eff["targets"].astype(str))}
truth = np.array([np.log(sums[pi[t]][gh] / sums[pi[t]].sum() * 1e6 + 1) - np.log(ctrl_cpm[gh] + 1) for t in T])
ours = np.array([eff["cap15"][ti[t], gc] for t in T])
pan_mean = pan_mean_all[gp]
u = lambda v: v / np.linalg.norm(v)
tm, om = truth.mean(0), ours.mean(0)
print(f"targets {len(T)}, genes {len(gh)}")
print(f"cos(truth mean, ours mean) {u(tm) @ u(om):+.4f}   cos(truth mean, Pan mean) {u(tm) @ u(pan_mean):+.4f}   "
      f"cos(ours mean, Pan mean) {u(om) @ u(pan_mean):+.4f}")
print(f"share of truth energy in its mean: {np.sum(tm**2) * len(T) / np.sum(truth**2):.3f}; ours: {np.sum(om**2) * len(T) / np.sum(ours**2):.3f}")


def sq_err(P):
    """per target, scale-optimal squared error / truth energy (1 = no skill), mean over targets."""
    a = np.sum(P * truth, 1) / np.maximum(np.sum(P * P, 1), 1e-12)
    return float(np.mean(np.sum((truth - np.maximum(a, 0)[:, None] * P) ** 2, 1) / np.sum(truth ** 2, 1)))


cosr = lambda P: float(np.mean(np.sum(P * truth, 1) / (np.linalg.norm(P, axis=1) * np.linalg.norm(truth, axis=1))))
oc = ours - om
print(f"{'prediction':44s} {'raw cos':>8s} {'rel sq err':>11s}")
for name, P in (("ours", ours), ("ours, shared part removed", oc),
                ("ours - shared + Pan mean (ours' mean norm)", oc + u(pan_mean) * np.linalg.norm(om)),
                ("ours - shared + Pan mean (x2)", oc + 2 * u(pan_mean) * np.linalg.norm(om)),
                ("ours - shared + TRUTH mean (oracle)", oc + u(tm) * np.linalg.norm(om))):
    print(f"{name:44s} {cosr(P):8.4f} {sq_err(P):11.4f}")
