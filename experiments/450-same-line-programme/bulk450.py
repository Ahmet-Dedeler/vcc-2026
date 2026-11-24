"""450c: a same-line bank's per-TARGET rows as a bulk donor (the October question, row 220 of DEF-RUNBOOK §3).

Pan 2025 hESC covers 30 of h1-test's 85 targets. For those, compare the h1-test truth effect (log pseudobulk CPM of the
target's 400 cells vs the controls, genes >= 10 control CPM, the 85 panel genes out) with:
  ours  = the cross-line transfer composition on H1 (323's `cap15`, K562-GWPS / HCT116 / HEK293T whitened, no H1 data);
  pan   = Pan's own row for the target (same line, different study);
  blend = unit(ours) + w * unit(pan).
Reports per-target cosine and a pds-like rank (is the prediction nearest its own truth among the 30).
"""
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
eff = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz")
genes_h1 = ref["genes"].astype(str); perts = ref["perts"].astype(str); T = ref["targets"].astype(str)
sums = ref["sums"].astype(float); b_ctrl = ref["b_ctrl"].astype(float)
ctrl_cpm = b_ctrl / b_ctrl.sum() * 1e6


def cat(d):
    return d["categories"][:].astype(str)[d["codes"][:]] if isinstance(d, h5py.Group) else d[:].astype(str)


with h5py.File(ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad", "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    sel = np.where((ds == "Pan2025_hESC") & np.isin(pg, T))[0]
    Xp = f["X"][np.sort(sel), :]
    Mp = f["layers/measure_mask"][np.sort(sel), :]
    pan_t = pg[np.sort(sel)]
Xp = np.where(Mp > 0, Xp, np.nan)
cg = eff["genes"].astype(str); ours_all = eff["cap15"]
panel = set(T)
common = [g for g in genes_h1 if g in set(cg) and g in set(vg) and g not in panel]
ih = {g: i for i, g in enumerate(genes_h1)}; ic = {g: i for i, g in enumerate(cg)}; ip = {g: i for i, g in enumerate(vg)}
gh = np.array([ih[g] for g in common]); gc = np.array([ic[g] for g in common]); gp = np.array([ip[g] for g in common])
keep = ctrl_cpm[gh] >= 10
gh, gc, gp = gh[keep], gc[keep], gp[keep]
covered = [t for t in T if t in set(pan_t)]
ti = {t: i for i, t in enumerate(T)}; pi = {t: i for i, t in enumerate(perts)}
truth, ours, pan = [], [], []
for t in covered:
    s = sums[pi[t]]
    cpm = s / s.sum() * 1e6
    truth.append(np.log(cpm[gh] + 1.0) - np.log(ctrl_cpm[gh] + 1.0))
    ours.append(ours_all[ti[t], gc])
    rows = np.where(pan_t == t)[0]
    pan.append(np.nanmean(Xp[rows][:, gp], axis=0))
truth, ours, pan = (np.nan_to_num(np.asarray(v)) for v in (truth, ours, pan))


def unit(M):
    return M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-12)


def cosrow(A, B):
    return np.sum(unit(A) * unit(B), axis=1)


def pds_like(P, Tr):
    """Fraction of targets whose prediction is closer (cosine) to its own truth than to the average other truth."""
    C = unit(P) @ unit(Tr).T
    own = np.diag(C)
    return float(np.mean([(C[i, np.arange(len(C)) != i] < own[i]).mean() for i in range(len(C))]))


# mean-centred versions (pds-style: remove the across-target mean of each side)
tc = truth - truth.mean(0); oc = ours - ours.mean(0); pc = pan - pan.mean(0)
print(f"targets covered by Pan: {len(covered)} of {len(T)}; genes {len(gh)} (>= 10 control CPM, panel out)")
print(f"{'predictor':22s} {'cos raw':>9s} {'cos centred':>12s} {'pds-like (centred)':>19s}")
for name, P, Pc in (("ours (cross-line)", ours, oc), ("pan (same line)", pan, pc)):
    print(f"{name:22s} {cosrow(P, truth).mean():9.4f} {cosrow(Pc, tc).mean():12.4f} {pds_like(Pc, tc):19.4f}")
for w in (0.25, 0.5, 1.0, 2.0, 4.0):
    B = unit(ours) + w * unit(pan); Bc = unit(oc) + w * unit(pc)
    print(f"{'blend w=' + str(w):22s} {cosrow(B, truth).mean():9.4f} {cosrow(Bc, tc).mean():12.4f} {pds_like(Bc, tc):19.4f}")
print("per-target cos (centred) ours vs pan:", np.round(cosrow(oc, tc), 3).tolist()[:10], np.round(cosrow(pc, tc), 3).tolist()[:10])
print("pan better on", int((cosrow(pc, tc) > cosrow(oc, tc)).sum()), "of", len(covered))
