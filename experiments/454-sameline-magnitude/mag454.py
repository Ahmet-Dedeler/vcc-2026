#!/usr/bin/env python3
"""454 — does a same-line bank know each target's EFFECT SIZE better than our donor key does? (H1, analytic, no emitter)

The same-line mode (451) blends only DIRECTION: a covered target keeps our row norm (the cap15 composition norm x the 410
key weight). This asks whether the bank's own per-target effect size tracks the truth's better than that product does.
Surface: h1-test truth (323's reference.npz: 85 targets, summed counts), Pan 2025 hESC (30 of 85 covered), the harness's
own bulk rows (435: `native_cap15` then `ps_bulk`, i.e. 410b's key with clip(0.3, 3)).

Truth effect per target: log1p(TS * pseudobulk CPM) difference vs the control pseudobulk, non-panel genes on the shared axis.
Energy is reported raw and jackknife-debiased (reference.npz `jk`; the scorer's own correction, per perturbation).

usage: RUN426_EMITTER=prod uv run python experiments/454-sameline-magnitude/mag454.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("RUN426_EMITTER", "prod")
sys.path.insert(0, str(ROOT / "experiments/435-zs-on-coarse"))
import run435 as M  # noqa: E402

TS = 5e4
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
rgenes = ref["genes"].astype(str)
perts = ref["perts"].astype(str)
sums = ref["sums"].astype(np.float64)
jk = ref["jk"].astype(np.float64)
ci = int(np.flatnonzero(perts == "non-targeting")[0])


def logpb(v):
    return np.log1p(TS * v / v.sum())


b = M.H.builder389()
eff, genes, targets = M.H.native_cap15(b)
eff_bulk = M.ps_bulk(eff, genes, targets)
genes = np.asarray(genes).astype(str)
targets = [str(t) for t in targets]
panel = set(targets)
gi = {g: i for i, g in enumerate(genes)}
shared = [g for g in rgenes if g in gi and g not in panel]
cols_ours = np.array([gi[g] for g in shared])
ri = {g: i for i, g in enumerate(rgenes)}
cols_ref = np.array([ri[g] for g in shared])

# Pan 2025's per-target rows on the challenge axis (450b's pan_rows, uncentred, NaN where unmeasured)
db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"


def _cat(g):
    if isinstance(g, h5py.Group) and "categories" in g:
        return np.asarray(g["categories"][:]).astype(str)[np.asarray(g["codes"][:])]
    return np.asarray(g[:]).astype(str)


with h5py.File(db, "r") as f:
    ds, pg, vg = _cat(f["obs/dataset"]), _cat(f["obs/perturbed_gene_name"]), _cat(f["var/gene_name"])
    sel = np.sort(np.where((ds == "Pan2025_hESC") & np.isin(pg, targets))[0])
    X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
    pt = pg[sel]
vi = {g: i for i, g in enumerate(vg)}
pcol = np.array([vi.get(g, -1) for g in shared])

rows = []
for i, t in enumerate(targets):
    pi = int(np.flatnonzero(perts == t)[0])
    d_true = (logpb(sums[pi]) - logpb(sums[ci]))[cols_ref]
    e_true = float(d_true @ d_true)
    ours = eff_bulk[i, cols_ours]
    rec = dict(target=t, true_e=e_true, true_e_jk=e_true - float(jk[pi]), jk=float(jk[pi]),
               ours_n=float(np.linalg.norm(ours)), native_n=float(np.linalg.norm(eff[i, cols_ours])),
               cos_ours=float(ours @ d_true / (np.linalg.norm(ours) * np.sqrt(e_true) + 1e-30)))
    k = np.flatnonzero(pt == t)
    if len(k):
        pr = X[k[0]]
        pv = np.where(pcol >= 0, pr[np.maximum(pcol, 0)], np.nan)
        m = np.isfinite(pv)
        rec.update(pan_n=float(np.linalg.norm(pv[m])), pan_meas=int(m.sum()),
                   pan_proj=float(pv[m] @ d_true[m] / np.sqrt(e_true)),
                   cos_pan=float(pv[m] @ d_true[m] / (np.linalg.norm(pv[m]) * np.linalg.norm(d_true[m]) + 1e-30)),
                   true_e_on_pan=float(d_true[m] @ d_true[m]))
    rows.append(rec)

cov = [r for r in rows if "pan_n" in r]
print(f"targets {len(rows)}, Pan-covered {len(cov)}; shared non-panel genes {len(shared)}")


def corr(name, x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    lx, ly = np.log(x[ok]), np.log(y[ok])
    print(f"  {name:44s} n {ok.sum():2d}  spearman {spearmanr(lx, ly)[0]:+.3f}  pearson(log) {pearsonr(lx, ly)[0]:+.3f}")


for lab, key in (("raw truth energy", "true_e"), ("jk-debiased truth energy", "true_e_jk")):
    print(f"\n[{lab}] vs predicted norm^2, Pan-covered targets")
    t = [r[key] for r in cov]
    corr("ours (cap15 x 410 key), norm^2", t, [r["ours_n"] ** 2 for r in cov])
    corr("native cap15 (no key), norm^2", t, [r["native_n"] ** 2 for r in cov])
    corr("Pan row norm^2 (measured genes)", t, [r["pan_n"] ** 2 for r in cov])
    corr("geo-mean ours x Pan", t, [r["ours_n"] * r["pan_n"] for r in cov])
    print(f"[{lab}] all 85 targets")
    corr("ours (cap15 x 410 key), norm^2", [r[key] for r in rows], [r["ours_n"] ** 2 for r in rows])
print("\nmedian cos(ours, truth) on covered:", np.median([r["cos_ours"] for r in cov]).round(4),
      " median cos(Pan, truth):", np.median([r["cos_pan"] for r in cov]).round(4))
(HERE / "work").mkdir(exist_ok=True)
(HERE / "work/mag454.json").write_text(json.dumps(rows, indent=1))
print("wrote work/mag454.json")
