#!/usr/bin/env python3
"""462l — the H1 six-member price of Arc's OWN programme as the block-sign prior (an UPPER reference for 458d/458: same
line AND lab as h1-test, so leaky by design; h1-test's 85 targets are excluded from the programme rows). Derived from
450's run450.py (the same wrapper around 435's harness); only the prior source differs.
  arcfull   data/donors/h1_2025_full (458d's bank: all public cells) P_down > 0.5 -> down; uncovered -> -1
  arcbank   data/donors/h1_2025 (458's bank: the 400-cell harness)
  clean4    the shipped prior (known answer: must equal 435's stats bit for bit)
usage: RUN426_EMITTER=prod uv run python run462.py <clean4|arcfull|arcbank> <run435 args...>   -> work/h1_<prior>/

--- 450's original docstring follows ---
450 — the H1 six-member price of replacing the CLEAN4 block-sign prior with a SAME-LINE bank's per-gene prior.

This is 435's harness (`run435.py`: production emitter `rc_emit_k` edb5b8dc, 410b ballast, coarse MDE) with ONE thing
swapped: `run400.prior_signs_full`, the per-gene sign the emitter steers the CPM block to. On H1 the same-line bank is
Pan 2025's hESC CRISPRi screen (GSE295214 pseudobulk database), a different lab, library and chemistry from the 2025 VCC
H1 data that h1-test comes from; h1-test's own 85 targets are excluded from its programme. This is the H1 analogue of the
builder's `clean4q80ja` (context A reads Nadig Jurkat's per-gene P_down instead of CLEAN4's).

usage: RUN426_EMITTER=prod uv run python run450.py <clean4|pan|pan_c4fill> <run435 args...>
  clean4      the shipped prior through this wrapper (known answer: must equal 435's stats bit for bit)
  pan         Pan hESC P_down > 0.5 -> down, < 0.5 -> up; genes Pan does not measure -> -1 (run400's uncovered default)
  pan_c4fill  Pan where measured, CLEAN4 elsewhere
Outputs go to work/h1_<prior>/ (435's files are never touched).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/450-same-line-programme"))
os.environ.setdefault("RUN426_EMITTER", "prod")
PRIOR = sys.argv.pop(1)
assert PRIOR in ("clean4", "arcfull", "arcbank"), PRIOR
sys.path.insert(0, str(ROOT / "experiments/435-zs-on-coarse"))
import run435 as M  # noqa: E402

WORK = HERE / "work" / f"h1_{PRIOR}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_clean4 = M.H.prior_signs_full
H1_TARGETS = set(np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz",
                         allow_pickle=True)["targets"].astype(str))


def _cat(ds):
    return ds["categories"][:].astype(str)[ds["codes"][:]] if isinstance(ds, h5py.Group) else ds[:].astype(str)


def pan_pdown():
    """Pan 2025 hESC P_down per gene (fraction of its perturbations with delta < 0), h1-test's targets excluded."""
    db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
    with h5py.File(db, "r") as f:
        dataset, pg, vg = _cat(f["obs/dataset"]), _cat(f["obs/perturbed_gene_name"]), _cat(f["var/gene_name"])
        rows = np.sort(np.where((dataset == "Pan2025_hESC") & ~np.isin(pg, list(H1_TARGETS)))[0])
        X = f["X"][rows, :]
        Mk = f["layers/measure_mask"][rows, :]
    X = np.where(Mk > 0, X, np.nan)
    with np.errstate(invalid="ignore"):
        pd = np.nanmean(np.where(np.isnan(X), np.nan, (X < 0).astype(float)), axis=0)
    return vg, pd, len(rows)


def arc_pdown(bank):
    """Arc H1 P_down per gene from a bank dir (delta.npy rows), h1-test's targets excluded."""
    d = ROOT / "data/donors" / bank
    delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    genes = np.asarray([l.strip() for l in open(d / "genes.txt")])
    rows = [i for i, p in enumerate(perts) if p not in H1_TARGETS]
    M = np.asarray(delta[rows], dtype=np.float64); nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    return genes, pd, len(rows)


def prior_signs_full(genes_full):
    base, cov_c4 = _clean4(genes_full)
    if PRIOR == "clean4":
        return base, cov_c4
    vg, pd, n = arc_pdown("h1_2025_full" if PRIOR == "arcfull" else "h1_2025")
    pos = {g: j for j, g in enumerate(vg)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    val = np.where(col >= 0, pd[np.maximum(col, 0)], np.nan)
    covered = np.isfinite(val)
    fill = np.full(len(genes_full), -1.0)
    out = np.where(covered, np.where(val > 0.5, -1.0, 1.0), fill)
    print(f"[462l] prior {PRIOR}: Arc rows {n}, covered {covered.mean():.4f} of the axis; differs from CLEAN4 on "
          f"{(out != base).mean():.4f}", flush=True)
    return out, float(covered.mean())


M.H.prior_signs_full = prior_signs_full

if __name__ == "__main__":
    M.main()
