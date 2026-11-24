#!/usr/bin/env python3
"""614 forward question (coordinator): which cancer-line Perturb-seq on disk has real power (>= 100 cells/target AND
>= 20k UMI/cell)? Measures each h5ad's perturbation label column, cells per perturbation (median over perturbations,
controls excluded), and median UMI per cell. Totals come from an obs count column if one exists, else from a contiguous
block of rows. No data is written except work/inventory614.json."""
import json
import sys
from pathlib import Path
import h5py
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
FILES = {
    "nadig_jurkat": "data/external/nadig2025/jurkat_raw_singlecell.h5ad",
    "nadig_hepg2": "data/external/nadig2025/hepg2_raw_singlecell.h5ad",
    "replogle_k562_essential": "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
    "replogle_k562_gwps": "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad",
    "replogle_rpe1": "data/external/scperturb/ReplogleWeissman2022_rpe1.h5ad",
    "frangieh_melanoma": "data/external/scperturb/FrangiehIzar2021_RNA.h5ad",
    "myllia_a549_ch1": "data/external/myllia-a549/GSE337804_A549_ch_1.RDS.h5ad",
    "altair_a549": "data/external/altair-a549/altair-scrna.h5ad",
}
LABELS = ("gene", "perturbation", "target_gene", "gene_target", "guide_target", "perturbation_name", "target", "gene_symbol")
COUNTS = ("ncounts", "UMI_count", "total_counts", "n_counts", "nCount_RNA", "nUMI", "UMI count")
CTRL = ("non-targeting", "control", "NT", "non_targeting", "Non-targeting", "nontargeting", "NTC", "safe-targeting")


def col(obs, names):
    for n in names:
        if n in obs:
            o = obs[n]
            if isinstance(o, h5py.Group) and "categories" in o:
                return n, o["categories"].asstr()[:][o["codes"][:]]
            if isinstance(o, h5py.Dataset):
                v = o[:]
                return n, (v.astype(str) if v.dtype.kind in "SOU" else v)
    return None, None


def umis(h, ncells):
    X = h["X"]
    if isinstance(X, h5py.Group):                       # CSR / CSC
        ip = X["indptr"][:]
        if len(ip) - 1 == ncells:
            lo, hi = 0, min(ncells, 5000)
            d = X["data"][ip[lo]:ip[hi]]
            rows = np.repeat(np.arange(hi - lo), np.diff(ip[lo:hi + 1]))
            return np.bincount(rows, weights=d, minlength=hi - lo)
        return None
    hi = min(ncells, 4000)                               # dense: one contiguous block
    return np.asarray(X[0:hi, :], dtype=np.float64).sum(1)


out = {}
for name, rel in FILES.items():
    p = ROOT / rel
    if not p.exists():
        out[name] = {"missing": True}; continue
    try:
        with h5py.File(p, "r") as h:
            obs = h["obs"]
            lname, lab = col(obs, LABELS)
            cname, cnt = col(obs, COUNTS)
            n = int(lab.shape[0]) if lab is not None else None
            u = cnt if cnt is not None else umis(h, n)
            labs = np.asarray(lab).astype(str)
            isc = np.isin(labs, CTRL) | np.char.startswith(np.char.lower(labs), "non") | (labs == "")
            vals, cts = np.unique(labs[~isc], return_counts=True)
            out[name] = {"label_col": lname, "count_col": cname, "cells": n, "perturbations": int(vals.size),
                         "controls": int(isc.sum()), "cells_per_pert_median": float(np.median(cts)) if cts.size else None,
                         "perts_with_ge_100_cells": int((cts >= 100).sum()),
                         "median_umi": float(np.median(np.asarray(u, float))) if u is not None else None,
                         "umi_source": "obs" if cnt is not None else "first rows"}
    except Exception as e:  # noqa: BLE001
        out[name] = {"error": f"{type(e).__name__}: {e}"[:200]}
    print(name, json.dumps(out[name]), flush=True)
(Path(__file__).resolve().parent / "work/inventory614.json").write_text(json.dumps(out, indent=1))
