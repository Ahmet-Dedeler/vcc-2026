"""Pseudobulk log2-CPM fold changes from Nadig 2025 raw single-cell counts.

HepG2 and Jurkat are essential-wide CRISPRi screens (GSE264667). X is dense
float32 integer counts. Cells are streamed in chunks so the full matrix is never
resident: we keep per-perturbation sums and sums of squares, then

    lfc = log2((CPM(mean_pert) + 1) / (CPM(mean_NTC) + 1))

Sampling variance is the delta method on the per-cell mean counts, using the
empirical cell-level variance (Poisson only when a perturbation has <2 cells).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from build_replogle import save_line, sha256_file, summarise

LN2 = float(np.log(2.0))
VAR_CAP = 25.0


def _decode(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return np.asarray(out, dtype=object)


def _delta_rel(mean: np.ndarray, var_mean: np.ndarray) -> np.ndarray:
    out = np.empty(mean.shape, dtype=np.float64)
    ok = mean > 1e-6
    out[ok] = var_mean[ok] / np.square(mean[ok])
    out[~ok] = 1e6
    return out


def accumulate(path: Path, chunk: int = 8192) -> dict:
    f = h5py.File(path, "r")
    try:
        n_cells_total, n_genes = f["X"].shape
        cats = _decode(f["obs/__categories/gene"][:])
        codes = np.asarray(f["obs/gene"][:], dtype=np.int64)
        name_cats = _decode(f["var/__categories/gene_name"][:])
        gene_names = name_cats[np.asarray(f["var/gene_name"][:], dtype=np.int64)].astype(str)
        n_cat = len(cats)
        if int(codes.min()) < 0 or int(codes.max()) >= n_cat:
            raise RuntimeError(f"gene codes out of range in {path.name}")
        ctrl_ids = np.flatnonzero(cats == "non-targeting")
        if len(ctrl_ids) != 1:
            raise RuntimeError(f"expected one non-targeting category, found {ctrl_ids.tolist()}")
        ctrl_id = int(ctrl_ids[0])

        sum_all = np.zeros((n_cat, n_genes), dtype=np.float64)
        sumsq_all = np.zeros((n_cat, n_genes), dtype=np.float64)
        n_all = np.zeros(n_cat, dtype=np.int64)
        for start in range(0, n_cells_total, chunk):
            end = min(n_cells_total, start + chunk)
            block = np.asarray(f["X"][start:end], dtype=np.float64)
            c = codes[start:end]
            np.add.at(sum_all, c, block)
            np.add.at(sumsq_all, c, np.square(block))
            np.add.at(n_all, c, 1)
            done = end
            if start == 0 or done == n_cells_total or (start // chunk) % 10 == 0:
                print(f"  {path.name} cells {done}/{n_cells_total}", flush=True)
    finally:
        f.close()

    # Collapse duplicate symbols (one column pair on these files).
    uniq, inv = np.unique(gene_names, return_inverse=True)
    if len(uniq) != len(gene_names):
        sum_c = np.zeros((len(uniq), n_cat), dtype=np.float64)
        sumsq_c = np.zeros_like(sum_c)
        np.add.at(sum_c, inv, sum_all.T)
        np.add.at(sumsq_c, inv, sumsq_all.T)
        sum_all = sum_c.T
        sumsq_all = sumsq_c.T
        gene_names = uniq
        n_genes = len(uniq)

    n = n_all.astype(np.float64)
    present = n >= 1
    mean = np.zeros_like(sum_all)
    mean[present] = sum_all[present] / n[present, None]
    ss = sumsq_all - np.square(sum_all) / np.maximum(n, 1.0)[:, None]
    s2 = np.maximum(ss, 0.0) / np.maximum(n[:, None] - 1.0, 1.0)
    var_mean = s2 / np.maximum(n, 1.0)[:, None]
    poisson = mean / np.maximum(n, 1.0)[:, None]
    var_mean = np.where(n[:, None] >= 2, var_mean, poisson)

    n_ctrl = float(n[ctrl_id])
    if n_ctrl < 100:
        raise RuntimeError(f"{path.name}: only {n_ctrl} control cells")
    mu_c = mean[ctrl_id]
    vm_c = var_mean[ctrl_id]
    cpm_c = mu_c / mu_c.sum() * 1e6
    log2_c = np.log2(cpm_c + 1.0)
    rel_c = _delta_rel(mu_c, vm_c)

    pert_ids = [i for i in range(n_cat) if i != ctrl_id and n[i] >= 1]
    perts = [str(cats[i]) for i in pert_ids]
    P = len(perts)
    G = n_genes
    lfc = np.empty((P, G), dtype=np.float64)
    var = np.empty((P, G), dtype=np.float64)
    n_out = np.empty(P, dtype=np.float64)
    for pi, i in enumerate(pert_ids):
        mu = mean[i]
        cpm = mu / max(float(mu.sum()), 1e-12) * 1e6
        lfc[pi] = np.log2(cpm + 1.0) - log2_c
        rel = _delta_rel(mu, var_mean[i]) + rel_c
        var[pi] = np.clip(rel / (LN2 * LN2), 1e-8, VAR_CAP)
        n_out[pi] = n[i]

    gindex = {g: j for j, g in enumerate(gene_names.tolist())}
    own = []
    for pi, pert in enumerate(perts):
        j = gindex.get(pert)
        if j is not None:
            own.append(lfc[pi, j])
    median_own = float(np.median(own)) if own else float("nan")
    if not np.isfinite(median_own) or median_own > -0.3:
        raise RuntimeError(f"{path.name}: median own-gene log2FC is {median_own}")

    return {
        "genes": np.asarray(gene_names, dtype=str),
        "perts": np.asarray(perts, dtype=str),
        "lfc": lfc,
        "var": var,
        "n_cells": n_out,
        "n_guides": np.ones(P, dtype=np.int32),
        "ctrl_log": np.log(cpm_c + 1.0),
        "n_ctrl": n_ctrl,
        "median_own_gene_lfc": median_own,
        "n_cells_total": int(n_cells_total),
        "min_n_cells": float(n_out.min()),
        "median_n_cells_targets": float(np.median(n_out)),
    }


def main() -> None:
    root = Path(__file__).resolve().parent
    data = root / "work" / "data"
    out = root / "work" / "responses"
    jobs = [
        ("hepg2", data / "GSE264667_hepg2_raw_singlecell_01.h5ad"),
        ("jurkat", data / "GSE264667_jurkat_raw_singlecell_01.h5ad"),
    ]
    for name, path in jobs:
        print(f"building {name}", flush=True)
        digest = sha256_file(path)
        rec = accumulate(path)
        meta = {
            "line_or_screen": name,
            "sha256": digest,
            "nbytes": path.stat().st_size,
            "source_file": path.name,
            "n_cells_total": rec["n_cells_total"],
            "min_n_cells": rec["min_n_cells"],
            **summarise(rec),
        }
        save_line(out / f"{name}.npz", rec, meta)
        print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
