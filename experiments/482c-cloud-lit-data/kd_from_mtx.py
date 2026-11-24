"""On-target knockdown from a MatrixMarket counts file plus cell/gene annotation CSVs.

Used for GSE345058 (A549 1,000-gene Cas9 CROP-seq). Assumes the mtx axes match the
annotation row order. Writes a JSON summary next to the annotations.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def load_genes(path: Path) -> list[str]:
    with open_text(path) as f:
        rows = list(csv.DictReader(f))
    # gene_annotations.csv: index, ensembl_id, gene_symbol, ...
    if "gene_symbol" in rows[0]:
        return [r["gene_symbol"] for r in rows]
    return [r[next(iter(r))] for r in rows]


def load_perts(path: Path) -> list[str]:
    with open_text(path) as f:
        rows = csv.DictReader(f)
        return [r["Perturbation"] for r in rows]


def mtx_header(path: Path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            nrows, ncols, nnz = map(int, line.split())
            return nrows, ncols, nnz
    raise RuntimeError("no mtx header")


def main():
    mtx = Path(sys.argv[1])
    cells_csv = Path(sys.argv[2])
    genes_csv = Path(sys.argv[3])
    out = Path(sys.argv[4])

    genes = load_genes(genes_csv)
    perts = load_perts(cells_csv)
    nrows, ncols, nnz = mtx_header(mtx)
    print(f"mtx {nrows} x {ncols} nnz {nnz}; genes {len(genes)}; cells {len(perts)}", flush=True)

    if (nrows, ncols) == (len(genes), len(perts)):
        orientation = "genes_x_cells"
    elif (nrows, ncols) == (len(perts), len(genes)):
        orientation = "cells_x_genes"
    else:
        raise SystemExit(f"dimension mismatch mtx {(nrows, ncols)} vs genes {len(genes)} cells {len(perts)}")

    # perturbation codes: 0 = nontargeting, 1..K = targets in first-seen order
    name_to_code: dict[str, int] = {}
    codes = np.empty(len(perts), dtype=np.int32)
    for i, p in enumerate(perts):
        key = p.strip()
        if key.lower() in {"nontargeting", "non-targeting", "ntc"}:
            codes[i] = 0
            continue
        if key not in name_to_code:
            name_to_code[key] = len(name_to_code) + 1
        codes[i] = name_to_code[key]
    K = len(name_to_code)
    code_to_name = {v: k for k, v in name_to_code.items()}
    n_cells = np.bincount(codes, minlength=K + 1)

    # gene row/col index for each target that exists in the matrix
    gene_index = {}
    for i, g in enumerate(genes):
        gene_index.setdefault(g, i)  # first occurrence

    # sums[code] of the target's own gene; ntc_sums[code] of that gene in NTC cells
    self_sum = np.zeros(K + 1, dtype=np.float64)
    ntc_sum = np.zeros(K + 1, dtype=np.float64)
    # map matrix gene index -> pert code (0 if not a target)
    if orientation == "genes_x_cells":
        gmax = nrows
    else:
        gmax = ncols
    gene_code = np.zeros(gmax, dtype=np.int32)
    missing = []
    for name, code in name_to_code.items():
        gi = gene_index.get(name)
        if gi is None or gi >= gmax:
            missing.append(name)
            continue
        gene_code[gi] = code

    print(f"targets {K} missing from gene axis {len(missing)} ntc cells {int(n_cells[0])}", flush=True)

    seen = 0
    with gzip.open(mtx, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            # skip the dimension line (three ints, no float)
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                a = int(parts[0])
                b = int(parts[1])
                v = float(parts[2])
            except ValueError:
                continue
            # MatrixMarket is 1-indexed
            a -= 1
            b -= 1
            if orientation == "genes_x_cells":
                gi, ci = a, b
            else:
                ci, gi = a, b
            if ci < 0 or ci >= len(codes) or gi < 0 or gi >= gmax:
                continue
            gc = gene_code[gi]
            if gc == 0:
                continue
            pc = codes[ci]
            if pc == 0:
                ntc_sum[gc] += v
            elif pc == gc:
                self_sum[gc] += v
            seen += 1
            if seen % 50_000_000 == 0:
                print(f"  target-gene entries {seen}", flush=True)

    rows = []
    for code in range(1, K + 1):
        nt = n_cells[code]
        nc = n_cells[0]
        if nt == 0 or nc == 0:
            continue
        ms = self_sum[code] / nt
        mn = ntc_sum[code] / nc
        rows.append((code_to_name[code], int(nt), ms, mn))

    expressed = [r for r in rows if r[3] > 0]
    down = [r for r in expressed if r[2] < r[3]]
    half = [r for r in expressed if r[2] <= 0.5 * r[3]]
    ratios = np.array([r[2] / r[3] for r in expressed]) if expressed else np.array([])

    summary = {
        "orientation": orientation,
        "mtx_shape": [nrows, ncols],
        "nnz": nnz,
        "n_cells": len(perts),
        "n_genes_axis": len(genes),
        "n_targets": K,
        "n_ntc_cells": int(n_cells[0]),
        "targets_missing_from_gene_axis": missing,
        "n_targets_with_ntc_mean_gt_0": len(expressed),
        "frac_mean_lt_ntc": (len(down) / len(expressed)) if expressed else None,
        "frac_mean_le_half_ntc": (len(half) / len(expressed)) if expressed else None,
        "median_residual_fraction": float(np.median(ratios)) if len(ratios) else None,
        "p10_residual_fraction": float(np.quantile(ratios, 0.1)) if len(ratios) else None,
        "p90_residual_fraction": float(np.quantile(ratios, 0.9)) if len(ratios) else None,
        "median_cells_per_target": int(np.median(n_cells[1:])),
        "min_cells_per_target": int(n_cells[1:].min()) if K else None,
        "max_cells_per_target": int(n_cells[1:].max()) if K else None,
    }
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
