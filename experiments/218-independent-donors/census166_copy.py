"""Metadata-only census of every perturbation h5ad on disk.

Rule 1 of the brief: get the target list and the gene axis from metadata, compute the panel
overlap, and only touch the big object if the overlap clears the bar. Everything here reads
`.obs` categorical codes and `.var` index through h5py -- no `.X` is ever opened, so a 8.8 GB
file costs the same as a 20 MB one.

Writes `work/census.json`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = set(pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"])
GENES = set(pd.read_csv(REPO / "data/challenge/gene_names.csv")["gene_name"])
assert len(PANEL) == 300 and len(GENES) == 18533, (len(PANEL), len(GENES))

# Column names used across scPerturb / Nadig / GEO conventions, in preference order.
PERT_COLS = ("perturbation", "gene", "target_gene", "gene_target", "guide_identity",
             "gene_symbol", "perturbation_name", "grna_target", "target")
LINE_COLS = ("cell_line", "cell_type", "celltype", "line", "tissue_type", "organism")


def _decode(a):
    if a.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return np.asarray(a)


def read_col(h, group, name):
    """Read one .obs column, resolving anndata categorical encodings."""
    g = h[group]
    if name not in g:
        return None
    node = g[name]
    if isinstance(node, h5py.Group):  # categorical
        cats = _decode(node["categories"][:])
        codes = node["codes"][:]
        return cats[np.clip(codes, 0, len(cats) - 1)]
    arr = node[:]
    if arr.dtype.kind in "iu" and "__categories" in g:
        # legacy anndata categorical
        cats = _decode(g["__categories"][name][:])
        return cats[np.clip(arr, 0, len(cats) - 1)]
    return _decode(arr)


def var_names(h):
    v = h["var"]
    key = v.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    if key in v:
        return _decode(v[key][:])
    return None


def inspect(path: Path) -> dict:
    rec: dict = {"file": str(path.relative_to(REPO)),
                 "size_gb": round(path.stat().st_size / 1e9, 3)}
    with h5py.File(path, "r") as h:
        obs = h["obs"]
        cols = list(obs.keys())
        rec["obs_columns"] = [c for c in cols if not c.startswith("_")]
        idxkey = obs.attrs.get("_index", "_index")
        if isinstance(idxkey, bytes):
            idxkey = idxkey.decode()
        try:
            rec["n_cells"] = int(obs[idxkey].shape[0])
        except Exception:
            rec["n_cells"] = None

        gv = var_names(h)
        if gv is not None:
            rec["n_genes"] = int(gv.size)
            rec["gene_axis_overlap_18533"] = int(len(GENES & set(gv.tolist())))

        pcol = next((c for c in PERT_COLS if c in cols), None)
        rec["pert_column"] = pcol
        if pcol is None:
            return rec
        p = read_col(h, "obs", pcol)
        # split combinatorial labels on the usual separators, keep atoms
        atoms: set[str] = set()
        for s in set(p.tolist()):
            for a in str(s).replace("+", "_").replace(",", "_").split("_"):
                a = a.strip()
                if a:
                    atoms.add(a)
        uniq = set(p.tolist())
        rec["n_pert_labels"] = len(uniq)
        rec["panel_overlap"] = int(len(PANEL & uniq))
        rec["panel_overlap_atoms"] = int(len(PANEL & atoms))
        vals, cnts = np.unique(p, return_counts=True)
        hit = np.isin(vals, list(PANEL))
        rec["panel_cells"] = int(cnts[hit].sum())
        rec["median_cells_per_label"] = int(np.median(cnts))
        rec["example_labels"] = sorted(uniq)[:8]

        lcol = next((c for c in LINE_COLS if c in cols), None)
        rec["line_column"] = lcol
        if lcol is not None:
            lv = read_col(h, "obs", lcol)
            u, c = np.unique(lv, return_counts=True)
            rec["lines"] = {str(a): int(b) for a, b in
                            sorted(zip(u, c), key=lambda t: -t[1])[:12]}
        for extra in ("perturbation_type", "nperts", "disease", "tissue", "organism"):
            if extra in cols:
                v = read_col(h, "obs", extra)
                u, c = np.unique(v, return_counts=True)
                rec[extra] = {str(a): int(b) for a, b in
                              sorted(zip(u, c), key=lambda t: -t[1])[:6]}
    return rec


def main(paths: list[str]) -> None:
    files: list[Path] = []
    for p in paths:
        q = Path(p)
        if not q.is_absolute():
            q = REPO / q
        files.extend(sorted(q.rglob("*.h5ad")) if q.is_dir() else [q])

    out = []
    for f in files:
        try:
            rec = inspect(f)
        except Exception as e:  # a broken or non-anndata file must not stop the census
            rec = {"file": str(f.relative_to(REPO)), "error": f"{type(e).__name__}: {e}"}
        out.append(rec)
        print(f"{rec.get('panel_overlap', '-'):>4}  {rec.get('size_gb', 0):>7.2f}GB  "
              f"{rec['file']}  perts={rec.get('n_pert_labels', '-')} "
              f"cells={rec.get('n_cells', '-')} lines={list((rec.get('lines') or {}))[:3]}",
              flush=True)

    (OUT / "census.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {OUT / 'census.json'} ({len(out)} files)")


if __name__ == "__main__":
    main(sys.argv[1:] or ["data/external"])
