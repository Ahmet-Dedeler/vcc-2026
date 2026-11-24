"""481: per-context control statistics, streamed from h5ad (peak RSS < 1.5 GB).

For each control file (or the non-targeting rows of a harness file) this writes work/ctrl_<tag>.npz with:
  genes     the file's gene axis (str)
  totals    pooled per-gene counts over the control cells (float64)
  libs      per-cell library sizes (float64)
  ntc       per-cell NTC id when present (str; empty otherwise)
These are the only inputs the builder's context-relative constants read: the render gate (pooled CPM > min_cpm), the
scorer's DE universe (pooled CPM > 5), the call-set CPM block (top 2,200 of the gate by pooled CPM), the library classes
(q25/q50/q75 of libs) and the MDE count law (ref CPM x median lib / 1e6).

usage: uv run python ctrl_stats481.py <h5ad> <tag> [--controls-only]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad")
    ap.add_argument("tag")
    ap.add_argument("--controls-only", action="store_true")
    a = ap.parse_args()
    with h5py.File(a.h5ad, "r") as h:
        genes = str_col(h, "var/_index")
        tg = str_col(h, "obs/target_gene")
        ntc = str_col(h, "obs/ntc_id") if "ntc_id" in h["obs"] else np.array([""] * tg.size)
        rows = np.flatnonzero(tg == "non-targeting") if a.controls_only else np.arange(tg.size)
        indptr = np.asarray(h["X"]["indptr"][:], dtype=np.int64)
        data, ind = h["X"]["data"], h["X"]["indices"]
        totals = np.zeros(genes.size, dtype=np.float64)
        libs = np.zeros(rows.size, dtype=np.float64)
        # stream contiguous runs of selected rows in chunks of <= 20M stored entries
        i = 0
        while i < rows.size:
            j = i
            lo = int(indptr[rows[i]])
            while j + 1 < rows.size and rows[j + 1] == rows[j] + 1 and int(indptr[rows[j + 1] + 1]) - lo <= 20_000_000:
                j += 1
            hi = int(indptr[rows[j] + 1])
            d = np.asarray(data[lo:hi], dtype=np.float64)
            k = np.asarray(ind[lo:hi], dtype=np.int64)
            totals += np.bincount(k, weights=d, minlength=genes.size)
            sub = indptr[rows[i]: rows[j] + 2] - lo
            libs[i: j + 1] = np.add.reduceat(d, sub[:-1]) if d.size else 0.0
            empty = sub[1:] == sub[:-1]
            libs[i: j + 1][empty] = 0.0
            i = j + 1
    out = HERE / "work" / f"ctrl_{a.tag}.npz"
    np.savez_compressed(out, genes=genes, totals=totals, libs=libs, ntc=ntc[rows])
    cpm = 1e6 * totals / totals.sum()
    print(f"{a.tag}: {rows.size} control cells, {genes.size} genes; median lib {np.median(libs):.0f}; "
          f"genes > 5 CPM {int((cpm > 5).sum())}, > 10 CPM {int((cpm > 10).sum())} -> {out}")


if __name__ == "__main__":
    main()
