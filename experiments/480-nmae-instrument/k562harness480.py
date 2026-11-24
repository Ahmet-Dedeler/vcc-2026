"""480: a K562 cell-level harness for the cancer fold whose truth is K562 (predicted from the Orion lines only).

Source: data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad (1,989,578 cells x 8,248 genes, dense float32,
chunks 3,886 x 33, gzip; cells in random order). Streams chunk-aligned row blocks (128 MB each), keeps 18,400 random
non-targeting cells and up to 150 random cells per A/B/C panel target with >= 30 cells, and writes a minimal CSR h5
(obs/target_gene, var/_index, X/{data,indices,indptr} + shape attr) that gates480.py reads.

The library of a cell is the sum over the file's 8,248 genes (98.6 % of its UMIs). Same shape as the Orion thinned
harnesses (150 cells per target, 18,400 controls); K562-GWPS median UMI ~11k vs Orion thinned 12-14k.

    uv run python k562harness480.py            # -> work/k562-gwps-panel.h5
    K562_N_CELLS=250 K562_MIN_CELLS=250 K562_OUT=work/k562-gwps-c250.h5 uv run python k562harness480.py
"""
from __future__ import annotations

import os
import time

import h5py
import numpy as np
import pandas as pd

from common480 import HERE, ROOT

SRC = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
N_CTRL, SEED = 18400, 0
N_CELLS = int(os.environ.get("K562_N_CELLS", 150))
MIN_CELLS = int(os.environ.get("K562_MIN_CELLS", 30))
OUT = HERE / os.environ.get("K562_OUT", "work/k562-gwps-panel.h5")


def main():
    t0 = time.time()
    panel = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str))
    h = h5py.File(SRC, "r")
    g = h["obs"]["gene"]
    lab = g["categories"][:].astype(str)[g["codes"][:]]
    genes = h["var"]["gene_name"][:].astype(str)
    X = h["X"]
    n, G = X.shape
    rng = np.random.default_rng(SEED)
    keep = []
    ctrl = np.flatnonzero(lab == "non-targeting")
    keep.append(np.sort(rng.choice(ctrl, N_CTRL, replace=False)))
    kept_t = []
    for t in sorted(panel):
        r = np.flatnonzero(lab == t)
        if r.size < MIN_CELLS:
            continue
        kept_t.append(t)
        keep.append(np.sort(rng.choice(r, min(N_CELLS, r.size), replace=False)))
    keep = np.sort(np.concatenate(keep))
    print(f"{len(kept_t)} panel targets, {keep.size} cells kept of {n} ({time.time() - t0:.0f}s)", flush=True)
    blk = X.chunks[0]
    with h5py.File(OUT, "w") as o:
        obs = o.create_group("obs")
        obs.create_dataset("target_gene", data=lab[keep].astype("S"))
        var = o.create_group("var")
        var.create_dataset("_index", data=genes.astype("S"))
        xg = o.create_group("X")
        xg.attrs["shape"] = (keep.size, G)
        dd = xg.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
        di = xg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
        indptr = [0]
        nnz = 0
        for lo in range(0, n, blk):
            hi = min(lo + blk, n)
            sel = keep[(keep >= lo) & (keep < hi)] - lo
            if sel.size == 0:
                continue
            D = X[lo:hi, :][sel]
            r_, c_ = np.nonzero(D)
            vals = D[r_, c_].astype(np.float32)
            dd.resize((nnz + vals.size,))
            di.resize((nnz + vals.size,))
            dd[nnz:] = vals
            di[nnz:] = c_.astype(np.int32)
            nnz += vals.size
            cnt = np.bincount(r_, minlength=sel.size)
            indptr.extend((indptr[-1] + np.cumsum(cnt)).tolist())
            if (lo // blk) % 50 == 0:
                print(f"  rows {hi}/{n}, kept nnz {nnz:,} ({time.time() - t0:.0f}s)", flush=True)
        xg.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))
        assert len(indptr) == keep.size + 1
    print(f"wrote {OUT} ({keep.size} cells, nnz {nnz:,}, {time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
