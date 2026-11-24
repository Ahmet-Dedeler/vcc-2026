"""196 — stream the SHIPPED per-(context, target, gene) panel totals off the a15 artifact.

`submissions/191-a15_local.h5ad` is byte-identical to the shipped a15 package (191 §6), so
these are a15's own delivered counts, not a proxy.  Also records, per block, the per-cell
total distribution and the stored-entry count, which are the two format-contract quantities
the depth sweep has to respect (AGENTS.md §2: 4.75e9 stored entries, 1e6 counts/cell).

    uv run python totals196.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
PRED = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "submissions/191-a15_local.h5ad"
TAG = sys.argv[2] if len(sys.argv) > 2 else "a15"
OUT = HERE / "work" / f"chal-totals-{TAG}.npz"


def main() -> None:
    t0 = time.time()
    A = ad.read_h5ad(PRED, backed="r")
    tgt = np.asarray(A.obs["target_gene"].astype(str).values)
    ctx = np.asarray(A.obs["context"].astype(str).values)
    genes = np.asarray([str(g) for g in A.var_names])
    A.file.close()
    n_cells, G = tgt.size, genes.size
    key = np.char.add(np.char.add(ctx, "|"), tgt)
    blocks, first = np.unique(key, return_index=True)
    order = np.argsort(first)
    blocks = blocks[order]
    starts = first[order]
    # every block must be 400 CONTIGUOUS rows, else the slice read below is wrong
    ends = np.append(starts[1:], n_cells)
    assert np.all(ends - starts == 400), "blocks are not contiguous 400-row runs"
    for b, s, e in zip(blocks, starts, ends):
        assert np.all(key[s:e] == b)
    B = blocks.size
    print(f"{PRED.name}: {n_cells} cells, {G} genes, {B} blocks [{time.time()-t0:.0f}s]",
          flush=True)

    totals = np.zeros((B, G), dtype=np.int64)
    cellstat = np.zeros((B, 4), dtype=np.float64)   # min, max, mean, sd of per-cell totals
    nnz_block = np.zeros(B, dtype=np.int64)
    with h5py.File(PRED, "r") as f:
        indptr = f["X"]["indptr"][:]
        dset_d, dset_i = f["X"]["data"], f["X"]["indices"]
        for k in range(B):
            s, e = int(starts[k]), int(ends[k])
            a, b = int(indptr[s]), int(indptr[e])
            data = dset_d[a:b]
            idx = dset_i[a:b]
            di = data.astype(np.int64)
            assert np.all(data == di) and np.all(di >= 0)
            np.add.at(totals[k], idx, di)
            counts = np.diff(indptr[s:e + 1])
            rowsum = np.add.reduceat(di, np.cumsum(np.append(0, counts))[:-1])
            rowsum = np.where(counts > 0, rowsum, 0)
            cellstat[k] = (rowsum.min(), rowsum.max(), rowsum.mean(), rowsum.std())
            nnz_block[k] = b - a
            if (k + 1) % 100 == 0:
                print(f"  block {k+1}/{B} [{time.time()-t0:.0f}s]", flush=True)
    np.savez_compressed(OUT, blocks=blocks, genes=genes, totals=totals,
                        cellstat=cellstat, nnz=nnz_block)
    print(f"wrote {OUT}  nnz {nnz_block.sum()}  "
          f"per-cell median {np.median(cellstat[:, 2]):.1f}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
