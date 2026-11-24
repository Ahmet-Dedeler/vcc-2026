"""217 falsifier 6 — stored-entry density against the 4.75e9 cap, and counts/cell against 1e6.

`vcc prep --dry-run` proves the format contract but does not print the density, so 210 read it
off the file.  This is that read, made a script so the known-answer row is reproducible: run it
on `submissions/210-cd4f034.h5ad` and it must print 210's published **1,577,729,493 nnz,
33.2154% of the cap, 4,382.6 nnz/cell** and the per-context stored entries
538,693,542 / 512,429,392 / 526,606,559.

Reads only HDF5 dataset shapes and the `indptr`/`codes` vectors, so peak memory is a few MB
regardless of the artifact's size.

    uv run python nnz217.py <pred.h5ad> [<pred.h5ad> ...]
"""
from __future__ import annotations

import sys

import h5py
import numpy as np

CAP = 4.75e9
CELL_CAP = 1_000_000


def main() -> None:
    for path in sys.argv[1:]:
        with h5py.File(path, "r") as handle:
            nnz = int(handle["X"]["data"].shape[0])
            shape = tuple(int(v) for v in handle["X"].attrs["shape"])
            group = handle["obs"]["context"]
            cats = np.asarray(group["categories"][:]).astype(str)
            codes = np.asarray(group["codes"][:])
            indptr = np.asarray(handle["X"]["indptr"][:])
            data_dtype = handle["X"]["data"].dtype
        per_cell = np.diff(indptr)
        print(f"{path}   dtype {data_dtype}")
        print(f"  nnz {nnz:,}   dense would be {shape[0] * shape[1]:,}   "
              f"% of the 4.75e9 cap {100 * nnz / CAP:.4f}%   nnz/cell {nnz / shape[0]:.1f}   "
              f"{'PASS' if nnz <= CAP else 'FAIL'}")
        for i, ctx in enumerate(cats):
            sel = codes == i
            print(f"    ctx {ctx}: cells {int(sel.sum()):,}   "
                  f"stored entries {int(per_cell[sel].sum()):,}")
        print(f"  max stored entries in one cell {int(per_cell.max()):,} "
              f"(the 1e6 cap is on COUNTS, measured in profile217.py)")


if __name__ == "__main__":
    main()
