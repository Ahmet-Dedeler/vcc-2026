"""480: a power-matched H1 control surface. h1-test at the cancer folds' power: 150 random cells per target (of 400),
all 18,400 controls, and every cell's counts binomially thinned by f (default 0.6: 20.2k -> ~12k median UMI, the
folds' 11.3-14.3k). Same line, same donors, same delivered predictions; only the truth's power changes.

    uv run python h1thin480.py [--cells 150] [--f 0.6] [--seed 0]     # -> work/h1-thin.h5 (then gates480.py)
"""
from __future__ import annotations

import argparse
import time

import h5py
import numpy as np

from common480 import HERE, ROOT
from gates480 import read_str

SRC = ROOT / "data/harness/h1-test.h5ad"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=150)
    ap.add_argument("--f", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=str(HERE / "work/h1-thin.h5"))
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(a.seed)
    h = h5py.File(SRC, "r")
    lab = read_str(h["obs"]["target_gene"])
    genes = read_str(h["var"]["_index"])
    X = h["X"]
    shape = tuple(int(v) for v in X.attrs["shape"])
    indptr = X["indptr"][:].astype(np.int64)
    ctrl = np.flatnonzero(lab == "non-targeting")
    keep = [ctrl]
    for t in sorted(set(lab) - {"non-targeting"}):
        r = np.flatnonzero(lab == t)
        keep.append(np.sort(rng.choice(r, min(a.cells, r.size), replace=False)))
    keep = np.sort(np.concatenate(keep))
    print(f"{keep.size} cells kept of {shape[0]} ({len(keep) - ctrl.size} perturbed), thinning f={a.f}", flush=True)
    with h5py.File(a.out, "w") as o:
        o.create_group("obs").create_dataset("target_gene", data=lab[keep].astype("S"))
        o.create_group("var").create_dataset("_index", data=genes.astype("S"))
        xg = o.create_group("X")
        xg.attrs["shape"] = (keep.size, shape[1])
        dd = xg.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
        di = xg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
        new_ptr = [0]
        nnz = 0
        BLK = 4096
        keep_set = np.zeros(shape[0], bool)
        keep_set[keep] = True
        for lo in range(0, shape[0], BLK):
            hi = min(lo + BLK, shape[0])
            rows = np.flatnonzero(keep_set[lo:hi]) + lo
            if rows.size == 0:
                continue
            b0, b1 = int(indptr[lo]), int(indptr[hi])
            V = X["data"][b0:b1]
            C = X["indices"][b0:b1]
            vals, idx, cnts = [], [], []
            for r in rows:
                a0, a1 = int(indptr[r]) - b0, int(indptr[r + 1]) - b0
                v = rng.binomial(np.round(V[a0:a1]).astype(np.int64), a.f)
                m = v > 0
                vals.append(v[m].astype(np.float32))
                idx.append(C[a0:a1][m].astype(np.int32))
                cnts.append(int(m.sum()))
            vv, ii = np.concatenate(vals), np.concatenate(idx)
            dd.resize((nnz + vv.size,))
            di.resize((nnz + vv.size,))
            dd[nnz:] = vv
            di[nnz:] = ii
            nnz += vv.size
            new_ptr.extend((new_ptr[-1] + np.cumsum(cnts)).tolist())
        print(f"  nnz {nnz:,} ({time.time() - t0:.0f}s)", flush=True)
        assert len(new_ptr) == keep.size + 1
        xg.create_dataset("indptr", data=np.asarray(new_ptr, dtype=np.int64))
    print(f"wrote {a.out} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
