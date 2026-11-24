#!/usr/bin/env python3
"""618 stage B: the h1-trainval harness truth (CSR h5 in 614's format, which ref614.py reads).

  perturbed: every h1-trainval cell of its 126 targets (400 each; disjoint from h1-test's 85);
  controls : h1-test's own 18,400 control cells, the reference control set of the validated H1 surface (528's truth528
             precedent: trainval's perturbations against h1-test's controls). Both files are the same 2025 Arc H1 Flex
             screen on the same 18,080-gene axis (asserted).
Streams both h5ads in row blocks; writes work/tv-truth.h5 (obs/target_gene, var/_index, X/{data float32, indices int32,
indptr int64} + shape) and work/tv-truth.json.
    uv run python truth618.py
"""
from __future__ import annotations

import json
import os
import sys

import anndata as ad
import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
TEST = os.path.join(REPO, "data/harness/h1-test.h5ad")
TV = os.path.join(REPO, "data/harness/h1-trainval.h5ad")
CTRL = "non-targeting"


def meta(path):
    a = ad.read_h5ad(path, backed="r")
    lab = a.obs["target_gene"].astype(str).to_numpy()
    genes = np.asarray(a.var_names.astype(str))
    a.file.close()
    return lab, genes


def rows_iter(path, rows, block=4000):
    """Yield (row_index, data, indices) for the given sorted rows, reading contiguous blocks."""
    with h5py.File(path, "r") as h:
        ip = h["X"]["indptr"][:]
        dd, di = h["X"]["data"], h["X"]["indices"]
        for s0 in range(0, len(rows), block):
            rs = rows[s0:s0 + block]
            r0, r1 = int(rs[0]), int(rs[-1]) + 1
            a0, a1 = int(ip[r0]), int(ip[r1])
            bd, bi = dd[a0:a1], di[a0:a1]
            for r in rs:
                s, e = int(ip[r]) - a0, int(ip[r + 1]) - a0
                yield int(r), bd[s:e], bi[s:e]


def main() -> None:
    lt, gt = meta(TEST)
    lv, gv = meta(TV)
    assert np.array_equal(gt, gv), "h1-test and h1-trainval gene axes differ"
    ctrl_rows = np.flatnonzero(lt == CTRL)
    assert len(ctrl_rows) == 18400, len(ctrl_rows)
    pert_rows = np.flatnonzero(lv != CTRL)
    targets = sorted(set(lv[pert_rows]))
    assert len(targets) == 126 and not (set(targets) & set(lt)), "trainval targets must be 126 and disjoint from h1-test"
    # perturbed rows grouped by target (sorted), file order inside a target; then the controls
    order_v = np.concatenate([pert_rows[lv[pert_rows] == t] for t in targets])
    labels = np.concatenate([lv[order_v], np.array([CTRL] * len(ctrl_rows))])
    out = os.path.join(HERE, "work", "tv-truth.h5")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with h5py.File(out + ".part", "w") as o:
        od = o.create_dataset("X/data", (0,), maxshape=(None,), dtype="f4", chunks=(1 << 20,))
        oi = o.create_dataset("X/indices", (0,), maxshape=(None,), dtype="i4", chunks=(1 << 20,))
        indptr, lib = [0], []

        def put(bufd, bufi):
            n0 = od.shape[0]
            d = np.concatenate(bufd) if bufd else np.zeros(0, "f4")
            i = np.concatenate(bufi) if bufi else np.zeros(0, "i4")
            od.resize((n0 + len(d),))
            oi.resize((n0 + len(i),))
            od[n0:] = d
            oi[n0:] = i

        for path, rows, order in ((TV, np.sort(order_v), order_v), (TEST, ctrl_rows, ctrl_rows)):
            # read sorted, emit in `order` (for TV the target-grouped order; for TEST the file order)
            got = {}
            for r, d, i in rows_iter(path, rows):
                got[r] = (d.astype(np.float32), i.astype(np.int32))
            bufd, bufi = [], []
            for r in order:
                d, i = got.pop(int(r))
                bufd.append(d)
                bufi.append(i)
                indptr.append(indptr[-1] + len(d))
                lib.append(float(d.sum()))
                if len(bufd) >= 4000:
                    put(bufd, bufi)
                    bufd, bufi = [], []
            put(bufd, bufi)
            del got
        o.create_dataset("X/indptr", data=np.asarray(indptr, dtype=np.int64))
        o["X"].attrs["shape"] = np.array([len(labels), len(gt)])
        o.create_dataset("obs/target_gene", data=labels.astype("S"))
        o.create_dataset("var/_index", data=gt.astype("S"))
    os.replace(out + ".part", out)
    lib = np.asarray(lib)
    isc = labels == CTRL
    info = dict(targets=len(targets), cells=int(len(labels)), perturbed=int((~isc).sum()), controls=int(isc.sum()),
                genes=int(len(gt)), nnz=int(indptr[-1]), library_median_pert=float(np.median(lib[~isc])),
                library_median_ctrl=float(np.median(lib[isc])),
                source={"perturbed": "data/harness/h1-trainval.h5ad", "controls": "data/harness/h1-test.h5ad"})
    json.dump(info, open(os.path.join(HERE, "work", "tv-truth.json"), "w"), indent=1)
    (open(os.path.join(HERE, "work", "targets_tv.txt"), "w")).write("\n".join(targets) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    sys.exit(main())
