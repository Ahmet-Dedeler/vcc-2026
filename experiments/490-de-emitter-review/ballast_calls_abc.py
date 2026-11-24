#!/usr/bin/env python3
"""490 F-B2: the ballast cells' footprint on the DELIVERED CALLS of real A/B/C blocks (453dr_final.h5ad, context D = A).

The emitter designs |z| on the n_car carriers (prescribed set up, complement capped at z 1.90). The scorer runs the
Wilcoxon on all 400 cells. Per sampled block: BH calls with the carriers alone (the emitter's frame, `rc_emit.de_from_block`)
vs with all 400 cells (the scorer's frame), both against the context's real controls on the scorer's DE universe.
Reports leaks (called at 400, not on carriers), losses, and the call-count change. Peak RSS ~3 GB (the sorted control
block is n_keep x 18,400 float64); run only with vm_stat free+inactive >= 33 GiB before 04:30Z.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026 import rc_emit  # noqa: E402

PRED = ROOT / "experiments/453-abc-dress-rehearsal/work/run/453dr_final.h5ad"
CTRL = {"D": "A", "E": "B", "F": "C"}


def categorical(h, name):
    grp = h[f"obs/{name}"]
    cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in grp["categories"][:]])
    return cats[grp["codes"][:]]


def main():
    c = sys.argv[1] if len(sys.argv) > 1 else "D"
    step = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    import anndata as ad
    ctl = ad.read_h5ad(ROOT / f"data/challenge/context_{CTRL[c]}.h5ad")
    X = ctl.X.tocsr() if not sp.isspmatrix_csr(ctl.X) else ctl.X
    control = rc_emit.control_block_from_counts(X, min_cpm=5.0)
    del ctl, X
    recs = []
    with h5py.File(PRED, "r") as h:
        ctx = categorical(h, "context")
        indptr = h["X/indptr"][:]
        g = int(h["X"].attrs["shape"][1])
        rows_c = np.flatnonzero(ctx == c)
        for r0 in rows_c[::400][::step]:
            a, b = int(indptr[r0]), int(indptr[r0 + 400])
            blk = sp.csr_matrix((h["X/data"][a:b].astype(np.float64), h["X/indices"][a:b],
                                 indptr[r0:r0 + 401] - a), shape=(400, g)).toarray()
            bal = blk.sum(1) > 9e5
            _, _, q_car, _ = rc_emit.de_from_block(blk[~bal], control)
            _, _, q_all, _ = rc_emit.de_from_block(blk, control)
            c_car, c_all = q_car < 0.05, q_all < 0.05
            recs.append(dict(K=int(bal.sum()), n_car_calls=int(c_car.sum()), n_all_calls=int(c_all.sum()),
                             leaks=int((c_all & ~c_car).sum()), losses=int((c_car & ~c_all).sum())))
            print(recs[-1], flush=True)
    K = np.array([x["K"] for x in recs])
    d = np.array([x["n_all_calls"] - x["n_car_calls"] for x in recs])
    lk = np.array([x["leaks"] for x in recs])
    ls = np.array([x["losses"] for x in recs])
    summ = dict(context=c, blocks=len(recs), K_q=[float(v) for v in np.quantile(K, [0, .5, 1])],
                dcalls_q=[float(v) for v in np.quantile(d, [0, .1, .5, .9, 1])],
                leaks_q=[float(v) for v in np.quantile(lk, [0, .5, .9, 1])],
                losses_q=[float(v) for v in np.quantile(ls, [0, .5, .9, 1])],
                n_car_calls_median=float(np.median([x["n_car_calls"] for x in recs])),
                corr_K_leaks=float(np.corrcoef(K, lk)[0, 1]) if K.std() > 0 and lk.std() > 0 else None)
    print(json.dumps(summ, indent=1))
    (HERE / "work" / f"ballast_calls_{c}.json").write_text(json.dumps(dict(summary=summ, blocks=recs), indent=1))


if __name__ == "__main__":
    main()
