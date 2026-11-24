#!/usr/bin/env python3
"""490 F-B on REAL A/B/C blocks: the ballast cells' footprint on the DE table's lfc (no truth; footprint only).

453dr_final.h5ad is the October-path rehearsal on A/B/C relabelled D/E/F, with the DE configuration 446b ships
(clean4q80, coarse MDE, zero repair -5, ratio-mode ballast; only the bulk load and names differ). For a sample of blocks
per context: carriers = rows with library < 9e5, ballast = the rest. lfc over the carriers (the emitter's own frame)
vs lfc over all 400 cells (the scorer's), on the scorer's DE universe (mean per-cell CPM > 5 on the context's controls).
Reads one 400-row block at a time with h5py; controls are streamed in row chunks. Peak RSS < 2 GB.
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
PRED = ROOT / "experiments/453-abc-dress-rehearsal/work/run/453dr_final.h5ad"
CTRL = {"D": "A", "E": "B", "F": "C"}          # 453's relabelled bundle; checked below against the logged classes
LOGGED_CLASSES = {"D": [13986.0, 20109.0, 27234.5], "E": [15721.5, 19946.0, 24708.8], "F": [14414.0, 20034.0, 26741.5]}


def control_stats(path):
    with h5py.File(path, "r") as h:
        indptr = h["X/indptr"][:]
        n, g = int(indptr.size - 1), int(h["X"].attrs["shape"][1])
        cpm_sum = np.zeros(g)
        libs = np.zeros(n)
        for lo in range(0, n, 2000):
            hi = min(lo + 2000, n)
            a, b = int(indptr[lo]), int(indptr[hi])
            m = sp.csr_matrix((h["X/data"][a:b].astype(np.float64), h["X/indices"][a:b], indptr[lo:hi + 1] - a),
                              shape=(hi - lo, g))
            s = np.asarray(m.sum(1)).ravel()
            libs[lo:hi] = s
            cpm_sum += np.asarray(m.multiply(1e6 / np.where(s > 0, s, 1.0)[:, None]).sum(0)).ravel()
    return cpm_sum / n, libs


def categorical(h, name):
    grp = h[f"obs/{name}"]
    cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in grp["categories"][:]])
    return cats[grp["codes"][:]]


def main():
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    out = {}
    with h5py.File(PRED, "r") as h:
        ctx = categorical(h, "context")
        tg = categorical(h, "target_gene")
        indptr = h["X/indptr"][:]
        g = int(h["X"].attrs["shape"][1])
        for c in ("D", "E", "F"):
            ref_lin, libs = control_stats(ROOT / f"data/challenge/context_{CTRL[c]}.h5ad")
            q = [round(float(v), 1) for v in np.quantile(libs, [0.25, 0.5, 0.75])]
            assert np.allclose(q, LOGGED_CLASSES[c], atol=0.6), (c, q)
            keep = ref_lin > 5.0
            rows_c = np.flatnonzero(ctx == c)
            starts = rows_c[::400]
            recs = []
            for bi, r0 in enumerate(starts[::step]):
                a, b = int(indptr[r0]), int(indptr[r0 + 400])
                blk = sp.csr_matrix((h["X/data"][a:b].astype(np.float64), h["X/indices"][a:b],
                                     indptr[r0:r0 + 401] - a), shape=(400, g)).toarray()
                lib = blk.sum(1)
                bal = lib > 9e5
                cpm = blk * (1e6 / lib[:, None])
                m_car = cpm[~bal].mean(0)[keep]
                m_all = cpm.mean(0)[keep]
                l_car = np.log2((m_car + 1e-9) / (ref_lin[keep] + 1e-9))
                l_all = np.log2((m_all + 1e-9) / (ref_lin[keep] + 1e-9))
                own = np.asarray(h["var/_index"][:]).astype(str)[keep] == tg[r0] if bi == 0 else None
                d = l_all - l_car
                top = np.zeros(int(keep.sum()), bool)
                top[np.argsort(-ref_lin[keep], kind="stable")[:2600]] = True
                recs.append(dict(K=int(bal.sum()), med_abs=float(np.median(np.abs(d))),
                                 q99_abs=float(np.quantile(np.abs(d), 0.99)),
                                 sign_flip_all=int((np.sign(l_all) != np.sign(l_car)).sum()),
                                 sign_flip_top2600=int(((np.sign(l_all) != np.sign(l_car)) & top).sum()),
                                 n_keep=int(keep.sum())))
            K = np.array([x["K"] for x in recs])
            flips = np.array([x["sign_flip_all"] for x in recs])
            ft = np.array([x["sign_flip_top2600"] for x in recs])
            out[c] = dict(blocks=len(recs), K_q=[float(v) for v in np.quantile(K, [0, .5, 1])],
                          median_abs_dlfc=float(np.median([x["med_abs"] for x in recs])),
                          q99_abs_dlfc_median=float(np.median([x["q99_abs"] for x in recs])),
                          sign_flips_per_block_mean=float(flips.mean()), sign_flips_max=int(flips.max()),
                          sign_flips_top2600_mean=float(ft.mean()), n_keep=recs[0]["n_keep"],
                          corr_K_flips=float(np.corrcoef(K, flips)[0, 1]) if K.std() > 0 else None)
            print(c, out[c], flush=True)
    (HERE / "work" / "ballast_abc.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
