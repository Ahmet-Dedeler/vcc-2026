#!/usr/bin/env python3
"""490: is the CLEAN4 block prior actually what the scorer reads on REAL A/B/C blocks (453dr_final.h5ad)?

Without the design key the head cannot be separated exactly, so the check runs on the 2,200 highest-CPM genes of the
DE universe (the CPM block is the top 2,200 by control CPM AMONG NON-HEAD genes, so these are block genes except for
the few head genes with high CPM). Expected: non-protected block genes (80 %) carry the prior sign exactly; protected
ones (20 %) carry E0's sign, which agrees with the prior ~half the time; so agreement ~0.9 per block, lower where the
prior's dose lambda < 1 fails to flip. lfc over all 400 cells (the scorer's frame).
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
sys.path.insert(0, str(HERE))
from ballast_abc import CTRL, PRED, categorical, control_stats  # noqa: E402


def main():
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    z = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
    psign = z["psign"].astype(float)
    out = {}
    with h5py.File(PRED, "r") as h:
        genes = np.asarray(h["var/_index"][:]).astype(str)
        assert np.array_equal(genes, z["genes"].astype(str))
        ctx = categorical(h, "context")
        tg = categorical(h, "target_gene")
        indptr = h["X/indptr"][:]
        g = genes.size
        for c in ("D", "E", "F"):
            ref_lin, _ = control_stats(ROOT / f"data/challenge/context_{CTRL[c]}.h5ad")
            keep = ref_lin > 5.0
            gi = np.flatnonzero(keep)
            top = gi[np.argsort(-ref_lin[gi], kind="stable")[:2200]]
            rows_c = np.flatnonzero(ctx == c)
            ag = []
            for r0 in rows_c[::400][::step]:
                a, b = int(indptr[r0]), int(indptr[r0 + 400])
                blk = sp.csr_matrix((h["X/data"][a:b].astype(np.float64), h["X/indices"][a:b],
                                     indptr[r0:r0 + 401] - a), shape=(400, g)).toarray()
                cpm = blk * (1e6 / blk.sum(1)[:, None])
                lfc = np.log2((cpm.mean(0)[top] + 1e-9) / (ref_lin[top] + 1e-9))
                tt = top[genes[top] != tg[r0]]
                lf = np.log2((cpm.mean(0)[tt] + 1e-9) / (ref_lin[tt] + 1e-9))
                ag.append(float((np.sign(lf) == psign[tt]).mean()))
            ag = np.array(ag)
            out[c] = dict(blocks=len(ag), agreement_q=[round(float(v), 4) for v in np.quantile(ag, [0, .1, .5, .9, 1])],
                          mean=round(float(ag.mean()), 4))
            print(c, out[c], flush=True)
    (HERE / "work" / "prior_delivery_abc.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
