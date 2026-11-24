"""End-to-end amplitude curve on h1-test: re-render, RE-EMIT, and score with cell_eval2.

Nothing is frozen here. At every `s` the profile is re-rendered from `b_ctrl + s*(b_hat - b_ctrl)`,
400 cells are drawn through the production `emit_counts(mode="quasi")`, the delete-1 jackknife is
recomputed on those cells, and `mse_unbiased_capped` / `distance_unbiased` are called. So the
`rho` bound, the #247 cap and the jackknife all move exactly as the scorer would see them.

    uv run python mech_alpha_h1.py [s-grid] [spread] [phi]
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix, vstack

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402
from decompose import terms  # noqa: E402
from msekit import CTRL, TS, group_stats, score  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
NCELL = 400


def main(grid, spread=0.25, phi=0.0, seed=0):
    rp, rs, rc, rj, rg = group_stats(REAL)
    pp, ps, pc, pj, pg = group_stats(PRED)
    assert np.array_equal(rg, pg)
    G = rg.size
    ci = int(np.flatnonzero(rp == CTRL)[0])
    bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
    pmeans = bulk_lognorm_means(ps.astype(np.float64), TS)
    dh = pmeans - bctrl[None, :]
    S = ps.sum(1).astype(np.float64)

    print(f"{'s':>7}{'mse':>10}{'A':>11}{'B':>11}{'claim':>10}{'rho':>7}{'X':>11}"
          f"{'sumD':>10}{'alpha':>8}   (spread={spread}, phi={phi})")
    out = []
    for s in grid:
        t0 = time.time()
        sums = np.zeros_like(ps)
        jks = np.zeros_like(pj)
        for i in range(pp.size):
            if s == 1.0:
                prof = ps[i].astype(np.float64)
            else:
                prof = np.expm1(bctrl + s * dh[i])
                np.clip(prof, 0.0, None, out=prof)
            tot = prof.sum()
            if tot <= 0:
                continue
            prof = prof / tot
            libs = np.full(NCELL, S[i] / NCELL)
            v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(seed * 7919 + i),
                                     mode="quasi", spread=spread, overdispersion=phi)
            X = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                           shape=(NCELL, G))
            sums[i] = np.asarray(X.sum(0)).ravel().astype(np.int64)
            jks[i] = jackknife_correction(X, np.zeros(NCELL, np.intp), 1, TS, chunk=64)[0]
        m, _, _ = score(pp, sums, np.full(pp.size, float(NCELL)), jks, rp, rs, rc, rj, rg)
        t = terms(pp, sums, np.full(pp.size, float(NCELL)), jks, rp, rs, rc, rj, rg)
        print(f"{s:>7.3f}{m:>10.5f}{t['A']:>11.6f}{t['budget']:>11.6f}{t['claim']:>10.6f}"
              f"{t['rho']:>7.4f}{t['X']:>11.6f}{t['sumD']:>10.6f}{1.0796*s:>8.4f}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        out.append((s, m, t["A"], t["budget"], t["claim"], t["rho"], t["X"], t["sumD"]))
    np.save(f"{ROOT}/experiments/096-mse-floor/work/h1_alpha_sp{spread}_phi{phi}.npy",
            np.asarray(out))
    return out


if __name__ == "__main__":
    g = [float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else \
        [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.85, 1.00]
    sp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    ph = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    main(g, sp, ph)
