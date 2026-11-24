"""Is the EMITTED DEPTH a lever on X?  An accident in a reconstruction says it might be.

`emit_counts(mode="quasi")` sets `L = int(median(lib_sizes))` and `totals = floor(n*L*profile)`.
Rebuilding the orion prediction with UNIFORM libs (median = mean = 13,843) instead of SAMPLED real
control libs (median 12,174) moved X from 0.008136 to 0.005072 -- a 38 % drop from a 14 % depth
change, with the per-gene count totals drifting only 0.27 % in L1.

The mechanism has to be the FLOOR.  A gene whose `n*L*p_g < 1` is emitted as exactly zero, so its
predicted dh_g = -b_ctrl_g; and the reference's own low-expression genes sit near zero too, so
d_g is also near -b_ctrl_g.  Emitting a zero where the reference has a zero is a genuine positive
contribution to X that has nothing to do with biology -- it is DEPTH MATCHING.  If that is what is
happening, the emitted depth is a free, buildable lever on the direction term.

Fast because X and A need only the group SUMS -- no jackknife.  Nothing about `claim` is measured
here, and `claim` is exactly what a depth change also moves, so this bounds a lever rather than
pricing one.
"""
from __future__ import annotations
import sys, time
import numpy as np
from scipy.sparse import csr_matrix
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget
from decompose import terms
from msekit import CTRL, TS, group_stats
from vcc2026.donor_transfer import emit_counts
import anndata as ad

which = sys.argv[1]
REAL, PRED = {
    "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
           f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
              f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
}[which]
rp, rs, rc, rj, rg = group_stats(REAL)
pp, ps, pc, pj, pg = group_stats(PRED)
G = rg.size
ci = int(np.flatnonzero(rp == CTRL)[0])
bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
dh0 = bulk_lognorm_means(ps.astype(np.float64), TS) - bctrl[None, :]
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
keep, cols, w, d, sD = t["keep"], t["cols"], t["w"], t["d"].copy(), t["sumD"]
for k, c in cols.items():
    d[k, c] = 0.0
R = ad.read_h5ad(REAL, backed="r")
rl = np.asarray(R.obs["target_gene"].astype(str).values)
libs_ctrl = np.asarray(R.X[np.flatnonzero(rl == CTRL)].sum(1)).ravel().astype(np.float64)
mp = np.flatnonzero(rl != CTRL)
libs_pert = np.asarray(R.X[mp[:4000]].sum(1)).ravel().astype(np.float64)
R.file.close()
L_ctrl, L_pert = float(np.median(libs_ctrl)), float(np.median(libs_pert))
print(f"# {which}: reference control median lib {L_ctrl:.0f}, reference PERTURBATION median lib "
      f"{L_pert:.0f} (ratio {L_pert/L_ctrl:.3f})", flush=True)
print(f"# cached prediction: X={t['X']:.6f} A={t['A']:.6f} kappa={t['X']/np.sqrt(t['budget']*sD):.5f}",
      flush=True)
NC = int(round(float(np.median(pc))))
print(f"\n{'L mult':>7}{'L used':>9}{'nonzero/gene-row':>18}{'X':>11}{'A':>11}{'B':>11}"
      f"{'kappa':>10}{'M':>8}", flush=True)
k_ref = t["X"] / np.sqrt(t["budget"] * sD)
for mult in (0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.5, 4.0):
    t0 = time.time()
    sums = np.zeros_like(ps)
    for i in range(pp.size):
        prof = np.expm1(bctrl + dh0[i]); np.clip(prof, 0.0, None, out=prof)
        tot = prof.sum()
        if tot <= 0:
            continue
        prof = prof / tot
        libs = np.full(NC, L_ctrl * mult)
        v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(7919 + i),
                                 mode="quasi", spread=0.25, overdispersion=0.0)
        M = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                       shape=(NC, G))
        sums[i] = np.asarray(M.sum(0)).ravel().astype(np.int64)
    pm = bulk_lognorm_means(sums.astype(np.float64), TS)
    dh = pm[keep] - bctrl[None, :]
    for k, c in cols.items():
        dh[k, c] = 0.0
    X = float((np.einsum("ij,ij->i", dh, d) * w).sum())
    A = float((np.einsum("ij,ij->i", dh, dh) * w).sum())
    B = float(_across_pert_budget(dh, cols, w))
    kap = X / np.sqrt(B * sD)
    nz = float((sums[keep] > 0).sum()) / keep.size
    print(f"{mult:>7.2f}{L_ctrl*mult:>9.0f}{nz:>18.0f}{X:>11.6f}{A:>11.6f}{B:>11.6f}"
          f"{kap:>10.5f}{kap/k_ref:>8.3f}   [{time.time()-t0:.0f}s]", flush=True)
