"""Does reference sampling noise attenuate the DECONFOUNDED cosine? Test, do not assume.

The claim under test: both harnesses read a lower cosine than the board because their reference
is noisier, and noise attenuates measured agreement. That is certainly true of the PLUG-IN cosine
X/sqrt(A*R), whose denominator carries C_p. It is NOT obviously true of X/sqrt(A*sumD), because
sumD = R - K - Kc has the reference's sampling variance subtracted off.

Direct test: thin the reference's cells per perturbation and re-measure. If the correction works,
the deconfounded cosine is flat in reference depth while the plug-in one falls.
"""
import sys
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix, issparse
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction
from cell_eval2.metrics.delta import _row_weights
from cell_eval2.prep import bulk_lognorm_means
from msekit import CTRL, TS, group_stats

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"

A = ad.read_h5ad(REAL)
lab = np.asarray(A.obs["target_gene"].astype(str).values)
genes = np.asarray([str(g) for g in A.var_names]); G = genes.size
X = A.X.tocsr() if issparse(A.X) else csr_matrix(A.X)
ci = np.flatnonzero(lab == CTRL)
Xc = csr_matrix(X[ci])
ctrl_sum = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
Cctrl = jackknife_correction(Xc, np.zeros(ci.size, np.intp), 1, TS, chunk=64)[0]
bctrl = bulk_lognorm_means(ctrl_sum[None, :], TS)[0]
del Xc

pp, ps, pc, pj, pg = group_stats(PRED)
pm = bulk_lognorm_means(ps.astype(np.float64), TS)
pidx = {str(p): i for i, p in enumerate(pp)}
gpos = {g: i for i, g in enumerate(genes)}
targets = [t for t in sorted(set(lab) - {CTRL}) if t in pidx]

print(f"h1-test: {len(targets)} perturbations, control {ci.size} cells, C_ctrl {Cctrl:.4f}")
print(f"{'ref cells/pert':>15}{'A':>11}{'X':>12}{'R':>11}{'K':>11}{'sumD':>11}"
      f"{'plug-in cos':>13}{'deconf cos':>12}")
rng = np.random.default_rng(0)
for ncell in (400, 200, 100, 50):
    rows_k, dh_l, d_l, Cp_l, cols = [], [], [], [], {}
    for k, t in enumerate(targets):
        r = np.flatnonzero(lab == t)
        if ncell < r.size:
            r = np.sort(rng.choice(r, ncell, replace=False))
        sub = csr_matrix(X[r])
        rs = np.asarray(sub.sum(0)).ravel().astype(np.float64)
        Cp_l.append(jackknife_correction(sub, np.zeros(r.size, np.intp), 1, TS, chunk=64)[0])
        d_l.append(bulk_lognorm_means(rs[None, :], TS)[0] - bctrl)
        dh_l.append(pm[pidx[t]] - bctrl)
        if t in gpos: cols[k] = gpos[t]
    dh = np.asarray(dh_l); d = np.asarray(d_l); Cp = np.asarray(Cp_l)
    w = _row_weights(len(targets), G, cols)
    a = np.einsum("ij,ij->i", dh, dh); x = np.einsum("ij,ij->i", dh, d)
    r2 = np.einsum("ij,ij->i", d, d)
    for k, c in cols.items():
        a[k] -= dh[k, c] ** 2; x[k] -= dh[k, c] * d[k, c]; r2[k] -= d[k, c] ** 2
    Aa = float((a * w).sum()); Xx = float((x * w).sum()); Rr = float((r2 * w).sum())
    Kk = float((Cp * w).sum()); Kc = float(np.full(len(targets), Cctrl) @ w)
    sD = Rr - Kk - Kc
    print(f"{ncell:>15}{Aa:>11.6f}{Xx:>12.6f}{Rr:>11.6f}{Kk:>11.6f}{sD:>11.6f}"
          f"{Xx/np.sqrt(Aa*Rr):>13.5f}{Xx/np.sqrt(Aa*sD):>12.5f}", flush=True)
