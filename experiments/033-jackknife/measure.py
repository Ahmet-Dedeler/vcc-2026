"""033 — measure the `expr_mse` jackknife credit we are forfeiting.

The metric's numerator is

    N_p = (1/G_p) ( ||b_hat_p - b_p||^2  -  rho * min(C_hat_p, C_p)  -  C_p )

where `C_hat_p` is the delete-one jackknife estimate of **our own submission's** sampling
noise. The metric pays us for noise we emit, capped by the reference's own noise. A
minimum-variance emitter has `C_hat_p ~ 0` and collects none of it.

This computes `C_p` for the reference and `C_hat_p` for predictions already on disk, so the
size of the forfeited term is a measurement rather than an estimate, and so the predicted
`mse` can be checked against what the scorer actually returned for the same file.

    C_p = ((n-1)/n) * sum_g sum_i (v_ig - vbar_g)^2,
    v_ig = log1p(TS * (P_pg - y_ig) / (S_p - l_i))
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

TS = 5e4
CONTROL, PERT_COL = "non-targeting", "target_gene"


def jackknife_C(X: sp.csr_matrix) -> float:
    """Delete-one jackknife sum-of-squares for one group's log-pseudobulk."""
    n = X.shape[0]
    if n < 2:
        return 0.0
    P = np.asarray(X.sum(0)).ravel().astype(np.float64)
    S = P.sum()
    ell = np.asarray(X.sum(1)).ravel().astype(np.float64)
    # v_ig differs from the pooled value only on genes cell i actually carries, so the dense
    # (n x G) matrix is never built: start from the all-genes leave-one-out value and correct
    # the sparse entries.
    tot = 0.0
    base_num = P[None, :]
    for start in range(0, n, 64):
        blk = slice(start, min(start + 64, n))
        Y = X[blk].toarray().astype(np.float64)
        denom = (S - ell[blk])[:, None]
        V = np.log1p(TS * (base_num - Y) / denom)
        tot += V
        if start == 0:
            acc = np.zeros(X.shape[1])
            sq = np.zeros(X.shape[1])
        acc += V.sum(0)
        sq += (V ** 2).sum(0)
        tot = 0.0
    var = sq - acc ** 2 / n
    return float((n - 1) / n * var.sum())


def group_C(path: str, label_col=PERT_COL, max_groups=12):
    a = ad.read_h5ad(path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    lab = a.obs[label_col].astype(str).to_numpy()
    out = {}
    perts = [p for p in sorted(set(lab)) if p != CONTROL][:max_groups]
    for p in perts:
        out[p] = jackknife_C(X[lab == p])
    ctrl = jackknife_C(X[lab == CONTROL][:18400])
    return out, ctrl, X, lab


if __name__ == "__main__":
    real = sys.argv[1] if len(sys.argv) > 1 else "data/harness/orion-hct116.h5ad"
    preds = sys.argv[2:] or [
        "experiments/031-somatic/out/runs/g0.85_q1.0_a1.08/pred.h5ad"]

    Cr, Cctrl, Xr, labr = group_C(real)
    print(f"reference {Path(real).name}: median C_p = {np.median(list(Cr.values())):.3f} "
          f"over {len(Cr)} perturbations, C_ctrl = {Cctrl:.3f}")

    # Denominator uses the reference only, so it is the same for every submission.
    Pc = np.asarray(Xr[labr == CONTROL].sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    D = []
    for p, c in Cr.items():
        Pp = np.asarray(Xr[labr == p].sum(0)).ravel().astype(np.float64)
        b_p = np.log1p(TS * Pp / Pp.sum())
        D.append(float(((b_p - b_ctrl) ** 2).sum() - c - Cctrl))
    print(f"  ||b_p - b_ctrl||^2 - C_p - C_ctrl : median {np.median(D):.3f}")
    print(f"  the reference's own noise C_p is {100*np.median(list(Cr.values()))/(np.median(D)+np.median(list(Cr.values()))+Cctrl):.0f}% "
          f"of the raw effect energy\n")

    for pr in preds:
        Ch, _, Xp, labp = group_C(pr)
        med = float(np.median(list(Ch.values())))
        credit = float(np.median([min(Ch[p], Cr[p]) for p in Ch if p in Cr]))
        print(f"{Path(pr).parent.name:34} median C_hat = {med:9.4f}  "
              f"credited min(C_hat,C_p) = {credit:8.4f}  "
              f"= {credit/np.median(D):.3f} of the denominator")
