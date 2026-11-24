"""155 — is there ANY target-specific direction signal in the control covariance?

`score.py` says no at k=20/50 PCs and says the k=0 answer is anti-correlated.  Both readings
are confounded, in opposite ways:

  * k=0 is dominated by COMPOSITION.  CPM normalisation fixes each cell's total, so
    sum_g cov(T, g) = 0 by construction and most genes' coefficients carry the same global
    sign.  That is a property of the normalisation, not of gene T.
  * k>=20 removes the top principal components, which is where real shared transcriptional
    programmes live.  Over-correcting can delete the signal along with the confound.

The clean separation does not need PCs at all: the confound is the part of beta[.,g] that is
the SAME for every target, so subtract the mean over the 132 targets (column centering).  What
survives is by construction target-specific.  Row centering removes each target's own
compositional offset.  This sweeps both, at every k, and stratifies by how well the target is
actually measured -- a target at 4 CPM and 12k UMI is ~0.05 counts/cell and cannot carry a
coefficient at all.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import head_purity, load_frame  # noqa: E402

s, inc, gdst, pok = load_frame()


def frame(M):
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def quick(P, pmask=None):
    """sign accuracy + head purity, optionally on a subset of perturbations."""
    s.pred_lfc = P
    if pmask is None:
        return s.sign_accuracy(), head_purity(s, P)[0]
    keep = np.zeros_like(s.real_sig)
    keep[pmask] = True
    rs, ad_ = s.real_sig.copy(), s.adj.copy()
    s.real_sig, s.adj = rs & keep, ad_ & keep
    out = (s.sign_accuracy(), head_purity(s, P)[0])
    s.real_sig, s.adj = rs, ad_
    return out


def centered(M, mode):
    M = M.copy()
    if mode in ("col", "double"):
        M -= M.mean(0, keepdims=True)        # remove the all-target consensus per gene
    if mode in ("row", "double"):
        M -= M.mean(1, keepdims=True)        # remove each target's compositional offset
    return M


def main() -> None:
    z20 = np.load(os.path.join(HERE, "work", "coef_k20.npz"), allow_pickle=True)
    cpm = z20["target_cpm"]

    print("A. centering sweep -- sign accuracy of the -beta candidate "
          "(0.5 = chance; +beta is 1-x)")
    print(f"{'k PCs':>6} {'raw':>18} {'row-centered':>18} {'col-centered':>18} "
          f"{'double-centered':>18}")
    print(f"{'':>6} " + " ".join(f"{'acc / head':>18}" for _ in range(4)))
    for k in (0, 5, 10, 20, 50):
        p = os.path.join(HERE, "work", f"coef_k{k}.npz")
        if not os.path.exists(p):
            continue
        B = np.load(p, allow_pickle=True)["beta"].astype(np.float64)
        cells = []
        for mode in (None, "row", "col", "double"):
            a, h = quick(frame(-centered(B, mode) if mode else -B))
            cells.append(f"{a:.4f} / {h:.4f}")
        print(f"{k:>6} " + " ".join(f"{c:>18}" for c in cells))

    print()
    print("B. stratified by how well the target itself is measured (k=20, double-centered)")
    B = np.load(os.path.join(HERE, "work", "coef_k20.npz"),
                allow_pickle=True)["beta"].astype(np.float64)
    P = frame(-centered(B, "double"))
    B0 = np.load(os.path.join(HERE, "work", "coef_k0.npz"),
                 allow_pickle=True)["beta"].astype(np.float64)
    P0 = frame(-centered(B0, "double"))
    order = np.argsort(-cpm)
    print(f"{'target stratum':<28} {'n':>4} {'CPM range':>16} "
          f"{'k20 acc/head':>18} {'k0 acc/head':>18}")
    for lo, hi, name in ((0, 33, "top quartile by CPM"), (33, 66, "2nd quartile"),
                         (66, 99, "3rd quartile"), (99, 132, "bottom quartile")):
        m = order[lo:hi]
        a20, h20 = quick(P, m)
        a0, h0 = quick(P0, m)
        print(f"{name:<28} {m.size:>4} {cpm[m].min():7.1f}-{cpm[m].max():<8.1f} "
              f"{f'{a20:.4f} / {h20:.4f}':>18} {f'{a0:.4f} / {h0:.4f}':>18}")

    print()
    print("C. how good could it be?  Spearman |corr| of the candidate against the reference "
          "lfc, per perturbation (k=20 double-centered)")
    from scipy.stats import spearmanr
    rs = []
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 20:
            continue
        rs.append(spearmanr(P[p, i], s.real_lfc[p, i]).statistic)
    rs = np.array(rs)
    print(f"  n={rs.size}  mean {rs.mean():+.4f}  median {np.median(rs):+.4f}  "
          f"frac>0 {(rs > 0).mean():.3f}  (a real signal would sit well above 0)")
    print(f"  one-sample t against 0: t = {rs.mean() / (rs.std(ddof=1) / np.sqrt(rs.size)):+.2f}")


if __name__ == "__main__":
    main()
