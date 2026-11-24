"""Does the master formula's A~=B shortcut hold for the Miller PANEL-MEAN family? Verify.

MSE(s) = 1 + [s^2 A' - 2 s X' + Kc - min(claim, s^2 B')]/sumD, minimised over s, has optimum at
s_c = sqrt(claim/B') with value

    1 + [ claim*(A'-B')/B'  +  Kc  -  2 sqrt(claim/B') X' ] / sumD

The shortcut MSE_min = 1 + (Kc - 2 sqrt(claim*W))/sumD, W = X'^2/A', drops the FIRST term by
assuming A' = B'. `_across_pert_budget` is a CENTRED sum of squares, so shrinking toward the panel
mean collapses B' while leaving A' -- exactly where the shortcut is least safe. Measured here.
"""
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.metrics.delta import _across_pert_budget
from cell_eval2.prep import bulk_lognorm_means
from decompose import terms
from msekit import TS, group_stats
from mech_miller import donor_key

def run(real, pred, name):
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    dh, d, w, cols = t["dh"], t["d"], t["w"], t["cols"]
    for k, c in cols.items(): dh[k, c] = 0.0; d[k, c] = 0.0
    bctrl, claim, Kc, sD = t["bctrl"], t["claim"], t["Kc"], t["sumD"]
    targets = [str(x) for x in pp[t["keep"]]]
    dk = donor_key(rg, targets)
    rank = np.argsort(np.argsort(-dk, axis=1), axis=1)
    dbar = dh.mean(0)

    def stats(D):
        A = float((np.einsum("ij,ij->i", D, D)*w).sum())
        X = float((np.einsum("ij,ij->i", D, d)*w).sum())
        B = _across_pert_babel(D) if False else _across_pert_budget(bctrl[None, :]+D, cols, w)
        return A, B, X

    def exact(A, B, X):
        """min over s of the true expression, both rho branches."""
        best = (None, 9e9)
        for s in np.linspace(0.02, 3.0, 1500):
            v = 1 + (s*s*A - 2*s*X + Kc - min(claim, s*s*B))/sD
            if v < best[1]: best = (s, v)
        return best

    A0, B0, X0 = stats(dh)
    s0, m0 = exact(A0, B0, X0)
    print(f"\n===== {name} =====")
    print(f"  UNSHRUNK   (A-B)/A = {(A0-B0)/A0:+.4f}   exact min mse {m0:.5f} at s={s0:.3f}")
    print(f"  {'family':<34}{'(A-B)/A':>10}{'exact min':>11}{'d vs unshrunk':>15}"
          f"{'shortcut':>10}{'sqrt(W/Wg)':>12}")
    Wg = X0*X0/A0
    for tgt, lbl in ((dbar, "toward PANEL MEAN"), (np.zeros_like(dbar), "toward CONTROL")):
        for g in (0.5, 1.0, 2.0):
            a = np.exp(-g*rank/2000.0)
            D = a*dh + (1-a)*tgt[None, :]
            A, B, X = stats(D)
            s, m = exact(A, B, X)
            W = X*X/A
            short = 1 + (Kc - 2*np.sqrt(claim*W))/sD
            print(f"  {lbl+f', gamma={g}':<34}{(A-B)/A:>10.4f}{m:>11.5f}{m-m0:>+15.5f}"
                  f"{short:>10.5f}{np.sqrt(W/Wg):>12.3f}")

for r, p, n in [(f"{ROOT}/data/harness/h1-test.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","h1-test"),
  (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","orion-thinned")]:
    run(r, p, n)
