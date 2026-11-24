"""Step 0: decompose X into a COMMON-MODE leg and a DIFFERENTIAL leg, on both harnesses.

The master identity is  MSE - 1 = (A - 2X + Kc - rho*claim)/sumD  and the only term a better
model moves is X = sum_p w_p <dh_p, d_p>.  Every route in this experiment is a re-weighting or
re-targeting of dh, so the first thing to establish is the geometry of dh and d.

Write the w-weighted means  dhbar = sum_p w_p dh_p / sum_p w_p   and   dbar likewise.  Then

    X   = Sw * <dhbar, dbar>  +  X_d          X_d = sum_p w_p <dh_p - dhbar, d_p - dbar>
    A   = Sw * ||dhbar||^2    +  B            B   = sum_p w_p ||dh_p - dhbar||^2   (~ the #348 budget)
    R   = Sw * ||dbar||^2     +  R_d

so the prediction has two independent jobs -- get the PANEL-MEAN response right, and get each
perturbation's DEVIATION from it right -- and our pipeline currently ties them to a single
scalar alpha.  This script measures how much of the reference's signal lives in each.
"""
from __future__ import annotations
import sys, os
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats

HARNESS = [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion-thinned", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
]


def load(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    assert np.array_equal(rg, pg)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    # zero each row's own target coordinate so plain einsums reproduce the scorer's sums
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0
        t["d"][k, c] = 0.0
    t["targets"] = [str(x) for x in pp[t["keep"]]]
    return t


def main():
    for name, real, pred in HARNESS:
        t = load(real, pred)
        dh, d, w = t["dh"], t["d"], t["w"]
        Sw = float(w.sum())
        dhbar = (w[:, None] * dh).sum(0) / Sw
        dbar = (w[:, None] * d).sum(0) / Sw
        ddh = dh - dhbar
        dd = d - dbar
        X = float((np.einsum("ij,ij->i", dh, d) * w).sum())
        A = float((np.einsum("ij,ij->i", dh, dh) * w).sum())
        R = float((np.einsum("ij,ij->i", d, d) * w).sum())
        Xc = Sw * float(dhbar @ dbar)
        Ac = Sw * float(dhbar @ dhbar)
        Rc = Sw * float(dbar @ dbar)
        Xd = float((np.einsum("ij,ij->i", ddh, dd) * w).sum())
        Ad = float((np.einsum("ij,ij->i", ddh, ddh) * w).sum())
        Rd = float((np.einsum("ij,ij->i", dd, dd) * w).sum())
        sD = t["sumD"]
        print(f"\n===== {name}: n={t['n']}  sumD={sD:.6g}  claim={t['claim']:.6g}  "
              f"budget B={t['budget']:.6g} =====")
        print(f"  identity checks: X {X:.6g} vs {t['X']:.6g} | A {A:.6g} vs {t['A']:.6g} "
              f"| R {R:.6g} vs {t['R']:.6g}")
        print(f"  split check: Xc+Xd={Xc+Xd:.6g}  Ac+Ad={Ac+Ad:.6g}  Rc+Rd={Rc+Rd:.6g}")
        print(f"  Ad vs scorer budget B: {Ad:.6g} vs {t['budget']:.6g} "
              f"(ratio {Ad/t['budget']:.4f})")
        print()
        print(f"  {'leg':<14}{'X':>12}{'A':>12}{'R':>12}{'cos=X/sqrt(AR)':>16}{'X^2/A':>12}")
        for lbl, xx, aa, rr in (("COMMON", Xc, Ac, Rc), ("DIFFERENTIAL", Xd, Ad, Rd),
                                ("total", X, A, R)):
            cs = xx / np.sqrt(aa * rr) if aa > 0 and rr > 0 else float("nan")
            print(f"  {lbl:<14}{xx:>12.6g}{aa:>12.6g}{rr:>12.6g}{cs:>16.4f}{xx*xx/aa:>12.6g}")
        print()
        print(f"  reference's common-mode share of R:  {100*Rc/R:.2f} %")
        print(f"  our  common-mode share of A:         {100*Ac/A:.2f} %")
        print(f"  common-mode share of X:              {100*Xc/X:.2f} %")
        print(f"  cos(dhbar, dbar)                     {float(dhbar@dbar)/np.sqrt((dhbar@dhbar)*(dbar@dbar)):+.4f}")
        print(f"  ||dhbar|| {np.linalg.norm(dhbar):.5f}   ||dbar|| {np.linalg.norm(dbar):.5f}"
              f"   optimal common scale c* = <dhbar,dbar>/||dhbar||^2 = "
              f"{float(dhbar@dbar)/float(dhbar@dhbar):.4f}")
        # what an ORACLE offset would be worth, in numerator units and in raw mse
        oracle_gain = Sw * float(dbar @ dbar)                   # numerator reduction
        build_gain = Sw * float(dhbar @ dbar) ** 2 / float(dhbar @ dhbar)
        print(f"  ORACLE constant-offset numerator gain  {oracle_gain:.6g}  "
              f"= {oracle_gain/sD:.4f} raw mse")
        print(f"  BUILDABLE (rescale our own dhbar)      {build_gain:.6g}  "
              f"= {build_gain/sD:.4f} raw mse   ({100*build_gain/oracle_gain:.1f} % of oracle)")
        np.savez(f"{ROOT}/experiments/097-direction/work/{name}-base.npz",
                 X=X, A=A, R=R, Xc=Xc, Ac=Ac, Rc=Rc, Xd=Xd, Ad=Ad, Rd=Rd, Sw=Sw,
                 sumD=sD, claim=t["claim"], budget=t["budget"], Kc=t["Kc"], mse=t["mse"])


if __name__ == "__main__":
    main()
