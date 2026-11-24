"""STEP 1 -- does the analytic surrogate rank the same as the realized emission?

`dh_analytic = log1p(TS*profile) - b_ctrl` assumes the quasi emitter's exact-total repair makes
the group sum equal `n*L*profile`.  If that holds, arms can be ranked without emitting a single
count.  The test: rebuild the SHIPPED config's `E` and compare its `dh_analytic` against the
`dh` the scorer actually reads off `pred.h5ad`.
"""
from __future__ import annotations

import numpy as np
from wkit import ROOT, Ctx, uniform


def main():
    for name in ("h1-test", "orion"):
        c = Ctx(name)
        print(f"\n########## {name}  n={len(c.perts)}  G={c.genes.size} ##########")
        print(f"  wcols rebuilt and matched cache: n={len(c.wcols)} "
              f"sha256={c.wcols_sha[:16]}...")
        print(f"  b_ctrl(build) vs b_ctrl(scorer): max|diff| = "
              f"{np.abs(c.b_ctrl - c.bctrl_scorer).max():.3e}")
        r = c.run(uniform)
        dhA, dhR = r["dh"], c.dh_real
        XA, AA, WA = r["X"], r["A"], r["W"]
        XR, AR, WR = c.XAW(dhR)

        nA = np.linalg.norm(dhA, axis=1)
        nR = np.linalg.norm(dhR, axis=1)
        cos = np.einsum("ij,ij->i", dhA, dhR) / np.maximum(nA * nR, 1e-30)
        print(f"\n  per-perturbation cosine(dh_analytic, dh_realized):")
        print(f"    median {np.median(cos):.6f}   min {cos.min():.6f}   "
              f"p05 {np.percentile(cos, 5):.6f}   max {cos.max():.6f}")
        print(f"  per-perturbation ||dh|| ratio analytic/realized: median "
              f"{np.median(nA/nR):.6f}  min {(nA/nR).min():.6f}  max {(nA/nR).max():.6f}")
        print(f"\n  {'':<12}{'analytic':>14}{'realized':>14}{'ratio':>10}")
        print(f"  {'X':<12}{XA:>14.8g}{XR:>14.8g}{XA/XR:>10.5f}")
        print(f"  {'A':<12}{AA:>14.8g}{AR:>14.8g}{AA/AR:>10.5f}")
        print(f"  {'W=X^2/A':<12}{WA:>14.8g}{WR:>14.8g}{WA/WR:>10.5f}")
        print(f"  sqrt(W_analytic / W_realized) = {np.sqrt(WA/WR):.5f}   "
              f"<-- the surrogate's bias on the deliverable metric")
        print(f"  aggregate cosine X/sqrt(A*R): analytic "
              f"{XA/np.sqrt(AA*c.t['R']):.5f}  realized {XR/np.sqrt(AR*c.t['R']):.5f}")
        print(f"  sumD={c.sumD:.6g} claim={c.claim:.6g} Kc={c.Kc:.6g}  "
              f"mse(scorer)={c.t['mse']:.6f}")
        np.savez(ROOT / f"experiments/097-direction/work/{name}-validate.npz",
                 XA=XA, AA=AA, WA=WA, XR=XR, AR=AR, WR=WR, cos=cos, nA=nA, nR=nR,
                 sumD=c.sumD, claim=c.claim, Kc=c.Kc)
        del c


if __name__ == "__main__":
    main()
