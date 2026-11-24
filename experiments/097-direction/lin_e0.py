"""e0 -- the incumbent, the two trivial baselines, and the mechanical coverage facts.

LEAK AUDIT for every row here:
  * `incumbent`      : reads the SHIPPED pred.h5ad only.  The reference enters through the
                       scorer (X) and nowhere else.
  * `panel mean`     : the incumbent's own w-weighted mean effect vector, broadcast to every
                       perturbation.  Our prediction only.
  * `raw donor s_p`  : donor delta minus the donor's own panel mean, gene-axis matched by
                       symbol.  Donor only, no whitener, no fit.
  * `raw donor m`    : the donor's panel-mean response vector, broadcast.  Donor only.
"""
from __future__ import annotations

import numpy as np

import lin_kit as K

DONORS = ("k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential")


def main():
    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P, G = t["dh"].shape
        k0, e0, X0, A0, B0, sD = K.shape(t["dh"], t)
        print(f"\n########## {name}: P={P} scored perturbations, G={G} genes, "
              f"banned donor = {t['banned']} ##########")
        print(f"  sumD={sD:.6g}  claim={t['claim']:.6g}  claim/B={t['claim']/B0:.4f}  "
              f"rho={t['rho']:.4f}  scorer mse={t['mse']:.6f}")
        print(f"  INCUMBENT  kappa={k0:.6f}  eps={e0:+.6f}   "
              f"(096's X/sqrt(A*sumD) = {X0/np.sqrt(A0*sD):.6f})")

        rows = []
        Sw = float(t["w"].sum())
        dhbar = (t["w"][:, None] * t["dh"]).sum(0) / Sw
        rows.append(("panel mean of our own dh (constant)",
                     np.repeat(dhbar[None, :], P, 0)))

        for dn in DONORS:
            if dn == t["banned"]:
                continue
            blk = K.donor_block(dn, t["genes"])
            idx = blk["index"]
            cov = np.array([str(x) in idx for x in t["targets"]])
            M = np.zeros((P, G), np.float64)
            for r, x in enumerate(t["targets"]):
                j = idx.get(str(x))
                if j is not None:
                    M[r] = blk["S"][j]
            rows.append((f"raw donor s_p  {dn} [cov {cov.sum()}/{P} perts, "
                         f"{blk['n_genes_hit']}/{G} genes]", M))
            rows.append((f"raw donor mean m  {dn} (constant)",
                         np.repeat(blk["m"][None, :], P, 0)))

        hdr = (f"\n  {'arm':<58}{'kappa':>9}{'M':>7}{'eps':>9}{'cos':>9}"
               f"{'B':>11}{'mse_min':>9}{'s*':>7}")
        print(hdr)
        r0 = K.shape_full(t["dh"], t)
        print(f"  {'INCUMBENT (shipped whitened donor transfer)':<58}"
              f"{r0['kappa']:>9.5f}{1.0:>7.3f}{r0['eps']:>+9.5f}{r0['cos']:>9.5f}"
              f"{r0['B']:>11.4g}{r0['mse_min']:>9.4f}{r0['s_star']:>7.3f}")
        for lbl, DH in rows:
            r = K.shape_full(DH, t)
            print(f"  {lbl:<58}{r['kappa']:>9.5f}{r['kappa']/k0:>7.3f}{r['eps']:>+9.5f}"
                  f"{r['cos']:>9.5f}{r['B']:>11.4g}{r['mse_min']:>9.4f}{r['s_star']:>7.3f}")


if __name__ == "__main__":
    main()
