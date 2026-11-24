"""e1's ceiling, re-read after removing the common mode -- the same test that killed e2.

The winning e1 family is `k562_gwps` embedding, NO intercept, k in {20, 50}, lambda >= 1.
At heavy ridge that model degenerates to a similarity-weighted average of the TRAINING half's
true responses, and a panel average carries a large perturbation-CONSTANT component.  A constant
raises X without touching `B`, so it inflates kappa for free.  On orion the e1 rows carry
eps = +0.165 against the incumbent's +0.012, which is exactly that signature.

So: split each held-out prediction into its w-weighted constant and its centred part, and read
the ceiling off the CENTRED part -- the only part the metric grants a budget for.
"""
from __future__ import annotations

import numpy as np

import lin_kit as K
from lin_e1 import embeddings

SPECS = [("k562_gwps", 20, 1e1), ("k562_gwps", 50, 1e2), ("k562_gwps", 10, 1e1),
         ("both", 20, 1e0), ("orion_hek293t", 50, 1e-2)]
NSPLIT = 30


def main():
    out = {}
    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P = t["dh"].shape[0]
        emb = embeddings(t)
        Y = t["d"]
        splits = []
        for s in range(NSPLIT):
            perm = np.random.default_rng(1000 + s).permutation(P)
            splits.append((perm[: P // 2], perm[P // 2:]))
        print(f"\n################ {name} (P={P}) ################")
        print(f"    {'arm':<34}{'kappa':>9}{'M':>7}{'eps':>10}"
              f"{'kappa_c':>10}{'M_c':>7}{'const cos':>11}")
        # incumbent, whole and centred, on the same held-out halves
        ik, ick = [], []
        for _, te in splits:
            ik.append(K.shape(t["dh"], t, te)[0])
            w = t["w"][te]
            c = (w[:, None] * t["dh"][te]).sum(0) / w.sum()
            ick.append(K.shape(t["dh"] - c[None, :], t, te)[0])
        ik, ick = np.array(ik), np.array(ick)
        print(f"    {'INCUMBENT':<34}{ik.mean():>9.5f}{1.0:>7.3f}"
              f"{K.shape_full(t['dh'], t)['eps']:>+10.5f}{ick.mean():>10.5f}{1.0:>7.3f}"
              f"{'-':>11}")
        for dn, k, lam in SPECS:
            Ef = emb[dn]["E"]
            kc = Ef.shape[1] // (2 if dn == "both" else 1)
            E = np.hstack([Ef[:, :k], Ef[:, kc:kc + k]]) if dn == "both" else Ef[:, :k]
            kk, kkc, ee, cc = [], [], [], []
            for tr, te in splits:
                Et, Yt = E[tr], Y[tr]
                scale = float(np.einsum("ij,ij->", Et, Et)) / max(E.shape[1], 1)
                W, c0 = K.ridge_fit(Et, Yt, lam * max(scale, 1e-30), False)
                DH = np.zeros_like(t["dh"])
                DH[te] = E[te] @ W + c0
                r = K.shape_full(DH, t, te)
                kk.append(r["kappa"])
                ee.append(r["eps"])
                w = t["w"][te]
                cst = (w[:, None] * DH[te]).sum(0) / w.sum()
                DHc = DH.copy()
                DHc[te] -= cst
                kkc.append(K.shape(DHc, t, te)[0])
                nc = np.linalg.norm(cst)
                dbar = (w[:, None] * t["d"][te]).sum(0) / w.sum()
                cc.append(float(cst @ dbar / max(nc * np.linalg.norm(dbar), 1e-30)))
            kk, kkc = np.array(kk), np.array(kkc)
            lbl = f"{dn} k={k} ic=F lam={lam:.0e}"
            print(f"    {lbl:<34}{kk.mean():>9.5f}{(kk/ik).mean():>7.3f}"
                  f"{np.mean(ee):>+10.5f}{kkc.mean():>10.5f}{(kkc/ick).mean():>7.3f}"
                  f"{np.mean(cc):>+11.4f}")
            out.setdefault(lbl, {})[name] = ((kk / ik).mean(), (kkc / ick).mean())
        K._CACHE.clear()

    print("\n===== e1 ceiling, cross-surface (rule 1: min over harnesses) =====")
    print(f"  {'arm':<34}{'M h1':>9}{'M ori':>9}{'M use':>9}"
          f"{'Mc h1':>9}{'Mc ori':>9}{'Mc use':>9}")
    for lbl, d in out.items():
        (m1, c1), (m2, c2) = d["h1-test"], d["orion"]
        print(f"  {lbl:<34}{m1:>9.3f}{m2:>9.3f}{min(m1,m2):>9.3f}"
              f"{c1:>9.3f}{c2:>9.3f}{min(c1,c2):>9.3f}")


if __name__ == "__main__":
    main()
