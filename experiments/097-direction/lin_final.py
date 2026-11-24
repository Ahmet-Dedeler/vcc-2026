"""The consolidated table: every arm, both harnesses, whole and common-mode-removed.

Ratios are quoted BOTH ways because the two answer different questions:
  `M`   = kappa(arm) / kappa(incumbent)                -- what the scorer would read
  `Mc`  = kappa(arm centred) / kappa(incumbent centred) -- what the arm knows about the
          perturbation, with the free constant taken off BOTH sides

`Mc` is the honest comparison for any arm whose `eps` differs from the incumbent's, because a
perturbation-constant raises X while leaving `_across_pert_budget` untouched, i.e. it buys kappa
for nothing.  Every arm in this experiment that beat the incumbent on `M` and lost on `Mc` was
buying exactly that.
"""
from __future__ import annotations

import numpy as np

import lin_kit as K
from lin_null import WIN, donor_S_on_targets, fit_arm

DON = ("k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential")


def centred(DH, t):
    w = t["w"]
    return DH - ((w[:, None] * DH).sum(0) / w.sum())[None, :]


def main():
    store = {}
    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P = t["dh"].shape[0]
        r0 = K.shape_full(t["dh"], t)
        k0 = r0["kappa"]
        k0c = K.shape(centred(t["dh"], t), t)[0]
        N = float((t["Chat"] * t["w"]).sum()) * (1 - 1 / P)
        print(f"\n################ {name}  P={P}  G={t['genes'].size} ################")
        print(f"  incumbent kappa {k0:.5f} (centred {k0c:.5f})  eps {r0['eps']:+.5f}  "
              f"B {r0['B']:.5f}  emitter noise in B {N:.5f} ({N/r0['B']*100:.1f}%)")
        print(f"  claim {t['claim']:.5f}  sumD {r0['sumD']:.5f}  claim/B {t['claim']/r0['B']:.3f}"
              f"   (challenge at the #247 cap: {K.CH['claim']/K.CH['B']:.3f})")
        print(f"\n  {'arm':<44}{'kappa':>10}{'M':>7}{'eps':>10}{'kap_c':>10}{'Mc':>7}"
              f"{'mse_min':>9}{'cov':>6}")
        print(f"  {'INCUMBENT (whitened g=0.65, emitted)':<44}{k0:>10.5f}{1.0:>7.3f}"
              f"{r0['eps']:>+10.5f}{k0c:>10.5f}{1.0:>7.3f}{r0['mse_min']:>9.4f}{P:>6}")

        arms = {}
        dhbar = ((t["w"][:, None] * t["dh"]).sum(0) / t["w"].sum())
        arms["panel mean of our own dh (constant)"] = (np.repeat(dhbar[None, :], P, 0), P)
        Ms = {}
        for dn in DON:
            if dn == t["banned"]:
                continue
            M, cov, med, _ = donor_S_on_targets(t, dn)
            Ms[dn] = (M / med, cov)
            arms[f"raw donor s_p  {dn}"] = (M, int(cov.sum()))
        acc = np.zeros((P, t["genes"].size))
        cnt = np.zeros(P)
        for dn in ("k562_gwps", "orion_hek293t"):
            Mn, cov = Ms[dn]
            acc[cov] += Mn[cov]
            cnt[cov] += 1
        acc[cnt > 0] /= cnt[cnt > 0, None]
        arms["UNWHITENED mean k562+hek293t (g=0)"] = (acc, int((cnt > 0).sum()))
        arms["e2 winner k562->hek293t k50 ic=T l=1"] = (fit_arm(t, WIN)[0], int((cnt > 0).sum()))
        arms["e2 winner, no intercept"] = (fit_arm(t, dict(WIN, icept=False))[0],
                                           int((cnt > 0).sum()))

        for lbl, (DH, cov) in arms.items():
            r = K.shape_full(DH, t)
            kc = K.shape(centred(DH, t), t)[0]
            kap = r["kappa"]
            print(f"  {lbl:<44}{kap:>10.5f}{kap/k0:>7.3f}{r['eps']:>+10.5f}"
                  f"{kc:>10.5f}{kc/k0c:>7.3f}{r['mse_min']:>9.4f}{cov:>6}")
            store.setdefault(lbl, {})[name] = (kap / k0, kc / k0c, r["eps"], r["mse_min"])
        K._CACHE.clear()

    print("\n\n===== CROSS-SURFACE (rule 1: min over harnesses; a sign flip kills the arm) =====")
    print(f"  {'arm':<44}{'M h1':>8}{'M ori':>9}{'M min':>8}"
          f"{'Mc h1':>8}{'Mc ori':>9}{'Mc min':>8}{'ch raw':>9}{'verdict':>14}")
    base = K.challenge_raw(1.0, 0.01424)
    print(f"  {'INCUMBENT':<44}{1.0:>8.3f}{1.0:>9.3f}{1.0:>8.3f}"
          f"{1.0:>8.3f}{1.0:>9.3f}{1.0:>8.3f}{base:>9.4f}{'reference':>14}")
    for lbl, d in store.items():
        (m1, c1, e1, _), (m2, c2, e2, _) = d["h1-test"], d["orion"]
        if not np.isfinite(m1) or not np.isfinite(m2):
            print(f"  {lbl:<44}{'B=0 -- kappa undefined for a perturbation-constant arm':>90}")
            continue
        raw = K.challenge_raw(min(c1, c2), max(e1, e2))
        v = ("DEAD sign" if np.sign(c1) != np.sign(c2)
             else ("BEATS e0" if min(c1, c2) > 1 else "below e0"))
        print(f"  {lbl:<44}{m1:>8.3f}{m2:>9.3f}{min(m1,m2):>8.3f}"
              f"{c1:>8.3f}{c2:>9.3f}{min(c1,c2):>8.3f}{raw:>9.4f}{v:>14}")
    print(f"\n  challenge raw mse is evaluated at the pessimistic pair (Mc = min, eps = max) "
          f"with exact.py's board constants; the incumbent's own row is {base:.4f}.")


if __name__ == "__main__":
    main()
