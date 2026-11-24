"""Where does e2's apparent win live?  Split every arm into its CONSTANT and CENTRED parts.

`_across_pert_budget` is invariant to adding a perturbation-constant vector, so an arm can raise
kappa = X/sqrt(B*sumD) purely by bolting on a constant that happens to align with the panel's
average reference effect -- B does not move, X does.  `eps = (A-B)/B` is the only thing that
sees it, and it is the term the challenge panel multiplies by `claim = 1.079 * sumD`.

So for each arm:
    c   = sum_p w_p DH_p / sum_p w_p      (the perturbation-constant part)
    DHc = DH - c                          (the differential part, the only part with a budget)
and kappa / eps / mse_min are reported for DH, for c alone, and for DHc alone -- on both
harnesses.  An arm whose differential part flips sign between harnesses is dead under rule 1
regardless of what its constant does.

Also: the permutation null repeated WITHOUT the intercept, which is the version of the null
that can actually reach zero.
"""
from __future__ import annotations

import numpy as np

import lin_kit as K
from lin_null import WIN, donor_S_on_targets, fit_arm

NPERM = 20


def split(DH, t):
    w = t["w"]
    c = (w[:, None] * DH).sum(0) / w.sum()
    return c, DH - c[None, :]


def line(lbl, DH, t, k0, corr):
    r = K.shape_full(DH, t)
    kap = r["kappa"]
    print(f"    {lbl:<50}{kap:>10.5f}{kap/k0:>8.3f}{kap/(k0*corr):>8.3f}"
          f"{r['eps']:>+11.5f}{r['cos']:>9.5f}{r['mse_min']:>9.4f}")
    return r


def main():
    store = {}
    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P = t["dh"].shape[0]
        r0 = K.shape_full(t["dh"], t)
        k0 = r0["kappa"]
        corr = np.sqrt(r0["B"] / (r0["B"] - float((t["Chat"] * t["w"]).sum()) * (1 - 1 / P)))
        print(f"\n################ {name} (P={P}) | incumbent kappa {k0:.5f}, "
              f"emission-corrected {k0*corr:.5f} ################")
        print(f"    {'arm':<50}{'kappa':>10}{'M':>8}{'M/corr':>8}{'eps':>11}"
              f"{'cos':>9}{'mse_min':>9}")

        arms = {}
        line("INCUMBENT (whitened, emitted)", t["dh"], t, k0, corr)

        Mk, ck, medk, _ = donor_S_on_targets(t, "k562_gwps")
        Mh, ch, medh, _ = donor_S_on_targets(t, "orion_hek293t")
        acc = np.zeros_like(Mk)
        cnt = np.zeros(P)
        for M, cv, md in ((Mk, ck, medk), (Mh, ch, medh)):
            acc[cv] += M[cv] / md
            cnt[cv] += 1
        acc[cnt > 0] /= cnt[cnt > 0, None]
        arms["plain renorm-mean k562+hek293t"] = acc
        arms["e2 winner (k=50 ic=T lam=1e0)"] = fit_arm(t, WIN)[0]
        arms["e2 winner, NO intercept"] = fit_arm(t, dict(WIN, icept=False))[0]

        for lbl, DH in arms.items():
            r = line(lbl, DH, t, k0, corr)
            c, DHc = split(DH, t)
            line(f"    ^ constant part only", np.repeat(c[None, :], P, 0), t, k0, corr)
            rc = line(f"    ^ centred (differential) part only", DHc, t, k0, corr)
            store.setdefault(lbl, {})[name] = (r["kappa"] / k0, rc["kappa"] / k0,
                                               r["eps"], rc["eps"])

        print(f"\n    permutation null WITHOUT the intercept ({NPERM} seeds) -- the version "
              f"that can reach zero")
        ba = K.donor_block(WIN["A"], t["genes"])
        bb = K.donor_block(WIN["B"], t["genes"])
        ns = len([p for p in ba["index"] if p in bb["index"]])
        rng = np.random.default_rng(11)
        ks = []
        for _ in range(NPERM):
            _, r, _, _, _ = fit_arm(t, dict(WIN, icept=False), rng.permutation(ns))
            ks.append(r["kappa"])
        ks = np.array(ks)
        real = K.shape_full(arms["e2 winner, NO intercept"], t)["kappa"]
        print(f"      shuffled kappa mean {ks.mean():+.5f} sd {ks.std():.5f} "
              f"[{ks.min():+.5f}, {ks.max():+.5f}]   real {real:+.5f}   "
              f"({int((ks >= real).sum())}/{NPERM} shuffles >= real)")
        K._CACHE.clear()

    print("\n\n===== cross-surface, the differential part only (rule 1: min, and sign) =====")
    print(f"  {'arm':<40}{'M h1':>9}{'M ori':>9}{'M use':>9}"
          f"{'Mc h1':>9}{'Mc ori':>9}{'Mc use':>9}{'verdict':>12}")
    for lbl, d in store.items():
        (m1, c1, e1, ce1), (m2, c2, e2, ce2) = d["h1-test"], d["orion"]
        v = "DEAD (sign)" if np.sign(c1) != np.sign(c2) else (
            "beats e0" if min(c1, c2) > 1 else "below e0")
        print(f"  {lbl:<40}{m1:>9.3f}{m2:>9.3f}{min(m1, m2):>9.3f}"
              f"{c1:>9.3f}{c2:>9.3f}{min(c1, c2):>9.3f}{v:>12}")


if __name__ == "__main__":
    main()
