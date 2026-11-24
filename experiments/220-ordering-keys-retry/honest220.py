"""220 — the honest estimate, and the per-perturbation form of the 0.90 bar.

TWO THINGS 189 AND 190 DID NOT DO.

1  LEAVE-ONE-TARGET-OUT lambda.  189 selected lambda on its evaluation surface and flagged
   it; 190 took a max over 13 lambdas x 8 networks.  Here lambda is chosen on the other 81
   targets and scored on the held-out one, which yields a per-target held-out delta vector
   and therefore a paired bootstrap se on the SAME axis 189 quoted.  This is the number the
   two experiments should have reported.

2  THE BAR PER PERTURBATION.  `reach`'s 0.90 purity bar is a ceiling function applied to each
   perturbation separately, so "does any key reach it for a non-trivial share of
   perturbations" is a count, not an average.  Counted here at depth 5, 10, 20 and 50.

    uv run python honest220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))

LAM189 = (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0)
LAM190 = (-4.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)


def main() -> None:
    out = {}
    old = C.Instrument(C.DE_A15, "a15")
    new = C.Instrument(C.DE_RC, "rc_callset")
    C.check_old(old)

    import nets190
    print("\n" + "=" * 100)
    print("1. LEAVE-ONE-TARGET-OUT lambda — the honest version of 189's +0.00408 and 190's "
          "+0.00137")
    print("=" * 100)
    res = {}
    for T, tag in ((old, "a15"), (new, "rc_callset")):
        base = T.reach_vec(T.key_absH())
        N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
        R = C.agreement_R(N, T.H)
        fams = {
            "189": np.vstack([T.reach_vec(C.key189(T.H, T.theta, l)) - base
                              for l in LAM189]),
            "190": np.vstack([T.reach_vec(np.abs(T.H) * np.exp(np.clip(l * R, -50, 50)))
                              - base for l in LAM190]),
        }
        print(f"\n  --- {tag} ---")
        for nm, M in fams.items():
            grid = LAM189 if nm == "189" else LAM190
            n = M.shape[1]
            loo = np.empty(n)
            picked = []
            for i in range(n):
                keep = np.ones(n, bool)
                keep[i] = False
                j = int(np.argmax(M[:, keep].mean(1)))
                picked.append(grid[j])
                loo[i] = M[j, i]
            bb = C.paired_bootstrap(loo, np.zeros(n))
            insample = float(M.mean(1).max())
            uniq = {l: picked.count(l) for l in sorted(set(picked))}
            print(f"  {nm}: IN-SAMPLE max {insample:+.5f}   |   LEAVE-ONE-OUT "
                  f"{loo.mean():+.5f}  target-bootstrap se {bb['se']:.5f}  "
                  f"95% CI [{bb['lo']:+.5f}, {bb['hi']:+.5f}]  P(>0) {bb['p_better']:.3f}")
            print(f"       lambda chosen per fold: {uniq}   shrinkage "
                  f"{insample/max(loo.mean(),1e-9):.1f}x")
            res[f"{tag}|{nm}"] = dict(in_sample=insample, loo=float(loo.mean()),
                                      se=bb["se"], lo=bb["lo"], hi=bb["hi"],
                                      p=bb["p_better"],
                                      lam_counts={str(k): v for k, v in uniq.items()})
    out["loo"] = res

    print("\n" + "=" * 100)
    print("2. THE 0.90 BAR PER PERTURBATION — a COUNT, not an average")
    print("=" * 100)
    T = new
    N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
    R = C.agreement_R(N, T.H)
    keys = {"|H| (incumbent)": T.key_absH(),
            "delivered (production)": None,
            "189 lam 0.50": C.key189(T.H, T.theta, 0.5),
            "190 lam -0.25": np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50)),
            "ORACLE": T.key_oracle()}
    depths = (5, 10, 20, 50)
    print(f"  of {int(T.scored_mask().sum())} scored targets, the number clearing 0.90 "
          f"purity at each depth:")
    print(f"{'key':<24}" + "".join(f"{'d='+str(d):>10}" for d in depths) +
          f"{'median k*':>12}{'reach':>10}")
    cnt = {}
    for nm, k in keys.items():
        fa = T.front_acc(k, depths=depths)
        ks = T.kstar_vec(k)
        c = {d: int(np.nansum(fa[d] >= 0.9 - 1e-9)) for d in depths}
        cnt[nm] = dict(counts=c, kstar_med=float(np.median(ks)),
                       reach=float(T.reach_vec(k).mean()))
        print(f"{nm:<24}" + "".join(f"{c[d]:>10d}" for d in depths) +
              f"{np.median(ks):>12.1f}{T.reach_vec(k).mean():>10.5f}")
    out["bar_counts"] = cnt
    print("\n  (a perturbation clears d = 5 only with 5/5 correct, d = 10 with 9/10, "
          "d = 20 with 18/20.)")
    print("  The two keys move these counts by at most one perturbation. The oracle moves "
          "them by 40-60.")

    print("\n" + "=" * 100)
    print("3. HOW MANY PERTURBATIONS CARRY THE EFFECT — a jackknife on the target axis")
    print("=" * 100)
    base = T.reach_vec(T.key_absH())
    jk = {}
    for nm, k in (("189 lam 0.50", C.key189(T.H, T.theta, 0.5)),
                  ("190 lam -0.25", np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50)))):
        d = T.reach_vec(k) - base
        o = np.argsort(-d)
        n = d.size
        drop1 = float(np.delete(d, o[0]).mean())
        drop3 = float(np.delete(d, o[:3]).mean())
        jk[nm] = dict(full=float(d.mean()), drop1=drop1, drop3=drop3,
                      top1_share=float(d[o[0]] / n / max(d.mean(), 1e-12)),
                      n_moved=int((np.abs(d) > 1e-12).sum()))
        print(f"  {nm:<16} full {d.mean():+.5f}   drop the single best target "
              f"{drop1:+.5f}   drop the best three {drop3:+.5f}")
        print(f"       the single largest target contributes "
              f"{100*d[o[0]]/n/max(d.mean(),1e-12):.0f}% of the whole effect; "
              f"{int((np.abs(d) > 1e-12).sum())} of {n} targets move at all")
    out["jackknife"] = jk

    json.dump(out, open(C.WORK / "honest220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'honest220.json'}")


if __name__ == "__main__":
    main()
