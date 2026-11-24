"""Cross-surface join for the linear-transfer arms: headline table + the transfer test.

RESEARCH-METHOD rule 1: the transferable estimate is the SMALLER of the two harnesses, and a
sign disagreement kills an arm.  Honesty requirement: a hyper-parameter tuned on one harness
must be reported APPLIED TO THE OTHER; the in-sample argmax is not a result.

Two guards this file exists to enforce:

* **kappa alone is gameable, and heavy ridge games it.**  A perturbation-CONSTANT prediction has
  `B = _across_pert_budget = 0` exactly, so kappa = X/sqrt(B*sumD) diverges while the arm is
  worth nothing.  Any arm whose common-mode share `eps = (A-B)/B` runs away is in that regime.
  Arms are therefore split into a SANE band (`eps <= EPS_SANE`) and a degenerate one, and the
  ranking is by the objective (`challenge raw mse`), never by kappa.
* **the cliff formula is an approximation once eps is large.**  `ch_raw` is exact.py's
  `1 + [Kc + claim*eps - 2 sqrt(claim*sumD) kappa]/sumD`; `ch_exact` re-minimises the true
  piecewise objective over the global amplitude.  They agree to 1e-4 in the sane band.
"""
from __future__ import annotations

import json

import numpy as np

import lin_kit as K

EPS_SANE = 1.0          # eps above this means the common mode exceeds the differential
BASE_EPS = 0.01424      # the incumbent's own eps (h1-test), exact.py's base row


def ch_exact(M, eps):
    """Challenge-panel MSE at the true amplitude optimum, not assuming the rho cliff."""
    B, sD, cl, Kc = K.CH["B"], K.CH["sumD"], K.CH["claim"], K.CH["Kc"]
    X = K.CH["kappa0"] * M * np.sqrt(B * sD)
    return K.mse_min(X, B * (1 + eps), B, sD, cl, Kc)[0]


def join(pairs):
    recs = []
    for key, d in pairs.items():
        if "h1-test" not in d or "orion" not in d:
            continue
        a, b = d["h1-test"], d["orion"]
        m1 = a.get("M", a.get("M_ratio"))
        m2 = b.get("M", b.get("M_ratio"))
        e1, e2 = a["eps"], b["eps"]
        if not all(np.isfinite(x) for x in (m1, m2, e1, e2)):
            continue
        M, E = min(m1, m2), max(e1, e2)
        recs.append(dict(key=key, m1=m1, m2=m2, M=M, e1=e1, e2=e2, E=E,
                         k1=a["kappa"], k2=b["kappa"], a=a, b=b,
                         raw=K.challenge_raw(M, E), exact=ch_exact(M, E),
                         sc=K.challenge_score(M, E),
                         sane=(abs(E) <= EPS_SANE) and (abs(e1) <= EPS_SANE),
                         sign_ok=(np.sign(m1) == np.sign(m2))))
    return recs


def table(recs, label, fmt, top=22, sane_only=True):
    base_raw, base_sc = K.challenge_raw(1.0, BASE_EPS), K.challenge_score(1.0, BASE_EPS)
    pool = [r for r in recs if r["sane"]] if sane_only else recs
    print(f"\n===== {label} =====")
    print(f"  {len(recs)} arms on BOTH harnesses; {len(pool)} in the sane band "
          f"(|eps| <= {EPS_SANE} on both), {len(recs)-len(pool)} degenerate "
          f"(prediction collapsing toward a perturbation-constant, B -> 0)")
    print(f"  INCUMBENT reference: M=1.000 eps=+{BASE_EPS:.5f} -> challenge raw {base_raw:.4f}"
          f"  score {base_sc:.4f}  overall {base_sc/6:.4f}")
    print(f"  {'arm':<50}{'M h1':>8}{'M ori':>9}{'M use':>8}{'eps h1':>9}{'eps ori':>9}"
          f"{'ch raw':>9}{'ch ex':>9}{'d/6':>9}{'sgn':>5}")
    for r in sorted(pool, key=lambda x: x["raw"])[:top]:
        print(f"  {fmt(r['key']):<50}{r['m1']:>8.3f}{r['m2']:>9.3f}{r['M']:>8.3f}"
              f"{r['e1']:>+9.4f}{r['e2']:>+9.4f}{r['raw']:>9.4f}{r['exact']:>9.4f}"
              f"{(r['sc']-base_sc)/6:>+9.4f}{'ok' if r['sign_ok'] else 'FLIP':>5}")
    return pool, base_sc


def transfer(pool, label, fmt, group):
    print(f"\n  ---- {label}: cross-surface hyper-parameter transfer (sane band only) ----")
    groups = {"(all sane arms)": pool}
    for r in pool:
        groups.setdefault(group(r["key"]), []).append(r)
    for gname, rs in groups.items():
        if len(rs) < 2:
            continue
        b1 = max(rs, key=lambda r: r["m1"])
        b2 = max(rs, key=lambda r: r["m2"])
        bb = min(rs, key=lambda r: r["raw"])
        print(f"    {gname}   [{len(rs)} arms, "
              f"{sum(1 for r in rs if not r['sign_ok'])} sign flips, "
              f"{sum(1 for r in rs if r['M'] > 1)} with min-M > 1]")
        print(f"      tuned on h1    : {fmt(b1['key'])}")
        print(f"          M h1 {b1['m1']:+.3f} (in-sample)  ->  M orion {b1['m2']:+.3f} "
              f"(THE NUMBER)   pessimistic ch raw {b1['raw']:.4f}")
        print(f"      tuned on orion : {fmt(b2['key'])}")
        print(f"          M orion {b2['m2']:+.3f} (in-sample)  ->  M h1 {b1['m1']:+.3f}"
              f" / {b2['m1']:+.3f} (THE NUMBER)   pessimistic ch raw {b2['raw']:.4f}")
        print(f"      best pessimistic: {fmt(bb['key'])}  M use {bb['M']:.3f} "
              f"eps use {bb['E']:+.5f}  ch raw {bb['raw']:.4f}")


def main():
    fmt1 = lambda k: f"{k[0]} k={k[1]} ic={str(k[2])[0]} lam={k[3]:.0e}"                  # noqa
    fmt2 = lambda k: f"{k[0]}|{k[1][:11]}->{k[2][:13]}|k={k[3]} ic={str(k[4])[0]} l={k[5]:.0e}"  # noqa

    rows1 = json.load(open("work/lin_e1.json"))
    p1 = {}
    for r in rows1:
        p1.setdefault((r["emb"], r["k"], r["icept"], r["lam"]), {})[r["harness"]] = r
    recs1 = join(p1)
    pool1, base_sc = table(recs1, "e1  SUPERVISED CEILING -- NOT SHIPPABLE", fmt1)
    transfer(pool1, "e1", fmt1, lambda k: f"emb={k[0]} icept={k[2]}")
    print("\n  ---- e1: the lambda curve, no intercept, whole sweep (not just the argmax) ----")
    for emb in ("k562_gwps", "orion_hek293t", "both"):
        for k in (10, 50):
            sub = sorted([r for r in recs1 if r["key"][0] == emb and r["key"][1] == k
                          and r["key"][2] is False], key=lambda r: r["key"][3])
            if not sub:
                continue
            print(f"    {emb} k={k} icept=False")
            print("      " + "  ".join(f"{r['key'][3]:.0e}:{r['M']:+.2f}" for r in sub))

    rows2 = json.load(open("work/lin_e2.json"))
    p2 = {}
    for r in rows2:
        p2.setdefault((r["kind"], r["A"], r["B"], r["k"], r["icept"], r["lam"]),
                      {})[r["harness"]] = r
    recs2 = join(p2)
    pool2, _ = table(recs2, "e2  SHIPPABLE BETWEEN-DONOR FIT -- the arm that decides it", fmt2,
                     top=25)
    transfer(pool2, "e2", fmt2, lambda k: f"kind={k[0]} icept={k[4]}")

    print("\n===== e3  VERDICT: does e2 beat e0 on BOTH harnesses? =====")
    win = [r for r in pool2 if r["M"] > 1.0 and r["sign_ok"]]
    best = max(pool2, key=lambda r: r["M"])
    bestobj = min(pool2, key=lambda r: r["raw"])
    print(f"  sane e2 arms with min(M) > 1 and no sign flip: {len(win)} / {len(pool2)}")
    print(f"  best min-M arm      : {fmt2(best['key'])}")
    print(f"      M = {best['M']:.3f}  (h1 {best['m1']:.3f}, orion {best['m2']:.3f})  "
          f"eps {best['e1']:+.4f}/{best['e2']:+.4f}  ch raw {best['raw']:.4f}")
    print(f"  best-objective arm  : {fmt2(bestobj['key'])}")
    print(f"      M = {bestobj['M']:.3f}  eps use {bestobj['E']:+.5f}  "
          f"ch raw {bestobj['raw']:.4f} vs incumbent "
          f"{K.challenge_raw(1.0, BASE_EPS):.4f}  "
          f"overall delta {(bestobj['sc']-base_sc)/6:+.4f}")
    if best["M"] < 1.0:
        print(f"  => e2 LOSES.  Shortfall on the transferable surface: "
              f"{(1-best['M'])*100:.1f}% of the incumbent's alignment.")


if __name__ == "__main__":
    main()
