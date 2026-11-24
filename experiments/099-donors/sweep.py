"""099 arm 2 -- how much of a candidate donor helps, if any.

`arms.py` gives every added donor an equal (1/3) share. That confounds "the new donor is bad"
with "diluting the incumbent is bad" -- which matters a great deal on `orion`, whose reference
IS HCT116 and which therefore prefers its same-atlas relative `orion_hek293t` (single-donor
hek293t reads 1.38x base2 at gamma 0).

Two designs that separate them:

  ADD    {k562:1, hek293t:1, X:w}   for w in 0 .. 1     -- if ANY w > 0 beats w = 0, the donor
                                                           carries something; if the optimum is
                                                           at w = 0 the negative is clean.
  SWAP   {hek293t:1, X:1} against {hek293t:1, k562:1}   -- is X a better second voter than the
                                                           incumbent K562, at equal weight?
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT

CAND = ("kolf2", "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "perturbmatch", "feng_ipsc",
        "orion_hct116", "h1_2025")
WS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--gamma", type=float, default=0.0)
    a = ap.parse_args()
    C = Ctx(a.harness)
    cands = [c for c in CAND if c not in C.banned]
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    union = tuple(dict.fromkeys(base + tuple(cands)))
    C.set_donors(union)
    cov = C.coverage()
    print(f"# {a.harness}, gamma={a.gamma}, base = {base}, P={C.P}")
    print("# coverage: " + ", ".join(f"{n} {cov[n]}/{C.P}" for n in union))

    def run(w):
        dh = C.dh_exact(C.effects(lambda _p: dict(w), gamma=a.gamma))
        X, A, W = C.XAW(dh)
        return dict(X=X, A=A, W=W, kappa=C.kappa(dh), pds=C.pds(dh),
                    nmae=C.nmae(dh), signdot=C.signdot(dh))

    w0 = {n: (1.0 if n in base else 0.0) for n in union}
    r0 = run(w0)
    print(f"\n# ADD:  base2 = {{{base[0]}:1, {base[1]}:1}} plus X at weight w")
    print(f"  {'candidate':<15}{'cov':>6}" + "".join(f"{'w='+str(w):>9}" for w in WS)
          + f"   {'argmax w':>9}")
    out = {"base": r0, "add": {}, "swap": {}, "gamma": a.gamma, "harness": a.harness,
           "coverage": cov, "base_donors": list(base)}
    for c in cands:
        ms, rows = [], {}
        for w in WS:
            ww = dict(w0); ww[c] = w
            r = run(ww)
            rows[w] = r
            ms.append(np.sign(r["X"]) * np.sqrt(r["W"] / r0["W"]))
        out["add"][c] = {str(k): v for k, v in rows.items()}
        print(f"  {c:<15}{cov[c]:>6}" + "".join(f"{m:>9.3f}" for m in ms)
              + f"   {WS[int(np.argmax(ms))]:>9}")

    print(f"\n# SWAP: {{{base[1]}:1, X:1}} -- is X a better partner than {base[0]}?")
    inc = {n: (1.0 if n == base[1] else 0.0) for n in union}
    print(f"  {'partner':<15}{'M':>9}{'kappa':>9}{'pds':>9}{'nmae':>9}")
    for c in [base[0]] + list(cands):
        ww = dict(inc); ww[c] = 1.0
        r = run(ww)
        out["swap"][c] = r
        print(f"  {c:<15}{np.sign(r['X'])*np.sqrt(r['W']/r0['W']):>9.3f}{r['kappa']:>9.5f}"
              f"{r['pds']:>9.5f}{r['nmae']:>9.5f}")
    print(f"\n  {'(base2 reference)':<15}{1.0:>9.3f}{r0['kappa']:>9.5f}"
          f"{r0['pds']:>9.5f}{r0['nmae']:>9.5f}")
    p = ROOT / f"experiments/099-donors/work/sweep_{a.harness}_g{a.gamma}.json"
    json.dump(out, open(p, "w"), indent=1)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
