"""099 arm 1 -- directional accuracy as a function of the DONOR SET, on one harness.

Every arm is the same pipeline with the same frozen whitener, frozen `L`, and the same
amplitude policy; the only difference is which donors vote in `MultiDonorModel.effect()`.

IMPLEMENTATION NOTE. The model is built ONCE over the UNION of every candidate donor and an
arm is expressed as a weight dict that puts 0 on the donors it excludes. That is exactly
equivalent to rebuilding the model on the subset up to a global scalar (`_renorm` differs only
by the common `ref`, `_global` only by a scalar, and the pipeline renormalises the median
effect norm at the end) -- which `--check` verifies against a genuinely rebuilt subset model.
"""
from __future__ import annotations
import argparse, itertools, json, sys, time
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, BASE_DONORS, ROOT

UNION = ("k562_gwps", "orion_hek293t", "orion_hct116", "kolf2", "perturbmatch", "feng_ipsc",
         "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "h1_2025")


def arms_for(harness: str, union: tuple[str, ...]):
    """(label, donor tuple) for every configuration legal on `harness`."""
    banned = {"orion": {"orion_hct116"}, "h1-test": {"h1_2025"},
              "hek293t": {"orion_hek293t"}}[harness]
    ok = [d for d in union if d not in banned]
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in banned)[:2]
    out = [("base2", base)]
    extras = [d for d in ok if d not in base]
    for e in extras:
        out.append((f"base2+{e}", base + (e,)))
    cd4 = tuple(d for d in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr") if d in ok)
    if len(cd4) == 3:
        out.append(("base2+cd4x3", base + cd4))
    if "kolf2" in ok and len(cd4) == 3:
        out.append(("base2+kolf2+cd4x3", base + ("kolf2",) + cd4))
    if "kolf2" in ok and "cd4_rest" in ok:
        out.append(("base2+kolf2+cd4_rest", base + ("kolf2", "cd4_rest")))
    if "orion_hct116" in ok:
        out.append(("PRODUCTION(base2+hct116)", base + ("orion_hct116",)))
        if "kolf2" in ok:
            out.append(("prod3+kolf2", base + ("orion_hct116", "kolf2")))
        if len(cd4) == 3:
            out.append(("prod3+cd4x3", base + ("orion_hct116",) + cd4))
    # single donors, for the "uniform beats the best single donor" check
    for d in base:
        out.append((f"single:{d}", (d,)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--gamma", type=float, default=0.65)
    a = ap.parse_args()

    t0 = time.time()
    C = Ctx(a.harness)
    banned = C.banned
    union = tuple(d for d in UNION if d not in banned)
    C.set_donors(union)
    cov = C.coverage()
    print(f"# harness {a.harness}: P={C.P} targets, {C.rows.size} with a usable reference DE "
          f"gate; banned donor(s) {sorted(banned) or 'none'}", flush=True)
    print(f"# donor coverage of this panel: "
          + ", ".join(f"{n} {cov[n]}/{C.P}" for n in union), flush=True)

    arms = arms_for(a.harness, UNION)
    if a.check:
        # equivalence of union+zero-weights against a genuinely rebuilt subset model
        base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in banned)[:2]
        w = {n: (1.0 if n in base else 0.0) for n in union}
        dh_u = C.dh_exact(C.effects(lambda _p: dict(w), gamma=a.gamma))
        C.set_donors(base)
        dh_s = C.dh_exact(C.effects(gamma=a.gamma))
        print(f"# CHECK union-with-zero-weights vs rebuilt subset: max|dh diff| = "
              f"{np.abs(dh_u - dh_s).max():.3e}, W {C.XAW(dh_u)[2]:.7g} vs {C.XAW(dh_s)[2]:.7g}")
        C.set_donors(union)

    ref = None
    print(f"\n  {'arm':<26}{'nd':>3}{'X':>11}{'A':>10}{'W':>12}{'M':>7}"
          f"{'kappa':>9}{'pds':>9}{'nmae':>9}{'signdot':>9}{'cosmed':>8}")
    res = {}
    for lbl, dons in arms:
        w = {n: (1.0 if n in dons else 0.0) for n in union}
        E = C.effects(lambda _p: dict(w), gamma=a.gamma)
        dh = C.dh_exact(E)
        X, A, W = C.XAW(dh)
        if ref is None:
            ref = W
        r = dict(donors=list(dons), X=X, A=A, W=W, M=float(np.sign(X) * np.sqrt(W / ref)),
                 kappa=C.kappa(dh), pds=C.pds(dh), nmae=C.nmae(dh), signdot=C.signdot(dh),
                 cosmed=float(np.median(C.per_pert_cos(dh))))
        res[lbl] = r
        print(f"  {lbl:<26}{len(dons):>3}{X:>11.5g}{A:>10.5g}{W:>12.6g}{r['M']:>7.3f}"
              f"{r['kappa']:>9.5f}{r['pds']:>9.5f}{r['nmae']:>9.5f}{r['signdot']:>9.5f}"
              f"{r['cosmed']:>8.4f}", flush=True)
    out = ROOT / f"experiments/099-donors/work/arms_{a.harness}_g{a.gamma}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"harness": a.harness, "gamma": a.gamma, "coverage": cov, "arms": res}, open(out, "w"), indent=1)
    print(f"\n# {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
