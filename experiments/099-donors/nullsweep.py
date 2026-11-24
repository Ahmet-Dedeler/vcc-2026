"""099 arm 3 -- is a candidate donor's marginal gain INFORMATION, or just dilution?

`sweep.py` shows `base2 + X` at w = 0.25 reading 1.01-1.03x on orion. But on orion the
incumbent `k562_gwps` alone reads 0.169x, so ANY third vector that displaces k562 weight raises
`W` for free. The null that separates the two is a **row permutation of the candidate**: keep
the donor's magnitude, gene axis, coverage pattern and mixture share exactly, and destroy only
the correspondence between its perturbation and ours.

  observed  = M(base2 + X at w)
  null      = M(base2 + shuffle(X) at w),  n draws

If the observed value sits inside the null, the donor contributes dilution, not direction.

Also reports every candidate ALONE, which is the most direct statement of whether it predicts
this reference at all.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT

CAND = ("kolf2", "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "perturbmatch", "feng_ipsc",
        "orion_hct116", "h1_2025")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--gamma", type=float, default=0.0)
    ap.add_argument("--w", type=float, default=0.25)
    ap.add_argument("--n-null", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    C = Ctx(a.harness)
    cands = [c for c in CAND if c not in C.banned]
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    union = tuple(dict.fromkeys(base + tuple(cands)))
    m = C.set_donors(union)
    cov = C.coverage()
    print(f"# {a.harness}, gamma={a.gamma}, mixture weight w={a.w}, {a.n_null} permutations")
    print(f"# base = {base}")

    def M(w, ref):
        dh = C.dh_exact(C.effects(lambda _p: dict(w), gamma=a.gamma))
        X, _, W = C.XAW(dh)
        return np.sign(X) * np.sqrt(W / ref)

    w0 = {n: (1.0 if n in base else 0.0) for n in union}
    dh0 = C.dh_exact(C.effects(lambda _p: dict(w0), gamma=a.gamma))
    ref = C.XAW(dh0)[2]

    print(f"\n# ALONE -- does the donor predict this reference by itself?")
    print(f"  {'donor':<15}{'cov':>6}{'M_alone':>10}{'kappa':>9}{'pds':>9}")
    alone = {}
    for c in list(base) + cands:
        ww = {n: (1.0 if n == c else 0.0) for n in union}
        dh = C.dh_exact(C.effects(lambda _p: dict(ww), gamma=a.gamma))
        _x,_a,_w = C.XAW(dh)
        alone[c] = dict(M=float(np.sign(_x)*np.sqrt(_w / ref)), kappa=C.kappa(dh), pds=C.pds(dh))
        print(f"  {c:<15}{cov[c]:>6}{alone[c]['M']:>10.3f}{alone[c]['kappa']:>9.5f}"
              f"{alone[c]['pds']:>9.5f}")

    print(f"\n# PERMUTATION NULL on the mixture, w={a.w}")
    print(f"  {'candidate':<15}{'observed':>10}{'null p50':>10}{'null p95':>10}"
          f"{'null max':>10}{'p':>8}")
    res = {}
    rng = np.random.default_rng(a.seed)
    for c in cands:
        ww = dict(w0); ww[c] = a.w
        obs = M(ww, ref)
        # permute the candidate's perturbation index map, keeping its coverage pattern
        orig = m.pert_index[c]
        nrow = m.s[c].shape[0]
        # `Donor.load` splits perts.txt on newlines, so a trailing blank maps to an
        # out-of-range row that is never queried in normal use.  Drop it before permuting.
        keys = [k for k in orig if k and orig[k] < nrow]
        vals = np.array([orig[k] for k in keys])
        nulls = []
        for _ in range(a.n_null):
            perm = rng.permutation(vals)
            m.pert_index[c] = {k: int(v) for k, v in zip(keys, perm)}
            nulls.append(M(ww, ref))
        m.pert_index[c] = orig
        nulls = np.array(nulls)
        p = float((nulls >= obs).mean())
        res[c] = dict(obs=float(obs), null_p50=float(np.median(nulls)),
                      null_p95=float(np.percentile(nulls, 95)),
                      null_max=float(nulls.max()), p=p, alone=alone[c])
        print(f"  {c:<15}{obs:>10.4f}{np.median(nulls):>10.4f}"
              f"{np.percentile(nulls,95):>10.4f}{nulls.max():>10.4f}{p:>8.3f}")
    out = ROOT / f"experiments/099-donors/work/null_{a.harness}_g{a.gamma}_w{a.w}.json"
    json.dump(dict(harness=a.harness, gamma=a.gamma, w=a.w, n_null=a.n_null,
                   base=list(base), coverage=cov, alone=alone, null=res), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
