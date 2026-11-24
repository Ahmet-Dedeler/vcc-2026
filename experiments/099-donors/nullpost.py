"""099 arm 5b -- permutation null for the post-whitening residual arms.

Adding any vector of the right magnitude to a whitened effect matrix changes `W`; the question
is whether the change is target-specific. The null shuffles the candidate's perturbation map,
so its gene content, magnitude, coverage and the exact size of the addition are preserved and
only the correspondence with our own rows is destroyed.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT, ALPHA

CD4 = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--gamma", type=float, default=0.65)
    ap.add_argument("--f", type=float, default=0.5)
    ap.add_argument("--n-null", type=int, default=60)
    ap.add_argument("--cands", default="kolf2,cd4_rest,cd4_stim8hr,cd4_stim48hr,"
                                       "perturbmatch,cd4_mean3_gn")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    C = Ctx(a.harness)
    cands = [c for c in a.cands.split(",") if c]
    need = sorted({x for c in cands for x in ([c] if not c.startswith("cd4_mean") else CD4)})
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    union = tuple(dict.fromkeys(base + tuple(n for n in need if n not in C.banned)))
    m = C.set_donors(union)
    w0 = {n: (1.0 if n in base else 0.0) for n in union}
    E0 = np.vstack([m.effect(p, 1.0, 0.0, weights=dict(w0)) for p in C.perts]).astype(np.float64)
    if a.gamma > 0:
        before = float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1)))
        E0[:, C.wcols] = C.W.apply(E0[:, C.wcols], a.gamma, 0.0, renorm=False)
        E0[:, C.wcols] *= before / (
            float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1))) or 1.0)
    nE = float(np.median(np.linalg.norm(E0, axis=1)))

    def finish(E):
        return C.dh_exact(E * (ALPHA * m.effect_scale /
                               (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)))

    def raw_of(c):
        V = np.zeros_like(E0)
        for i, p in enumerate(C.perts):
            j = m.pert_index[c].get(p)
            if j is None or j >= m.s[c].shape[0]:
                continue
            V[i, m.tcols[c]] = m.s[c][j][m.dcols[c]]
        return V * (nE / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0))

    def vec(c):
        if c == "cd4_mean3_gn":
            V = sum(raw_of(x) for x in CD4) / 3.0
            sc = np.abs(V).mean(0); live = sc[sc > 0]
            sc = np.maximum(sc, np.percentile(live, 25.0))
            b_ = float(np.median(np.linalg.norm(V, axis=1)))
            V = V / sc
            V *= b_ / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0)
        else:
            V = raw_of(c)
        return V * (nE / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0))

    dh0 = finish(E0)
    X0, A0, W0 = C.XAW(dh0)
    print(f"# {a.harness} gamma={a.gamma} f={a.f}  base kappa={C.kappa(dh0):.5f} "
          f"pds={C.pds(dh0):.5f} nmae={C.nmae(dh0):.5f}")
    print(f"  {'candidate':<15}{'obs M':>8}{'p50':>8}{'p95':>8}{'max':>8}{'p':>7}"
          f"{'d_pds':>9}{'d_nmae':>9}")
    rng = np.random.default_rng(a.seed)
    out = {}
    for c in cands:
        dh = finish(E0 + a.f * vec(c))
        X, A, W = C.XAW(dh)
        obs = float(np.sign(X) * np.sqrt(W / W0))
        dp = C.pds(dh) - C.pds(dh0); dn = C.nmae(dh) - C.nmae(dh0)
        parts = CD4 if c.startswith("cd4_mean") else (c,)
        origs = {x: m.pert_index[x] for x in parts}
        nulls = []
        for _ in range(a.n_null):
            for x in parts:
                o = origs[x]; nrow = m.s[x].shape[0]
                ks = [k for k in o if k and o[k] < nrow]
                vs = rng.permutation(np.array([o[k] for k in ks]))
                m.pert_index[x] = {k: int(v) for k, v in zip(ks, vs)}
            dhn = finish(E0 + a.f * vec(c))
            Xn, An, Wn = C.XAW(dhn)
            nulls.append(float(np.sign(Xn) * np.sqrt(Wn / W0)))
        for x in parts:
            m.pert_index[x] = origs[x]
        nulls = np.array(nulls)
        out[c] = dict(obs=obs, p50=float(np.median(nulls)), p95=float(np.percentile(nulls, 95)),
                      mx=float(nulls.max()), p=float((nulls >= obs).mean()),
                      d_pds=float(dp), d_nmae=float(dn))
        print(f"  {c:<15}{obs:>8.3f}{np.median(nulls):>8.3f}{np.percentile(nulls,95):>8.3f}"
              f"{nulls.max():>8.3f}{(nulls>=obs).mean():>7.3f}{dp:>+9.5f}{dn:>+9.5f}")
    p = ROOT / f"experiments/099-donors/work/nullpost_{a.harness}_g{a.gamma}_f{a.f}.json"
    json.dump(dict(harness=a.harness, gamma=a.gamma, f=a.f, out=out), open(p, "w"), indent=1)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
