"""227 — PLOLO replication: the ceiling ladder and the separability gate on the real 300-panel.

Four folds, three cancer lines pooled and `kolf2` reported SEPARATELY (it is pluripotent and
inflates).  Pool construction is 188's validated one: `|lfc| * sqrt(cpm)` -- the Wilcoxon-power
proxy, which reproduced the reference's own BH set at 0.7219 overlap on h1-test -- taken at
h1-test's median gate size of 89 genes per row.  188 sec 0 established that a top-N-by-|lfc| pool
INVERTS `nmae`'s sign, so the pool is not a free choice and this is the one that transfers.

Endpoints on every arm: `nmae`, `pds` (162's exact rank statistic), and 197's energy-weighted
aggregate cosine `cos_agg = sqrt(1 - sse)`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E188 = os.path.join(REPO, "experiments/188-magnitude-calibration")
sys.path.insert(0, E188)
sys.path.insert(0, os.path.join(REPO, "experiments/162-retrodicting-validator"))
from budget188 import FOLDS, deliver, load     # noqa: E402
from mapfit import to_lfc                      # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
POOL_N = 89
RNG = np.random.default_rng(227_004)


def topn(key, n):
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        k = min(n, key.shape[1])
        out[i, np.argpartition(-key[i], k - 1)[:k]] = True
    return out


def endpoints(P, T):
    """162's `pds` + 197's aggregate cosine, on the full scored matrix."""
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-12)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S)
    pds = float(np.mean(((S < own[:, None]).sum(1)
                         + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)))
    a, x, d = float((P * P).sum()), float((P * T).sum()), float((T * T).sum())
    return pds, float(np.sqrt(max(0.0, x * x / (a * d)))) * (1 if x >= 0 else -1)


def nmae_rows(P, T, pool):
    v = []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        den = np.abs(T[i, s]).sum()
        if s.size < 10 or not np.isfinite(den) or den == 0:
            continue
        v.append(np.abs(P[i, s] - T[i, s]).sum() / den)
    return np.array(v)


def qw(P, T, pool):
    v = []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        w = np.abs(T[i, s])
        if s.size < 10 or w.sum() == 0:
            continue
        v.append(w[np.sign(P[i, s]) == np.sign(T[i, s])].sum() / w.sum())
    return float(np.mean(v))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    z = np.load(os.path.join(E188, "work/chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    bA, gA = z["bA"], z["gA"]
    res = {}
    for held in FOLDS:
        F = load(held)
        mk = np.zeros(cg.size, bool)
        mk[F["scored"]] = True
        c = np.flatnonzero(mk & cg)
        bh = F["b"]
        cpmf = 1e6 * np.expm1(bh[c]) / float(np.expm1(bh).sum())
        T = to_lfc(F["truth"], bh, gA)[:, c]
        P0 = to_lfc(deliver(F["pred"], bA, gA), bA, gA)[:, c]
        pool = topn(np.abs(T) * np.sqrt(cpmf)[None, :], POOL_N)
        base = float(nmae_rows(P0, T, pool).mean())
        pds0, cos0 = endpoints(P0, T)
        rec = dict(n_rows=int(T.shape[0]), n_cols=int(c.size), nmae=base, pds=pds0,
                   cos_agg=cos0, q_w=qw(P0, T, pool), arms={})

        def add(name, P):
            nm = float(nmae_rows(P, T, pool).mean())
            p, cs = endpoints(P, T)
            rec["arms"][name] = dict(nmae=nm, d_nmae=base - nm, pds=p, d_pds=p - pds0,
                                     cos_agg=cs, d_cos_pct=100 * (cs / cos0 - 1) if cos0 else 0.0,
                                     q_w=qw(P, T, pool))

        # --- the two channels Phase 1/2 found, replicated -----------------------------------
        ss = np.linspace(0.2, 3.0, 141)
        add("best global scalar",
            P0 * float(ss[int(np.argmin([nmae_rows(P0 * x, T, pool).mean() for x in ss]))]))
        d = np.array([np.median(T[i, np.flatnonzero(pool[i])] - P0[i, np.flatnonzero(pool[i])])
                      for i in range(T.shape[0])])
        add("ORACLE per-row offset (pool)", P0 + d[:, None])
        # realisable estimator of that offset, from our own profile only
        mh = np.array([np.median(P0[i]) for i in range(T.shape[0])])
        add("our own median shift x1", P0 + mh[:, None])
        pr = RNG.permutation(T.shape[0])
        add("SHUFFLE oracle offset", P0 + d[pr][:, None])

        # --- the separability gate: cheapest-in-cosine sign fixes (oracle ordering) ---------
        for f in (0.05, 0.10, 0.20, 0.40):
            P = P0.copy()
            for i in range(T.shape[0]):
                s = np.flatnonzero(pool[i])
                bad = s[np.sign(P0[i, s]) != np.sign(T[i, s])]
                if bad.size == 0:
                    continue
                order = bad[np.argsort(np.abs(P0[i, bad]), kind="stable")]
                k = int(round(f * bad.size))
                if k:
                    P[i, order[:k]] = -P0[i, order[:k]]
            add(f"fix lo-|h| wrong signs, f={f:g}", P)

        res[held] = rec
        print(f"[{held:14s}] rows {rec['n_rows']:3d} cols {rec['n_cols']:5d}  "
              f"nmae {base:.6f}  pds {pds0:.5f}  cos_agg {cos0:.6f}  q_w {rec['q_w']:.4f}",
              flush=True)

    def macro(folds, path):
        return float(np.mean([_get(res[f], path) for f in folds]))

    def _get(d, path):
        for p in path:
            d = d[p]
        return d

    groups = {"3 cancer folds": CANCER, "kolf2 (separate)": ("kolf2",)}
    summary = {}
    for gname, folds in groups.items():
        s = dict(nmae=macro(folds, ["nmae"]), pds=macro(folds, ["pds"]),
                 cos_agg=macro(folds, ["cos_agg"]), q_w=macro(folds, ["q_w"]), arms={})
        for a in res[folds[0]]["arms"]:
            s["arms"][a] = {k: macro(folds, ["arms", a, k])
                            for k in ("nmae", "d_nmae", "pds", "d_pds", "cos_agg",
                                      "d_cos_pct", "q_w")}
        summary[gname] = s
    out = dict(pool=f"|lfc|*sqrt(cpm) top-{POOL_N}", per_fold=res, summary=summary)
    json.dump(out, open(os.path.join(WORK, "plolo227.json"), "w"), indent=1)

    for gname, s in summary.items():
        print(f"\n== {gname} ==  base nmae {s['nmae']:.6f}  pds {s['pds']:.5f}  "
              f"cos_agg {s['cos_agg']:.6f}  q_w {s['q_w']:.4f}")
        print(f"  {'arm':32s} {'nmae':>9s} {'d nmae':>9s} {'q_w':>7s} {'d pds':>9s} "
              f"{'cos_agg':>9s} {'d cos %':>9s}")
        for a, d in s["arms"].items():
            print(f"  {a:32s} {d['nmae']:9.6f} {d['d_nmae']:+9.6f} {d['q_w']:7.4f} "
                  f"{d['d_pds']:+9.6f} {d['cos_agg']:9.6f} {d['d_cos_pct']:+9.3f}")


if __name__ == "__main__":
    main()
