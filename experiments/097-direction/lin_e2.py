"""e2 -- the SHIPPABLE Ahlmann-Eltze variant: W fitted entirely BETWEEN donors.

    embed p with the top-k PCA scores of donor A's effect matrix
    -> ridge-solve W (+ optional intercept) to predict donor B's gene-specific effects on the
       shared perturbations, on the harness gene axis
    -> apply that same W to the harness's targets and score kappa on the harness

Nothing in the fit touches the harness reference.  This is the only arm of experiment 097's
linear route that could actually ship, so it is the one that decides the question.

Also reported, as controls on the same rows:
  * `A raw`      : dh_p = s_A[p] straight through (no fit, no whitener) -- e0's row.
  * `A rank-k`   : dh_p = s_A[p] projected onto donor A's own top-k PCA basis.  This isolates
                   how much of e2's behaviour is the RANK CONSTRAINT rather than the transfer.

Every (k, lambda, intercept) combination for one (A, B) pair is served from three Gram
quantities -- `Gee = Ec'Ec`, `Gey = Ec'Yc`, `||Yc||^2` -- accumulated in float64 over row
chunks, so donor B's 8,000 x 18,000 block is touched twice, not 88 times.

LEAK AUDIT
  * embedding donor A : `data/donors/<A>/delta.npy` only.
  * target donor B    : `data/donors/<B>/delta.npy` only, on the A-cap-B perturbations.
  * harness reference : enters ONLY through the scorer (X, sumD).  Never in the fit.
  * banned per harness: `h1_2025` on h1-test, `orion_hct116` on orion (each IS its harness).
"""
from __future__ import annotations

import itertools
import json

import numpy as np

import lin_kit as K

KS = (5, 10, 20, 50)
LAMS = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)
POOL = ("k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential", "h1_2025",
        "orion_hct116")
KMAX = max(KS)
CHUNK = 1024


def build(t):
    P = t["dh"].shape[0]
    D = {}
    for dn in POOL:
        if dn == t["banned"]:
            continue
        blk = K.donor_block(dn, t["genes"])
        US, V, evr, tot = K.embed_pca(blk["S"], KMAX)
        idx = blk["index"]
        E = np.zeros((P, KMAX))
        cov = np.zeros(P, bool)
        for r, x in enumerate(t["targets"]):
            j = idx.get(str(x))
            if j is not None:
                E[r] = US[j]
                cov[r] = True
        D[dn] = dict(blk=blk, US=US, V=V, evr=evr, E=E, cov=cov)
        print(f"    [{dn}] {blk['S'].shape[0]} perts x {blk['n_genes_hit']}/{t['genes'].size} "
              f"matched genes | evr cum @5 {evr[:5].sum():.4f} @10 {evr[:10].sum():.4f} "
              f"@50 {evr[:50].sum():.4f} | covers {cov.sum()}/{P} harness targets", flush=True)
    return D


def grams(EA, Yrows, icept):
    """`(Gee[k,k], Gey[k,G], sy2, ebar, ybar, ss_tot)` in float64, two streamed passes."""
    n, kk = EA.shape
    G = Yrows.shape[1]
    ybar = np.zeros(G)
    ss_tot = 0.0
    for lo in range(0, n, CHUNK):
        Y = Yrows[lo:lo + CHUNK].astype(np.float64)
        ybar += Y.sum(0)
        ss_tot += float(np.einsum("ij,ij->", Y, Y))
    ybar /= n
    ebar = EA.mean(0) if icept else np.zeros(kk)
    if not icept:
        ybar = np.zeros(G)
    Ec = EA - ebar
    Gee = Ec.T @ Ec
    Gey = np.zeros((kk, G))
    sy2 = 0.0
    for lo in range(0, n, CHUNK):
        Y = Yrows[lo:lo + CHUNK].astype(np.float64) - ybar
        Gey += Ec[lo:lo + CHUNK].T @ Y
        sy2 += float(np.einsum("ij,ij->", Y, Y))
    return Gee, Gey, sy2, ebar, ybar, ss_tot


def solve(Gee, Gey, sy2, ebar, ybar, ss_tot, k, lam_mult):
    """Ridge at rank `k`; returns `(W[k,G], c[G], r2_fit, |W|)` with no second pass over Y."""
    Ge = Gee[:k, :k]
    scale = float(np.trace(Ge)) / k
    W = np.linalg.solve(Ge + lam_mult * max(scale, 1e-30) * np.eye(k), Gey[:k])
    c = ybar - ebar[:k] @ W
    resid = sy2 - 2 * float(np.einsum("kg,kg->", W, Gey[:k])) + float(
        np.einsum("kg,kj,jg->", W, Ge, W))
    return W, c, 1.0 - max(resid, 0.0) / max(ss_tot, 1e-30), float(np.linalg.norm(W))


def target_matrix(t, D, dn):
    P = t["dh"].shape[0]
    blk = D[dn]["blk"]
    M = np.zeros((P, t["genes"].size))
    for r, x in enumerate(t["targets"]):
        j = blk["index"].get(str(x))
        if j is not None:
            M[r] = blk["S"][j]
    return M


def _pick(r):
    return dict(kappa=r["kappa"], eps=r["eps"], cos=r["cos"], Bud=r["B"], mse_min=r["mse_min"])


def run(name):
    t = K.prep(name)
    P = t["dh"].shape[0]
    k0 = K.shape(t["dh"], t)[0]
    print(f"\n########## e2 SHIPPABLE BETWEEN-DONOR FIT -- {name} (P={P}) "
          f"incumbent kappa={k0:.5f} ##########", flush=True)
    D = build(t)
    names = list(D)
    rows = []

    for dn in names:
        if D[dn]["cov"].sum() == 0:
            continue
        M = target_matrix(t, D, dn)
        r = K.shape_full(M, t)
        rows.append(dict(kind="A raw", A=dn, B="-", k=0, lam=0.0, icept=False,
                         cov=int(D[dn]["cov"].sum()), nfit=0, **_pick(r),
                         M_ratio=r["kappa"] / k0, r2_fit=float("nan"), wnorm=float("nan")))
        for k in KS:
            V = D[dn]["V"][:, :k]
            r = K.shape_full((M @ V) @ V.T, t)
            rows.append(dict(kind="A rank-k", A=dn, B="-", k=k, lam=0.0, icept=False,
                             cov=int(D[dn]["cov"].sum()), nfit=0, **_pick(r),
                             M_ratio=r["kappa"] / k0, r2_fit=float("nan"),
                             wnorm=float("nan")))

    for A, B in itertools.permutations(names, 2):
        if D[A]["cov"].sum() < 5:
            continue
        ia, ib = D[A]["blk"]["index"], D[B]["blk"]["index"]
        shared = [p for p in ia if p in ib]
        if len(shared) < 50:
            print(f"    SKIP {A} -> {B}: only {len(shared)} shared perturbations", flush=True)
            continue
        ra = np.array([ia[p] for p in shared])
        rb = np.array([ib[p] for p in shared])
        EA = D[A]["US"][ra]
        YB = D[B]["blk"]["S"][rb]
        for icept in (True, False):
            g = grams(EA, YB, icept)
            for k in KS:
                for lm in LAMS:
                    W, c, r2f, wn = solve(*g, k, lm)
                    DH = D[A]["E"][:, :k] @ W + c
                    DH[~D[A]["cov"]] = 0.0
                    r = K.shape_full(DH, t)
                    rows.append(dict(kind="A->B", A=A, B=B, k=k, lam=lm, icept=icept,
                                     cov=int(D[A]["cov"].sum()), nfit=len(shared), **_pick(r),
                                     M_ratio=r["kappa"] / k0, r2_fit=r2f, wnorm=wn))
        print(f"    {A:>15} -> {B:<15} n_shared={len(shared)}", flush=True)

    for A in names:
        if D[A]["cov"].sum() < 5:
            continue
        for B, C in itertools.combinations([x for x in names if x != A], 2):
            ia = D[A]["blk"]["index"]
            ib, ic = D[B]["blk"]["index"], D[C]["blk"]["index"]
            shared = [p for p in ia if p in ib and p in ic]
            if len(shared) < 50:
                continue
            ra = np.array([ia[p] for p in shared])
            YB = D[B]["blk"]["S"][np.array([ib[p] for p in shared])].copy()   # float32
            YB += D[C]["blk"]["S"][np.array([ic[p] for p in shared])]
            YB *= 0.5
            g = grams(D[A]["US"][ra], YB, False)
            for k in (10, 50):
                for lm in (1e-4, 1e-2, 1e0, 1e2):
                    W, c, r2f, wn = solve(*g, k, lm)
                    DH = D[A]["E"][:, :k] @ W + c
                    DH[~D[A]["cov"]] = 0.0
                    r = K.shape_full(DH, t)
                    rows.append(dict(kind="A->mean(B,C)", A=A, B=f"{B}+{C}", k=k, lam=lm,
                                     icept=False, cov=int(D[A]["cov"].sum()),
                                     nfit=len(shared), **_pick(r), M_ratio=r["kappa"] / k0,
                                     r2_fit=r2f, wnorm=wn))
            del YB
    for r in rows:
        r["harness"] = name
        r["inc_kappa"] = k0
    D.clear()
    K._CACHE.clear()
    return rows


def show(rows, name, top=30):
    print(f"\n  ---- {name}: e2 arms ranked by kappa (top {top}) ----")
    print(f"  {'kind':<14}{'A':<15}{'B':<28}{'k':>4}{'lam':>8}{'ic':>4}{'cov':>5}"
          f"{'nfit':>6}{'kappa':>9}{'M':>7}{'eps':>9}{'r2fit':>9}")
    ok = [r for r in rows if np.isfinite(r["kappa"])]
    for r in sorted(ok, key=lambda x: -x["kappa"])[:top]:
        print(f"  {r['kind']:<14}{r['A']:<15}{r['B']:<28}{r['k']:>4}{r['lam']:>8.0e}"
              f"{str(r['icept'])[0]:>4}{r['cov']:>5}{r['nfit']:>6}{r['kappa']:>9.5f}"
              f"{r['M_ratio']:>7.3f}{r['eps']:>+9.5f}{r['r2_fit']:>9.4f}")


if __name__ == "__main__":
    allr = []
    for nm in ("h1-test", "orion"):
        rr = run(nm)
        show(rr, nm)
        allr += rr
    with open("work/lin_e2.json", "w") as f:
        json.dump(allr, f, default=float)
    print(f"\n  wrote work/lin_e2.json ({len(allr)} rows)")
