"""e1 -- SUPERVISED CEILING for the Ahlmann-Eltze linear form.  **NOT SHIPPABLE.**

Construction (Ahlmann-Eltze, Huber & Anders 2025, "the linear model with a cross-cell-line
pretrained perturbation embedding"):

    embed p with the top-k PCA scores of a DONOR's perturbation x gene effect matrix
    -> ridge-solve W (plus an optional control-centred intercept) against an L2 objective
    -> predict dh_p = e_p @ W + c

Here W is fitted on a random HALF of the harness's own perturbations against the harness's own
reference `d_p`, and kappa is read on the HELD-OUT half.  That is deliberately a leak with
respect to shipping -- the challenge contexts carry zero perturbation data, so there is nothing
to fit W on there.  It answers only: *is the functional form better than ours when it is handed
target-context supervision?*

LEAK AUDIT
  * embedding      : donor delta matrices only (k562_gwps / orion_hek293t), never the harness.
  * W and c        : fitted on `t["d"][train]`, the harness reference, TRAIN rows only.
  * score          : kappa / eps / mse_min on TEST rows only, and the incumbent is re-scored on
                     the SAME test rows so M is a like-for-like ratio.
  * no hyper-parameter is chosen using test rows inside a split; the (k, lambda) grid is printed
    whole and the cross-harness transfer is reported separately in `lin_report.py`.
"""
from __future__ import annotations

import json
import sys

import numpy as np

import lin_kit as K

KS = (5, 10, 20, 50)
LAMS = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)
NSPLIT = 30
EMB = ("k562_gwps", "orion_hek293t", "both")


def embeddings(t, kmax=max(KS)):
    """`{donor: (E[P, kmax], covered[P], evr, n_cov_donor)}` on this harness's targets."""
    out = {}
    for dn in ("k562_gwps", "orion_hek293t"):
        blk = K.donor_block(dn, t["genes"])
        US, V, evr, tot = K.embed_pca(blk["S"], kmax)
        idx = blk["index"]
        P = t["dh"].shape[0]
        E = np.zeros((P, kmax))
        cov = np.zeros(P, bool)
        for r, x in enumerate(t["targets"]):
            j = idx.get(str(x))
            if j is not None:
                E[r] = US[j]
                cov[r] = True
        out[dn] = dict(E=E, cov=cov, evr=evr, ndonor=blk["S"].shape[0],
                       ngene=blk["n_genes_hit"], V=V)
        print(f"    [{dn}] PCA on {blk['S'].shape[0]} donor perts x {blk['n_genes_hit']} "
              f"matched genes; evr top-5 {np.round(evr[:5], 4).tolist()}  "
              f"cum@5 {evr[:5].sum():.4f} @10 {evr[:10].sum():.4f} @20 {evr[:20].sum():.4f} "
              f"@50 {evr[:50].sum():.4f}; harness coverage {cov.sum()}/{P}")
    A, B = out["k562_gwps"], out["orion_hek293t"]
    out["both"] = dict(E=np.hstack([A["E"], B["E"]]), cov=A["cov"] | B["cov"],
                       evr=None, ndonor=None, ngene=None, V=None, concat=True)
    return out


def run(name):
    t = K.prep(name)
    P = t["dh"].shape[0]
    print(f"\n########## e1 SUPERVISED CEILING -- {name} (P={P}) ##########")
    emb = embeddings(t)
    Y = t["d"]                                   # target coords already zeroed by K.prep
    rows = []
    rng0 = np.random.default_rng(20260828)
    splits = []
    for s in range(NSPLIT):
        perm = np.random.default_rng(1000 + s).permutation(P)
        splits.append((perm[: P // 2], perm[P // 2:]))
    # incumbent, re-scored on each test half
    inc = np.array([K.shape(t["dh"], t, te)[0] for _, te in splits])
    print(f"    incumbent kappa on the same held-out halves: "
          f"{inc.mean():.5f} +/- {inc.std():.5f}  (full panel {K.shape(t['dh'], t)[0]:.5f})")

    for dn in EMB:
        Efull = emb[dn]["E"]
        kcap = Efull.shape[1] // (2 if dn == "both" else 1)
        for k in KS:
            if dn == "both":
                E = np.hstack([Efull[:, :k], Efull[:, kcap:kcap + k]])
            else:
                E = Efull[:, :k]
            for icept in (True, False):
                for lam_mult in LAMS:
                    kap, eps, r2i, r2o, wn, mm = [], [], [], [], [], []
                    for tr, te in splits:
                        Et, Yt = E[tr], Y[tr]
                        scale = float(np.einsum("ij,ij->", Et, Et)) / max(E.shape[1], 1)
                        lam = lam_mult * max(scale, 1e-30)
                        W, c = K.ridge_fit(Et, Yt, lam, icept)
                        DH = np.zeros_like(t["dh"])
                        DH[tr] = Et @ W + c
                        DH[te] = E[te] @ W + c
                        r2i.append(K.r2(Yt, DH[tr]))
                        r2o.append(K.r2(Y[te], DH[te]))
                        wn.append(float(np.linalg.norm(W)))
                        rr = K.shape_full(DH, t, te)
                        kap.append(rr["kappa"])
                        eps.append(rr["eps"])
                        mm.append(rr["mse_min"])
                    kap = np.array(kap)
                    Mr = float((kap / inc).mean())
                    ev = float(np.mean(eps))
                    rows.append(dict(harness=name, emb=dn, k=k, icept=icept, lam=lam_mult,
                                     kappa=float(kap.mean()), kappa_sd=float(kap.std()),
                                     M=Mr, M_sd=float((kap / inc).std()),
                                     eps=ev, eps_sd=float(np.std(eps)),
                                     r2_in=float(np.mean(r2i)), r2_out=float(np.mean(r2o)),
                                     wnorm=float(np.mean(wn)),
                                     mse_half=float(np.mean(mm)),
                                     ch_raw=float(K.challenge_raw(Mr, ev)),
                                     ch_score=float(K.challenge_score(Mr, ev)),
                                     inc_kappa=float(inc.mean())))
    return rows


def show(rows, name):
    print(f"\n  ---- {name}: full (k, lambda) grid, mean over {NSPLIT} splits ----")
    print(f"  {'emb':<14}{'k':>4}{'icpt':>6}{'lam':>9}{'kappa':>9}{'sd':>8}"
          f"{'M':>7}{'Msd':>7}{'eps':>9}{'R2 in':>9}{'R2 out':>9}{'|W|':>9}"
          f"{'mse1/2':>9}{'ch_raw':>9}{'ch_sc':>8}")
    for r in rows:
        print(f"  {r['emb']:<14}{r['k']:>4}{str(r['icept'])[0]:>6}{r['lam']:>9.0e}"
              f"{r['kappa']:>9.5f}{r['kappa_sd']:>8.5f}{r['M']:>7.3f}{r['M_sd']:>7.3f}"
              f"{r['eps']:>+9.5f}{r['r2_in']:>9.4f}{r['r2_out']:>9.4f}{r['wnorm']:>9.3g}"
              f"{r['mse_half']:>9.4f}{r['ch_raw']:>9.4f}{r['ch_score']:>8.4f}")


if __name__ == "__main__":
    allr = []
    for nm in ("h1-test", "orion"):
        rr = run(nm)
        show(rr, nm)
        allr += rr
    with open("work/lin_e1.json", "w") as f:
        json.dump(allr, f)
    print(f"\n  wrote work/lin_e1.json ({len(allr)} rows)", file=sys.stderr)
