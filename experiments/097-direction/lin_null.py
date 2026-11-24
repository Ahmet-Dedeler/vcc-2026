"""Falsification battery for e2's positive.  Four controls, each able to kill the result.

e2's winner is `k562_gwps -> orion_hek293t, k=50, intercept, lambda=1e0`, M = 1.90 (h1) /
2.03 (orion).  A positive on 24 of 543 arms is exactly the shape a multiple-comparison artefact
takes, and RESEARCH-METHOD rule 3 says a mechanism that was never printed is not a result.  So:

(N1) PERMUTATION NULL.  Refit W with donor B's rows SHUFFLED against donor A's, breaking the
     perturbation correspondence and nothing else.  If the shuffled map still scores, the gain
     is not in the learned map.

(N2) THE "JUST DROP THE WHITENER" CONTROL.  The incumbent's whitening exponent gamma = 0.65 was
     tuned on `pds`, a RETRIEVAL metric, long before `mse` left its dead zone.  A plain
     per-donor-renormalised average of the SAME two donors, unwhitened, is the cheapest possible
     competitor.  If that already reaches M ~ 1.9, the Ahlmann-Eltze form adds nothing.

(N3) EMISSION-NOISE CORRECTION.  The incumbent's kappa is read off an EMITTED counts file, so
     its budget `B` carries the emitter's own sampling variance, which every analytic arm here
     is free of.  `sum_p w_p Chat_p` is precisely that variance (it is the prediction's own
     jackknife trace, which the scorer already computes), so the incumbent's kappa can be
     corrected upward by `sqrt(B / (B - sum_p w_p Chat_p (1 - 1/n)))` and the comparison redone
     on equal terms.

(N4) PER-PERTURBATION GEOMETRY.  Distribution, not summary statistic: the per-row cosine
     between prediction and reference, for the incumbent and for the winner.
"""
from __future__ import annotations

import numpy as np

import lin_kit as K
from lin_e2 import grams, solve

WIN = dict(A="k562_gwps", B="orion_hek293t", k=50, icept=True, lam=1e0)
ALSO = [dict(A="k562_gwps", B="orion_hek293t", k=20, icept=True, lam=1e0),
        dict(A="k562_gwps", B="orion_hek293t", k=10, icept=True, lam=1e-1),
        dict(A="k562_gwps", B="orion_hek293t", k=50, icept=False, lam=1e0)]
NPERM = 20


def donor_S_on_targets(t, dn):
    blk = K.donor_block(dn, t["genes"])
    P = t["dh"].shape[0]
    M = np.zeros((P, t["genes"].size))
    cov = np.zeros(P, bool)
    for r, x in enumerate(t["targets"]):
        j = blk["index"].get(str(x))
        if j is not None:
            M[r] = blk["S"][j]
            cov[r] = True
    med = float(np.median(np.linalg.norm(blk["S"][:, blk["hit"]], axis=1))) or 1.0
    return M, cov, med, blk


def fit_arm(t, spec, perm=None):
    """Fit e2's W between two donors (optionally with B's rows permuted) and score it."""
    ba = K.donor_block(spec["A"], t["genes"])
    bb = K.donor_block(spec["B"], t["genes"])
    shared = [p for p in ba["index"] if p in bb["index"]]
    ra = np.array([ba["index"][p] for p in shared])
    rb = np.array([bb["index"][p] for p in shared])
    if perm is not None:
        rb = rb[perm]
    US, V, evr, tot = K.embed_pca(ba["S"], 50)
    g = grams(US[ra], bb["S"][rb], spec["icept"])
    W, c, r2f, wn = solve(*g, spec["k"], spec["lam"])
    P = t["dh"].shape[0]
    E = np.zeros((P, 50))
    cov = np.zeros(P, bool)
    for r, x in enumerate(t["targets"]):
        j = ba["index"].get(str(x))
        if j is not None:
            E[r] = US[j]
            cov[r] = True
    DH = E[:, :spec["k"]] @ W + c
    DH[~cov] = 0.0
    return DH, K.shape_full(DH, t), r2f, wn, len(shared)


def main():
    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P = t["dh"].shape[0]
        r0 = K.shape_full(t["dh"], t)
        k0 = r0["kappa"]
        print(f"\n################ {name}  (P={P}) ################")
        print(f"  incumbent kappa {k0:.5f}  eps {r0['eps']:+.5f}  B {r0['B']:.5f}")

        # -------- N3 --------------------------------------------------------------------
        noise = float((t["Chat"] * t["w"]).sum()) * (1 - 1 / P)
        Btrue = r0["B"] - noise
        corr = np.sqrt(r0["B"] / Btrue) if Btrue > 0 else np.nan
        print(f"\n  (N3) EMISSION-NOISE CORRECTION")
        print(f"    sum_p w_p Chat_p (1-1/n) = {noise:.5f}   "
              f"= {noise/r0['B']*100:.1f}% of the incumbent's budget B")
        print(f"    signal-only budget B_true = {Btrue:.5f}  ->  kappa correction x{corr:.4f}")
        print(f"    noise-corrected incumbent kappa = {k0*corr:.5f}   "
              f"(every analytic arm below must clear THIS, not {k0:.5f})")

        # -------- N2 --------------------------------------------------------------------
        print(f"\n  (N2) PLAIN UNWHITENED DONOR AVERAGES (no fit, no whitener, no PCA)")
        print(f"    {'arm':<52}{'kappa':>9}{'M':>7}{'M/corr':>8}{'eps':>10}{'cov':>6}")
        pieces = {}
        for dn in ("k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential"):
            if dn == t["banned"]:
                continue
            M, cov, med, _ = donor_S_on_targets(t, dn)
            pieces[dn] = (M / med, cov)
        combos = [("k562_gwps",), ("orion_hek293t",), ("k562_gwps", "orion_hek293t"),
                  ("k562_gwps", "orion_hek293t", "nadig_jurkat"),
                  ("k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential")]
        for cb in combos:
            cb = tuple(x for x in cb if x in pieces)
            if not cb:
                continue
            acc = np.zeros((P, t["genes"].size))
            cnt = np.zeros(P)
            for dn in cb:
                Mn, cov = pieces[dn]
                acc[cov] += Mn[cov]
                cnt[cov] += 1
            acc[cnt > 0] /= cnt[cnt > 0, None]
            r = K.shape_full(acc, t)
            print(f"    {'renorm-mean of ' + '+'.join(cb):<52}{r['kappa']:>9.5f}"
                  f"{r['kappa']/k0:>7.3f}{r['kappa']/(k0*corr):>8.3f}{r['eps']:>+10.5f}"
                  f"{int((cnt>0).sum()):>6}")

        # -------- the winner and its neighbours, then N1 ---------------------------------
        print(f"\n  (N1) THE WINNING ARM, THEN THE SAME ARM WITH THE PAIRING SHUFFLED")
        print(f"    {'arm':<52}{'kappa':>9}{'M':>7}{'M/corr':>8}{'eps':>10}{'r2fit':>9}{'|W|':>9}")
        for spec in [WIN] + ALSO:
            DH, r, r2f, wn, ns = fit_arm(t, spec)
            lbl = (f"{spec['A']}->{spec['B']} k={spec['k']} "
                   f"ic={str(spec['icept'])[0]} lam={spec['lam']:.0e}")
            print(f"    {lbl:<52}{r['kappa']:>9.5f}{r['kappa']/k0:>7.3f}"
                  f"{r['kappa']/(k0*corr):>8.3f}{r['eps']:>+10.5f}{r2f:>9.4f}{wn:>9.3g}")
            if spec is WIN:
                ns_shared = ns
                win_DH = DH
        rng = np.random.default_rng(7)
        ks, es, ms = [], [], []
        for s in range(NPERM):
            perm = rng.permutation(ns_shared)
            _, r, _, _, _ = fit_arm(t, WIN, perm)
            ks.append(r["kappa"])
            es.append(r["eps"])
            ms.append(r["mse_min"])
        ks, es = np.array(ks), np.array(es)
        print(f"    {'^^ SHUFFLED pairing, ' + str(NPERM) + ' seeds':<52}"
              f"{np.nanmean(ks):>9.5f}{np.nanmean(ks)/k0:>7.3f}"
              f"{np.nanmean(ks)/(k0*corr):>8.3f}{np.nanmean(es):>+10.5f}")
        print(f"      shuffled kappa: mean {np.nanmean(ks):+.5f}  sd {np.nanstd(ks):.5f}  "
              f"min {np.nanmin(ks):+.5f}  max {np.nanmax(ks):+.5f}   "
              f"({int((ks >= K.shape(win_DH, t)[0]).sum())}/{NPERM} >= the real fit)")
        print(f"      shuffled eps  : mean {np.nanmean(es):+.5f}  max {np.nanmax(es):+.5f}"
              f"      shuffled mse_min mean {np.nanmean(ms):.4f} "
              f"(real {K.shape_full(win_DH, t)['mse_min']:.4f}, "
              f"incumbent {r0['mse_min']:.4f})")

        # -------- N4 --------------------------------------------------------------------
        def rowcos(H):
            a = np.linalg.norm(H, axis=1)
            b = np.linalg.norm(t["d"], axis=1)
            return np.einsum("ij,ij->i", H, t["d"]) / np.maximum(a * b, 1e-30)
        ci, cw = rowcos(t["dh"]), rowcos(win_DH)
        print(f"\n  (N4) PER-PERTURBATION COSINE, distribution not summary")
        for lbl, c in (("incumbent", ci), ("e2 winner", cw)):
            q = np.percentile(c, [5, 25, 50, 75, 95])
            print(f"    {lbl:<12} mean {c.mean():+.4f}  q05..q95 "
                  f"[{q[0]:+.4f} {q[1]:+.4f} {q[2]:+.4f} {q[3]:+.4f} {q[4]:+.4f}]  "
                  f"frac>0 {np.mean(c > 0):.3f}")
        print(f"    rank correlation between the two per-row cosines: "
              f"{np.corrcoef(np.argsort(np.argsort(ci)), np.argsort(np.argsort(cw)))[0,1]:+.3f}")
        K._CACHE.clear()


if __name__ == "__main__":
    main()
