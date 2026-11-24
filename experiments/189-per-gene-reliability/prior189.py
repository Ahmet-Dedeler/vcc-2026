"""189 Phase 1c — the prior itself: hierarchical fit, LOLO scaling, and what shrinkage does.

THE MODEL
---------
`r_g` is a binomial proportion with a known `n_g`, so the honest estimate is a posterior mean.
Two estimators are built and both are carried through every downstream table:

  RAW      r_g = (sum_p k_gp) / (sum_p n_gp)                        -- the pilot's statistic
  EB       beta-binomial posterior mean, (a, b) by marginal MLE     -- pooled counts
  HIER     logit(q_gp) = m_p + theta_g,  theta_g ~ N(0, tau^2)      -- the delivered prior

HIER exists because the ten cross-LINE pairs do NOT share a base rate: pooled agreement runs
0.5221 (HEK293T|K562) to 0.6147 (K562|JURKAT), a 0.093 spread, and genes are carried by
different SUBSETS of pairs (10,602 genes in the Orion pair against 6,414 in K562|RPE1).  Pooling
raw counts therefore ranks a gene partly by which pairs happen to carry it -- a coverage
artefact with nothing to do with reliability.  A per-pair intercept removes it exactly, and the
Gaussian prior on `theta_g` is the same shrinkage the beta-binomial applies, on the logit scale.

`tau` is fitted by the standard Laplace EB update, so the amount of shrinkage is estimated from
the data rather than chosen.

LEAVE-ONE-LINE-OUT
------------------
For a held-out LINE L the fit uses only pairs with neither endpoint equal to L, and is scored
against the observed per-gene agreement of the pairs that DO involve L.  Five folds
(HCT116 / HEK293T / K562 / RPE1 / JURKAT).  Deployment on h1-test is leave-H1-out by
construction: H1 is not a fitting line and `data/donors/h1_2025` is never opened.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

from common189 import WORK, decile_spread, eb_shrink, spearman

LINES = ["HCT116", "HEK293T", "K562", "RPE1", "JURKAT"]
MIN_N = 50          # the pilot's support threshold


# --------------------------------------------------------------------------------------
def load_pairs(nonpanel: bool = False, name: str = "pairs189"):
    z = np.load(os.path.join(WORK, f"{name}.npz"), allow_pickle=True)
    pre = "kNP::" if nonpanel else "k::"
    names = sorted(k[len(pre):] for k in z.files if k.startswith(pre))
    K = np.column_stack([z[pre + p] for p in names])
    N = np.column_stack([z[("nNP::" if nonpanel else "n::") + p] for p in names])
    genes = np.array([str(g) for g in z["genes"]])
    b = {L: z[f"b::{L}"] for L in LINES}
    return genes, names, K, N, b


def fit_hier(K: np.ndarray, N: np.ndarray, cols, iters: int = 60, tau0: float = 0.3):
    """logit(q_gp) = m_p + theta_g, theta_g ~ N(0, tau^2).  Returns (theta, m, tau).

    Alternating Newton on (theta | m) and (m | theta) with a Laplace EB update for tau.
    """
    k = K[:, cols]
    n = N[:, cols]
    G, P = k.shape
    live = n.sum(1) > 0
    m = np.zeros(P)
    for j in range(P):
        p = k[:, j].sum() / max(n[:, j].sum(), 1)
        m[j] = np.log(p / (1 - p))
    th = np.zeros(G)
    tau = tau0
    for _ in range(iters):
        for _ in range(4):                                    # theta | m
            e = np.exp(np.clip(m[None, :] + th[:, None], -30, 30))
            q = e / (1 + e)
            grad = (k - n * q).sum(1) - th / tau ** 2
            hess = -(n * q * (1 - q)).sum(1) - 1.0 / tau ** 2
            th = th - grad / hess
            th = np.clip(th, -6, 6)
        e = np.exp(np.clip(m[None, :] + th[:, None], -30, 30))
        q = e / (1 + e)
        info = (n * q * (1 - q)).sum(1)
        var = 1.0 / (info + 1.0 / tau ** 2)
        tau_new = float(np.sqrt(np.mean(th[live] ** 2 + var[live])))
        for _ in range(4):                                    # m | theta
            e = np.exp(np.clip(m[None, :] + th[:, None], -30, 30))
            q = e / (1 + e)
            g_ = (k - n * q).sum(0)
            h_ = -(n * q * (1 - q)).sum(0)
            m = m - g_ / np.where(h_ == 0, -1.0, h_)
        if abs(tau_new - tau) < 1e-9:
            tau = tau_new
            break
        tau = tau_new
    th[~live] = 0.0
    return th, m, float(tau)


def raw_rate(K, N, cols):
    k, n = K[:, cols].sum(1), N[:, cols].sum(1)
    return np.where(n > 0, k / np.maximum(n, 1), np.nan), n


# --------------------------------------------------------------------------------------
def main() -> None:
    genes, names, K, N, B = load_pairs()
    P = {p: j for j, p in enumerate(names)}
    print(f"{len(genes)} working genes, {len(names)} cross-LINE pairs\n")

    def cols_excluding(L):
        return [j for p, j in P.items() if L not in p.split("|")]

    def cols_involving(L):
        return [j for p, j in P.items() if L in p.split("|")]

    # ---- 1. the delivered prior: every pair, H1 held out by construction ---------------
    th_all, m_all, tau_all = fit_hier(K, N, list(range(len(names))))
    raw_all, n_all = raw_rate(K, N, list(range(len(names))))
    ok = n_all > MIN_N
    eb_all, a, b, mom = eb_shrink(K.sum(1), N.sum(1))
    print(f"HIER on all 10 pairs: tau {tau_all:.4f}   pair offsets "
          + " ".join(f"{p.split('|')[0][:3]}/{p.split('|')[1][:3]}:{1/(1+np.exp(-mm)):.3f}"
                     for p, mm in zip(names, m_all)))
    print(f"EB beta-binomial on pooled counts: a {a:.1f} b {b:.1f} "
          f"(prior mean {a/(a+b):.4f}, prior n {a+b:.0f}; MoM start {mom[0]:.1f}/{mom[1]:.1f})")
    print(f"raw r_g over {int(ok.sum())} supported genes: mean {np.nanmean(raw_all[ok]):.4f} "
          f"sd {np.nanstd(raw_all[ok], ddof=1):.4f}")
    print(f"EB  r_g: mean {eb_all[ok].mean():.4f} sd {eb_all[ok].std(ddof=1):.4f}   "
          f"shrink factor {eb_all[ok].std(ddof=1)/np.nanstd(raw_all[ok], ddof=1):.3f}")
    print(f"spearman(raw, EB) {spearman(raw_all[ok], eb_all[ok]):+.4f}   "
          f"spearman(raw, HIER) {spearman(raw_all[ok], th_all[ok]):+.4f}   "
          f"spearman(EB, HIER) {spearman(eb_all[ok], th_all[ok]):+.4f}")

    # ---- 2. LOLO folds ---------------------------------------------------------------
    print(f"\n{'=' * 96}\nLEAVE-ONE-LINE-OUT  (fit on pairs without L, score on pairs with L)")
    print(f"{'held-out L':<10}{'fit pairs':>10}{'tau':>8}   "
          f"{'eval pair':<18}{'rho':>8}{'d10 spread':>12}{'raw d10':>10}{'EB d10':>9}"
          f"{'expr d10':>10}")
    fold_rows = []
    for L in LINES:
        fc = cols_excluding(L)
        th, m, tau = fit_hier(K, N, fc)
        rawL, nfit = raw_rate(K, N, fc)
        ebL, *_ = eb_shrink(K[:, fc].sum(1), N[:, fc].sum(1))
        for j in cols_involving(L):
            other = [x for x in names[j].split("|") if x != L][0]
            kk, nn = K[:, j], N[:, j]
            m_ok = (nn > MIN_N) & (nfit > MIN_N)
            y = np.where(nn > 0, kk / np.maximum(nn, 1), np.nan)
            expr = B[L]                                  # the RECIPIENT's own control expression
            sel = m_ok & np.isfinite(expr)
            rho = spearman(th[sel], y[sel])
            _, _, d_h = decile_spread(th[sel], y[sel])
            _, _, d_r = decile_spread(rawL[sel], y[sel])
            _, _, d_e = decile_spread(ebL[sel], y[sel])
            _, _, d_x = decile_spread(expr[sel], y[sel])
            fold_rows.append(dict(held=L, other=other, n_genes=int(sel.sum()), rho=rho,
                                  d_hier=d_h, d_raw=d_r, d_eb=d_e, d_expr=d_x, tau=tau,
                                  fit_pairs=len(fc)))
            print(f"{L:<10}{len(fc):>10}{tau:>8.4f}   {L + '|' + other:<18}{rho:>+8.4f}"
                  f"{d_h:>+12.4f}{d_r:>+10.4f}{d_e:>+9.4f}{d_x:>+10.4f}")
    mr = lambda k_: float(np.mean([r[k_] for r in fold_rows]))
    print(f"{'MEAN':<10}{'':>10}{'':>8}   {'':<18}{mr('rho'):>+8.4f}{mr('d_hier'):>+12.4f}"
          f"{mr('d_raw'):>+10.4f}{mr('d_eb'):>+9.4f}{mr('d_expr'):>+10.4f}")

    # ---- 3. how the prior scales with the number of fitting LINES ---------------------
    print(f"\n{'=' * 96}\nSCALING: held-out predictive power vs number of fitting LINES")
    print(f"{'fit lines':>10}{'fit pairs':>11}{'subsets':>9}{'mean rho':>10}"
          f"{'mean d10 spread':>17}{'mean tau':>10}")
    scale_rows = []
    for nfit_lines in (2, 3, 4):
        rhos, spreads, taus, nsub = [], [], [], 0
        for L in LINES:
            rest = [x for x in LINES if x != L]
            for sub in itertools.combinations(rest, nfit_lines):
                fc = [P[p] for p in names
                      if p.split("|")[0] in sub and p.split("|")[1] in sub]
                if not fc:
                    continue
                nsub += 1
                th, m, tau = fit_hier(K, N, fc)
                nfitv = N[:, fc].sum(1)
                taus.append(tau)
                for j in cols_involving(L):
                    kk, nn = K[:, j], N[:, j]
                    sel = (nn > MIN_N) & (nfitv > MIN_N)
                    y = np.where(nn > 0, kk / np.maximum(nn, 1), np.nan)
                    rhos.append(spearman(th[sel], y[sel]))
                    spreads.append(decile_spread(th[sel], y[sel])[2])
        scale_rows.append(dict(lines=nfit_lines, pairs=nfit_lines * (nfit_lines - 1) // 2,
                               subsets=nsub, rho=float(np.mean(rhos)),
                               spread=float(np.mean(spreads)), tau=float(np.mean(taus))))
        print(f"{nfit_lines:>10}{nfit_lines*(nfit_lines-1)//2:>11}{nsub:>9}"
              f"{np.mean(rhos):>+10.4f}{np.mean(spreads):>+17.4f}{np.mean(taus):>10.4f}")

    # ---- 4. what shrinkage does to the decile spread ---------------------------------
    print(f"\n{'=' * 96}\nWHAT SHRINKAGE DOES  (fold-mean held-out decile spread)")
    print(f"  raw fraction              {mr('d_raw'):+.4f}")
    print(f"  beta-binomial EB          {mr('d_eb'):+.4f}")
    print(f"  hierarchical (delivered)  {mr('d_hier'):+.4f}")
    print(f"  recipient expression      {mr('d_expr'):+.4f}")

    np.savez_compressed(os.path.join(WORK, "prior189.npz"), genes=genes, pairs=np.array(names),
                        theta=th_all, m=m_all, tau=tau_all, raw=raw_all, eb=eb_all,
                        n_all=n_all, K=K, N=N,
                        **{f"b::{L}": B[L] for L in LINES})
    json.dump(dict(tau_all=tau_all, a=a, b=b, folds=fold_rows, scaling=scale_rows,
                   pair_names=names, pair_offsets=list(map(float, m_all))),
              open(os.path.join(WORK, "prior189.json"), "w"), indent=2)
    print("\nwrote work/prior189.npz + work/prior189.json")


if __name__ == "__main__":
    main()
