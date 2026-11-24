"""189 Phase 3 — every falsifier the mandate names, on one set of lines of code.

F1  PERMUTATION NULL.  `theta_g` shuffled across genes, re-run on BOTH surfaces: the donor LOLO
    folds (where the prior is fitted and evaluated) and h1-test (where it is deployed).

F2  EXPRESSION-ONLY CONTROL.  The whole pipeline with `theta_g` replaced by the recipient's own
    control expression.  If expression does as well, the finding is "use expression".

F3  THE h1-test TRAP.  `data/donors/h1_2025` is BANNED and never opened.  Beyond that, the two
    Orion lines cover 84/85 of h1-test's targets and K562 covers 82/85, so the STRICT control is
    a prior refitted with those 85 targets deleted from every donor pair
    (`build189.py --drop-h1test-targets`).  If the h1-test numbers move, the prior was reading
    target identity through the gene-level statistic.

F4  PANEL DISJOINTNESS.  259 of the 300 challenge panel targets sit in the three-way donor pool,
    so the pilot's fitting pool is NOT disjoint from the deployment set the way 161's core5 was.
    The control is the prior refitted on non-panel targets only (already built as the `NP`
    columns of `pairs189.npz`), read on the same endpoints.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common189 import WORK, decile_spread, spearman
from harness189 import Surface189, check, paired_bootstrap
from members189 import precision, topn_call
from prior189 import LINES, MIN_N, fit_hier, load_pairs, raw_rate

LAM = 0.5           # the h1-selected best lambda from screen189.py


def variants():
    """(tag, theta, genes) for every prior variant the falsifiers need."""
    out = []
    for tag, kw in (("MAIN: all targets, all 10 pairs", dict()),
                    ("F3 strict: h1-test's 85 targets deleted", dict(name="pairs189.noh1")),
                    ("F4 panel-free: 300 panel targets deleted", dict(nonpanel=True))):
        genes, names, K, N, B = load_pairs(**kw)
        th, m, tau = fit_hier(K, N, list(range(len(names))))
        raw, n = raw_rate(K, N, list(range(len(names))))
        out.append((tag, genes, th, raw, n, tau))
    return out


def main() -> None:
    T = Surface189()
    check(T)
    base = np.log(np.maximum(np.abs(T.H), 1e-30))
    ref = T.reach_vec(np.abs(T.H))
    res = {}

    # ---- F3 / F4 on h1-test ----------------------------------------------------------
    print(f"\n{'=' * 112}\nF3 / F4 — prior variants read on h1-test  "
          f"(|H| incumbent 0.23533)")
    print(f"{'variant':<44}{'tau':>7}{'rho vs MAIN':>13}{'reach@lam=0.5':>15}"
          f"{'d vs |H|':>10}{'SE':>8}{'prec@50':>9}{'q_w':>9}")
    K0, q0 = T.mass_ratio(), T.qw()
    main_th = None
    for tag, genes, th, raw, n, tau in variants():
        thh = T.map_gene(genes, th)
        if main_th is None:
            main_th = thh
        key = np.exp(base + LAM * thh)
        r = T.reach(key)
        d, se, lo, hi, pb = paired_bootstrap(T.reach_vec(key), ref)
        pr, _ = precision(T, topn_call(T, key, 50))
        w = np.exp(LAM * thh)
        w = w * (K0 / T.mass_ratio(w))
        print(f"{tag:<44}{tau:>7.4f}{spearman(main_th, thh):>+13.4f}{r:>15.5f}"
              f"{d:>+10.5f}{se:>8.5f}{pr:>9.4f}{T.qw(w):>9.5f}")
        res[tag] = dict(tau=tau, reach=r, d=d, se=se, prec50=pr, qw=T.qw(w),
                        rho_vs_main=spearman(main_th, thh))

    # ---- F1 / F2 on the donor LOLO folds ---------------------------------------------
    print(f"\n{'=' * 112}\nF1 / F2 — donor LOLO folds: prior vs shuffled prior vs expression")
    genes, names, K, N, B = load_pairs()
    P = {p: j for j, p in enumerate(names)}
    rows = []
    for L in LINES:
        fc = [j for p, j in P.items() if L not in p.split("|")]
        th, m, tau = fit_hier(K, N, fc)
        nfit = N[:, fc].sum(1)
        for j in [j for p, j in P.items() if L in p.split("|")]:
            kk, nn = K[:, j], N[:, j]
            sel = (nn > MIN_N) & (nfit > MIN_N) & np.isfinite(B[L])
            y = kk[sel] / np.maximum(nn[sel], 1)
            real = decile_spread(th[sel], y)[2]
            nul = [decile_spread(np.random.default_rng(70 + s).permutation(th)[sel], y)[2]
                   for s in range(20)]
            expr = decile_spread(B[L][sel], y)[2]
            rows.append(dict(held=L, pair=names[j], real=real, null_mean=float(np.mean(nul)),
                             null_sd=float(np.std(nul, ddof=1)), expr=expr,
                             z=(real - np.mean(nul)) / max(np.std(nul, ddof=1), 1e-9)))
    print(f"{'held-out':<10}{'eval pair':<20}{'prior d10':>11}{'null d10':>10}{'null sd':>9}"
          f"{'z':>8}{'expression d10':>16}")
    for r in rows:
        print(f"{r['held']:<10}{r['pair']:<20}{r['real']:>+11.4f}{r['null_mean']:>+10.4f}"
              f"{r['null_sd']:>9.4f}{r['z']:>+8.1f}{r['expr']:>+16.4f}")
    mr = lambda k_: float(np.mean([r[k_] for r in rows]))
    print(f"{'MEAN':<30}{mr('real'):>+11.4f}{mr('null_mean'):>+10.4f}{mr('null_sd'):>9.4f}"
          f"{mr('z'):>+8.1f}{mr('expr'):>+16.4f}")

    # ---- F1 / F2 on h1-test, with a proper null distribution --------------------------
    print(f"\n{'=' * 112}\nF1 / F2 — h1-test: 50-seed permutation null for the ordering key")
    th = T.theta
    real = T.reach(np.exp(base + LAM * th))
    nul = [T.reach(np.exp(base + LAM * np.random.default_rng(2000 + s).permutation(th)))
           for s in range(50)]
    ze = (T.b_ctrl - T.b_ctrl.mean()) / T.b_ctrl.std()
    best_expr = max((T.reach(np.exp(base + mu * ze)), mu)
                    for mu in (-1, -0.5, -0.25, 0.25, 0.5, 1))
    print(f"  |H|                                     0.23533")
    print(f"  |H| x exp(+0.50 theta)                  {real:.5f}  "
          f"(+{real - 0.23533:.5f} over |H|)")
    print(f"  50-seed shuffled-theta null             mean {np.mean(nul):.5f} "
          f"sd {np.std(nul, ddof=1):.5f}  range [{min(nul):.5f}, {max(nul):.5f}]")
    print(f"  P(null >= real)                         "
          f"{float(np.mean(np.asarray(nul) >= real)):.3f}   "
          f"z = {(real - np.mean(nul)) / np.std(nul, ddof=1):+.2f}")
    print(f"  best expression tilt over the same grid {best_expr[0]:.5f} "
          f"at mu = {best_expr[1]:+.2f}  ({best_expr[0] - 0.23533:+.5f} over |H|)")

    res["lolo"] = rows
    res["h1_null"] = dict(real=real, null=nul, best_expr=list(best_expr))
    json.dump(res, open(os.path.join(WORK, "falsify189.json"), "w"), indent=2)
    print("\nwrote work/falsify189.json")


if __name__ == "__main__":
    main()
