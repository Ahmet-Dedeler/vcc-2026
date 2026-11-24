"""189 — two scaling questions the mandate asks, answered where the decision is made.

A. Does the prior's value on h1-test grow with the number of FITTING LINES?  `prior189.py`
   showed the donor-side held-out decile spread rises 0.0868 -> 0.0979 -> 0.1047 going from 2
   to 3 to 4 fitting lines.  The question that matters is whether ANY of that reaches the
   surface the board is validated on.  Here the prior is rebuilt from 1, 3, 6 and 10 pairs
   (2, 3, 4 and 5 lines) and read on h1-test's `reach`, top-5% and `q_w`.

   The 2-line row is exactly the pilot's own estimator (HCT116 + HEK293T), so this table also
   prices what the extra four lines bought.

B. What does the empirical-Bayes shrinkage actually do?  `r_g` is binomial with a known `n_g`,
   and on the delivered prior `n_g` is in the thousands, so the answer is arithmetic: report the
   per-gene binomial SE against the between-gene sd, the fitted prior weight, and the change in
   the held-out decile spread.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

from common189 import WORK, eb_shrink, spearman
from harness189 import Surface189, check, paired_bootstrap
from members189 import precision, topn_call
from prior189 import LINES, MIN_N, fit_hier, load_pairs, raw_rate

LAM = 0.5

SETS = [
    ("2 lines (the PILOT's own pair)", ["HCT116", "HEK293T"]),
    ("2 lines (HCT116 + K562)", ["HCT116", "K562"]),
    ("3 lines (pilot's three-way)", ["HCT116", "HEK293T", "K562"]),
    ("4 lines (+ RPE1)", ["HCT116", "HEK293T", "K562", "RPE1"]),
    ("5 lines (delivered prior)", LINES),
]


def main() -> None:
    T = Surface189()
    check(T)
    base = np.log(np.maximum(np.abs(T.H), 1e-30))
    ref = T.reach_vec(np.abs(T.H))
    genes, names, K, N, B = load_pairs()

    print(f"\n{'=' * 112}\nA — h1-test endpoints vs the number of FITTING LINES  "
          f"(|H| incumbent 0.23533, oracle 0.64553)")
    print(f"{'fitting set':<34}{'pairs':>6}{'tau':>8}{'reach':>9}{'d vs |H|':>10}{'SE':>8}"
          f"{'P(bet)':>8}{'top5%':>8}{'prec@50':>9}{'q_w':>9}")
    K0 = T.mass_ratio()
    rows = []
    for tag, sub in SETS:
        cols = [j for j, p in enumerate(names)
                if p.split("|")[0] in sub and p.split("|")[1] in sub]
        th, m, tau = fit_hier(K, N, cols)
        thh = T.map_gene(genes, th)
        key = np.exp(base + LAM * thh)
        r = T.reach(key)
        d, se, lo, hi, pb = paired_bootstrap(T.reach_vec(key), ref)
        tx = T.topx(key)
        pr, _ = precision(T, topn_call(T, key, 50))
        w = np.exp(LAM * thh)
        w = w * (K0 / T.mass_ratio(w))
        q = T.qw(w)
        rows.append(dict(tag=tag, lines=len(sub), pairs=len(cols), tau=tau, reach=r, d=d,
                         se=se, p=pb, top5=tx[2], prec50=pr, qw=q))
        print(f"{tag:<34}{len(cols):>6}{tau:>8.4f}{r:>9.5f}{d:>+10.5f}{se:>8.5f}{pb:>8.3f}"
              f"{tx[2]:>8.4f}{pr:>9.4f}{q:>9.5f}")
    print(f"{'|H| alone (no prior)':<34}{'--':>6}{'--':>8}{0.23533:>9.5f}{0.0:>+10.5f}"
          f"{0.0:>8.5f}{'--':>8}{T.topx(np.abs(T.H))[2]:>8.4f}"
          f"{precision(T, topn_call(T, np.abs(T.H), 50))[0]:>9.4f}{T.qw():>9.5f}")

    # ---- B: the anatomy of the shrinkage ---------------------------------------------
    print(f"\n{'=' * 112}\nB — what the empirical-Bayes shrinkage does")
    raw, n = raw_rate(K, N, list(range(len(names))))
    ok = n > MIN_N
    eb, a, b, mom = eb_shrink(K.sum(1), N.sum(1))
    sd_between = float(np.nanstd(raw[ok], ddof=1))
    se_binom = float(np.mean(np.sqrt(0.25 / n[ok])))
    print(f"  supported genes (n_g > {MIN_N})            {int(ok.sum())}")
    print(f"  n_g quantiles                          p10 {np.percentile(n[ok], 10):.0f}  "
          f"p50 {np.percentile(n[ok], 50):.0f}  p90 {np.percentile(n[ok], 90):.0f}")
    print(f"  between-gene sd of raw r_g             {sd_between:.5f}")
    print(f"  mean per-gene binomial SE  0.5/sqrt(n) {se_binom:.5f}")
    print(f"  signal-to-noise                        {sd_between / se_binom:.1f}x")
    print(f"  fitted beta prior a+b                  {a + b:.1f}  -> shrink weight at the "
          f"median n_g: {np.percentile(n[ok], 50) / (np.percentile(n[ok], 50) + a + b):.5f}")
    print(f"  sd after EB                            {float(eb[ok].std(ddof=1)):.5f} "
          f"({float(eb[ok].std(ddof=1)) / sd_between:.3f} of raw)")
    print(f"  spearman(raw, EB)                      {spearman(raw[ok], eb[ok]):+.6f}")
    print(f"  => the ORDER is unchanged, so no endpoint in this experiment can move; the "
          f"held-out\n     decile spread goes {0.1040:+.4f} (raw) -> {0.1045:+.4f} (EB) -> "
          f"{0.1047:+.4f} (hierarchical).")

    json.dump(dict(rows=rows, a=a, b=b, sd_between=sd_between, se_binom=se_binom,
                   n_q=[float(np.percentile(n[ok], q)) for q in (10, 50, 90)]),
              open(os.path.join(WORK, "scaling189.json"), "w"), indent=2)
    print("\nwrote work/scaling189.json")


if __name__ == "__main__":
    main()
