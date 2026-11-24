"""293 -- from the measured ladder to the two ceilings a 400-cell challenge truth implies.

Variance model (the standard attenuation model, stated so it can be checked):  a construct k's true
response is  c_k = g + u_k  (g the gene-level mean over the construct distribution, u_k the construct-
specific part, orthogonal in expectation), and the scored truth is a noisy 400-cell measurement
t_k = c_k + e.  With s2g / s2u / s2e the three variances and S = s2g + s2u + s2e:

  (a)  same construct, two independent measurements   = (s2g + s2u) / S
  (b)  two constructs of one gene                     = s2g / S
  rho_guide = (b)/(a)                                 = s2g / (s2g + s2u)

  perfect CONSTRUCT-level prediction vs the scored truth   cos = sqrt((a))
  perfect GENE-level prediction vs the scored truth        cos = sqrt((b))     [the zero-shot ceiling:
      no 2026 construct is on disk, so no model can know which construct produced the truth]

(a) and (b) are read at n = 100 and 200 per side and Spearman-Brown-extrapolated to the challenge's
400-cell budget (k = 400/n applied to the measured half-split value; SB of a 200/200 split with k = 2 is
the reliability of a 400-cell measurement against another 400-cell measurement).

    PYTHONDONTWRITEBYTECODE=1 uv run python ceiling293.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib293 as X  # noqa: E402
from lib293 import L  # noqa: E402

P = L.Logger(X.WORK / "ceiling293.log")
TAGS = {"native (2025 files, 54,256 median UMI/cell)": "",
        "CHALLENGE REGIME (thinned to 19,965 UMI/cell; h1-test 19,965, vcc2026 controls 20,340)": "_thin0.368"}


def boot(fn, arrs, n=L.NBOOT, seed=L.SEED):
    """Paired target bootstrap of an arbitrary statistic of several aligned per-target arrays."""
    m = np.ones(len(arrs[0]), dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    arrs = [np.asarray(a)[m] for a in arrs]
    rng = np.random.default_rng(seed)
    ii = rng.integers(0, arrs[0].size, size=(n, arrs[0].size))
    d = np.asarray([fn(*[a[ii[j]] for a in arrs]) for j in range(n)])
    return dict(point=float(fn(*arrs)), lo=float(np.nanpercentile(d, 2.5)),
                hi=float(np.nanpercentile(d, 97.5)), n=int(arrs[0].size))


def main() -> None:
    t0 = time.time()
    P(f"293 ceiling293.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    out: dict = {}
    for label, tag in TAGS.items():
        A = json.loads((X.WORK / f"agg293{tag}.json").read_text())
        P(f"\n=== {label}")
        out[tag or "native"] = {"label": label, "rungs": {}}
        for n in (100, 200):
            r = A["rungs"][str(n)]
            rows = list(csv.DictReader(open(X.WORK / f"construct293{tag}.csv")))
            rows = [x for x in rows if int(x["n"]) == n]
            a = np.asarray([float(x["a_cos"]) for x in rows])
            b = np.asarray([float(x["b_cos"]) for x in rows])
            tg = [x["target"] for x in rows]
            k = 400.0 / n

            def sbm(v):
                return float(np.median([L.spearman_brown(z, k) for z in v]))

            res = dict(
                n=n, n_targets=len(tg), targets=tg,
                a_measured=L.boot_median(a), b_measured=L.boot_median(b),
                a_400=boot(lambda x: sbm(x), [a]),
                b_400=boot(lambda y: sbm(y), [b]),
                rho_guide=boot(lambda x, y: float(np.median(y) / np.median(x)), [a, b]),
                # the two ceilings, propagated through the same bootstrap resample
                ceil_construct=boot(lambda x: float(np.sqrt(max(sbm(x), 0.0))), [a]),
                ceil_gene=boot(lambda y: float(np.sqrt(max(sbm(y), 0.0))), [b]),
                ceil_gene_via_rho=boot(
                    lambda x, y: float(np.sqrt(max(sbm(x) * (np.median(y) / np.median(x)), 0.0))), [a, b]),
            )
            out[tag or "native"]["rungs"][n] = res
            P(f"  n={n}/side, {len(tg)} targets")
            P(f"    (a) {L.fmt(res['a_measured'])}  ->  400/400 SB {L.fmt(res['a_400'])}")
            P(f"    (b) {L.fmt(res['b_measured'])}  ->  400/400 SB {L.fmt(res['b_400'])}")
            P(f"    rho_guide {L.fmt(res['rho_guide'])}")
            P(f"    CEILING vs a 400-cell truth:  perfect CONSTRUCT-level prediction cos "
              f"{L.fmt(res['ceil_construct'])}   perfect GENE-level (zero-shot) cos {L.fmt(res['ceil_gene'])} "
              f"[via rho: {L.fmt(res['ceil_gene_via_rho'])}]")

    # per-target table at the challenge regime, n = 200 and 100
    P("\n=== per-target, CHALLENGE REGIME")
    for n in (100, 200):
        rows = [x for x in csv.DictReader(open(X.WORK / "construct293_thin0.368.csv")) if int(x["n"]) == n]
        rows.sort(key=lambda x: -float(x["b_cos"]))
        P(f"\n  n={n}/side, sorted by (b):")
        P(f"    {'target':9s} {'split':11s} {'cells 1/2':13s} {'(a)':>8s} {'(b)':>8s} {'(d)':>8s} {'b/a':>8s} "
          f"{'pds a':>7s} {'pds b':>7s} {'kd big':>7s} {'kd small':>8s}")
        for x in rows:
            P(f"    {x['target']:9s} {x['split']:11s} {x['cells_big']:>6s}/{x['cells_small']:<6s} "
              f"{float(x['a_cos']):+8.4f} {float(x['b_cos']):+8.4f} "
              f"{(float(x['d_cos']) if x['d_cos'] not in ('', 'nan') else float('nan')):+8.4f} "
              f"{float(x['rho_ratio']):+8.3f} {float(x['a_pds']):7.3f} {float(x['b_pds']):7.3f} "
              f"{float(x['kd_big']):+7.3f} {float(x['kd_small']):+8.3f}")

    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, X.WORK / "ceiling293.json")
    P(f"\ndone {out['seconds']}s")
    P.save()


if __name__ == "__main__":
    main()
