"""293 -- POST-HOC (PREREGISTRATION sec 5, labelled as such): does the weaker construct's on-target
knockdown explain the between-target spread in construct agreement?

The knockdown is the construct's OWN gene column, which is removed from every scored statistic, so the
stratification is not circular.  Threshold: |lfc| >= 1.0 on the log1p(50,000 p) scale for BOTH constructs
("both constructs work"), plus a median split as a threshold-free check.

    PYTHONDONTWRITEBYTECODE=1 uv run python strat293.py
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

P = L.Logger(X.WORK / "strat293.log")
TAGS = {"native": "", "challenge": "_thin0.368"}


def block(rows, name, n):
    a = np.asarray([float(x["a_cos"]) for x in rows])
    b = np.asarray([float(x["b_cos"]) for x in rows])
    a2 = np.asarray([float(x["a2_cos"]) if x.get("a2_cos") not in (None, "", "nan") else np.nan
                     for x in rows])
    tg = [x["target"] for x in rows]
    if a.size == 0:
        return None
    k = 400.0 / n
    sb_a = float(np.median([L.spearman_brown(v, k) for v in a]))
    sb_b = float(np.median([L.spearman_brown(v, k) for v in b]))
    rng = np.random.default_rng(L.SEED)
    ii = rng.integers(0, a.size, size=(L.NBOOT, a.size))
    rr = np.median(b[ii], 1) / np.median(a[ii], 1)
    d = dict(name=name, n_targets=len(tg), targets=tg,
             a=L.boot_median(a), b=L.boot_median(b),
             rho=dict(point=float(np.median(b) / np.median(a)), lo=float(np.percentile(rr, 2.5)),
                      hi=float(np.percentile(rr, 97.5))),
             a_400=sb_a, b_400=sb_b,
             a2=L.boot_median(a2),
             rho_dis=(float(np.median(b) / np.sqrt(max(np.median(a), 0.0) * max(np.nanmedian(a2), 0.0)))
                      if np.isfinite(np.nanmedian(a2)) and np.nanmedian(a2) > 0 and np.median(a) > 0 else None),
             ceil_construct=float(np.sqrt(max(sb_a, 0.0))), ceil_gene=float(np.sqrt(max(sb_b, 0.0))),
             pds_a=float(np.median([float(x["a_pds"]) for x in rows])),
             pds_b=float(np.median([float(x["b_pds"]) for x in rows])))
    P(f"    {name:34s} n={d['n_targets']:2d}  (a) {L.fmt(d['a'])}  (b) {L.fmt(d['b'])}  "
      f"rho_guide {d['rho']['point']:+.3f} [{d['rho']['lo']:+.3f}, {d['rho']['hi']:+.3f}]  "
      f"pds (a) {d['pds_a']:.3f} -> (b) {d['pds_b']:.3f}")
    P(f"    {'':34s}    (a2) other construct's halves {d['a2']['point']:+.4f}; disattenuated rho "
      f"{('%+.3f' % d['rho_dis']) if d['rho_dis'] is not None else '   n/a'}; SB400 (a) {sb_a:+.4f} (b) "
      f"{sb_b:+.4f}  ->  ceiling construct-level {d['ceil_construct']:.3f}, GENE-level (zero-shot) "
      f"{d['ceil_gene']:.3f}")
    return d


def main() -> None:
    t0 = time.time()
    P(f"293 strat293.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    out: dict = {}
    for lab, tag in TAGS.items():
        out[lab] = {}
        allrows = list(csv.DictReader(open(X.WORK / f"construct293{tag}.csv")))
        for n in (100, 200):
            rows = [x for x in allrows if int(x["n"]) == n]
            for x in rows:
                x["kd_weak"] = min(abs(float(x["kd_big"])), abs(float(x["kd_small"])))
                x["kd_strong"] = max(abs(float(x["kd_big"])), abs(float(x["kd_small"])))
                x["a2v"] = float(x["a2_cos"]) if x.get("a2_cos") not in (None, "", "nan") else np.nan
            med = float(np.median([x["kd_weak"] for x in rows]))
            P(f"\n  {lab}, n={n}/side  (median |lfc| of the weaker construct {med:.3f})")
            out[lab][n] = dict(
                all=block(rows, "ALL", n),
                both_work=block([x for x in rows if x["kd_weak"] >= 1.0], "both constructs |lfc| >= 1.0", n),
                one_fails=block([x for x in rows if x["kd_weak"] < 1.0], "weaker construct |lfc| < 1.0", n),
                upper_half=block([x for x in rows if x["kd_weak"] >= med], "weaker construct, upper half", n),
                lower_half=block([x for x in rows if x["kd_weak"] < med], "weaker construct, lower half", n),
                a2_measurable=block([x for x in rows if np.isfinite(x["a2v"]) and x["a2v"] >= 0.05],
                                    "BOTH constructs measurable (a2>=0.05)", n),
                a2_unmeasurable=block([x for x in rows if not (np.isfinite(x["a2v"]) and x["a2v"] >= 0.05)],
                                      "2nd construct not measurable (a2<0.05)", n),
                median_kd_weak=med,
                kd_table={x["target"]: dict(kd_big=float(x["kd_big"]), kd_small=float(x["kd_small"]),
                                            kd_weak=x["kd_weak"], a=float(x["a_cos"]), b=float(x["b_cos"]))
                          for x in rows})
    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, X.WORK / "strat293.json")
    P(f"\ndone {out['seconds']}s")
    P.save()


if __name__ == "__main__":
    main()
