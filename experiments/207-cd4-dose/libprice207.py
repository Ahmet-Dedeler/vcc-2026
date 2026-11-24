"""207 Phase 2/3 — price the FORM sweep (library, residual, composed) against the incumbent.

Baseline for every delta here is the SHIPPED configuration: `f0.1` in `sweep207.json`
(post-whitening residual at 099's `f = 0.10`, no library member, alpha 0.982).

Every library arm is priced AT ITS OWN RE-SOLVED ALPHA (`@a*`), because adding a fourth member
to the consensus lowers the pre-whitening median row norm (measured: x0.84 at w = 0.5) and the
operator preserves that norm, so the delivered amplitude would otherwise move.  The cost of the
re-solve is reported separately, exactly as 203 had to.

    uv run python libprice207.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from price207 import (CALIB, FOLDS, SPAN_NMAE, SPAN_PDS, W_HI, OURS_COS,  # noqa: E402
                      mse_at, score_mse)

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
NB = 4000


def main() -> None:
    by = {(r["fold"], r["arm"]): r
          for r in json.loads((WORK / "sweep207.json").read_text())}
    by.update({(r["fold"], r["arm"]): r
               for r in json.loads((WORK / "library207.json").read_text())})
    k_cos = json.loads((WORK / "price207.json").read_text())["k_cos"]
    base = "f0.1"

    def m(arm, key):
        return float(np.mean([by[(f, arm)][key] for f in FOLDS]))

    def boot(arm, key):
        rng = np.random.default_rng(20260903)
        d = []
        for f in FOLDS:
            a = np.asarray(by[(f, arm)][key])
            b = np.asarray(by[(f, base)][key])
            idx = rng.integers(0, len(a), (NB, len(a)))
            d.append(np.array([(a[i] - b[i]).mean() for i in idx]))
        d = np.mean(d, axis=0)
        return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float((d > 0).mean())

    cos_ship = m(base, "cos_agg")
    arms = [a for (f, a) in by if f == FOLDS[0] and (a.startswith("lib"))]
    arms = [a for a in arms if all((f, a) in by for f in FOLDS)]
    arms.sort()

    rows = []
    for arm in arms + ["f0", base, "f0.2", "f0.4"]:
        if not all((f, arm) in by for f in FOLDS):
            continue
        d_pds = m(arm, "pds") - m(base, "pds")
        rel = (m(arm, "cos_agg") - cos_ship) / cos_ship
        d_qw = m(arm, "a15_sgn_w") - m(base, "a15_sgn_w")
        o_pds = d_pds / SPAN_PDS / 6.0 * CALIB
        o_nmae = (2.0 * W_HI * d_qw) / SPAN_NMAE / 6.0
        out = {}
        for tag, k in (("k1", 1.0), ("kcal", k_cos)):
            bc = OURS_COS * (1.0 + k * rel)
            out[tag] = (score_mse(mse_at(bc)) - score_mse(mse_at(OURS_COS))) / 6.0
        lo, hi, p = boot(arm, "_cos_per_target_a23")
        lop, hip, pp = boot(arm, "_pds_per_target_a23")
        a = by[(FOLDS[0], arm)]["audit"]
        rows.append(dict(arm=arm, w=a.get("lib_w"), f=a.get("residual_f", a.get("effective_f")),
                         alpha=a.get("alpha"), prenorm=a.get("pre_norm_ratio"),
                         cos_agg=m(arm, "cos_agg"), rel=rel, pds=m(arm, "pds"),
                         d_pds=d_pds, d_pds_ci=[lop, hip], d_pds_p=pp,
                         d_cos=m(arm, "cos_agg") - cos_ship, d_cos_ci=[lo, hi], d_cos_p=p,
                         q_w=m(arm, "a15_sgn_w"), d_qw=d_qw,
                         o_pds=o_pds, o_nmae=o_nmae, o_mse_k1=out["k1"],
                         o_mse_kcal=out["kcal"],
                         total_k1=o_pds + o_nmae + out["k1"],
                         total_kcal=o_pds + o_nmae + out["kcal"],
                         rownorm=m(arm, "median_row_norm_wcols")))

    rows.sort(key=lambda r: -r["total_k1"])
    print(f"{'arm':<22}{'alpha':>7}{'prenorm':>8}{'cos_agg':>10}{'Δcos':>10}"
          f"{'Δcos 95% CI':>24}{'P>0':>6}{'Δpds':>10}{'P>0':>6}"
          f"{'ovr pds':>9}{'ovr mse':>9}{'TOTAL':>9}{'TOT cal':>9}")
    print("-" * 139)
    for r in rows:
        print(f"{r['arm']:<22}{(r['alpha'] or 0.982):>7.3f}{(r['prenorm'] or 1):>8.4f}"
              f"{r['cos_agg']:>10.6f}{r['d_cos']:>+10.6f}"
              f"  [{r['d_cos_ci'][0]:+.6f},{r['d_cos_ci'][1]:+.6f}]{r['d_cos_p']:>6.3f}"
              f"{r['d_pds']:>+10.6f}{r['d_pds_p']:>6.3f}"
              f"{r['o_pds']:>+9.5f}{r['o_mse_k1']:>+9.5f}{r['total_k1']:>+9.5f}"
              f"{r['total_kcal']:>+9.5f}")
    (WORK / "libprice207.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {WORK/'libprice207.json'}")


if __name__ == "__main__":
    main()
