"""216 — digest the library-form re-test and price the three forms.

Known-answer rows carried inside this analysis, both from 207:
  * `lib0` must reproduce `f0` (CD4 OFF) exactly -- the library weight is a no-op at w = 0;
  * `lib0.5@a0` on `orion_hct116` must reproduce 207 §5.5's published `cos_agg` 0.024304,
    the one fold on which 207 found the library form NOT bad.

Pricing.  The `@adel` policy holds the DELIVERED median row norm at the incumbent's, so and
only so is the arm comparable to a residual arm through 096/197's `X`-only identity.  At
`@a0` and `@apre` the delivered amplitude is 0.26-0.60x and 211's standing rule applies --
`claim` sits at 100.0000% of its cap in the binding context, so **no amplitude credit may be
taken**, and those rows are reported without an `mse` price at all.

    uv run python formanalyse216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
SPAN_PDS, CALIB, SIGMA_PDS = 0.452, 1.011, 0.001494
OURS_COS = 0.079249
SUM_D, KC, CLAIM = 0.22914103411323902, 0.0054102686242422674, 0.24614814267726884
A_TERM, B_TERM = 0.2331583274832182, 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089


def mse_at(c):
    q, eps = CLAIM / SUM_D, A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw):
    return float(min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0))


def main() -> None:
    rows = json.loads((WORK / "form216.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in rows}
    sw = {(r["fold"], r["arm"]): r for r in json.loads((WORK / "sweep216.json").read_text())}
    price = json.loads((WORK / "price216.json").read_text())
    k_cos = price["k_cos"]
    ship_cos = price["shipped_cos_agg"]
    ship_pds = float(np.mean([sw[(h, "f0.1")]["pds"] for h in FOLDS]))
    inc_cos = float(np.mean([sw[(h, "f0.34")]["cos_agg"] for h in FOLDS]))
    inc_pds = float(np.mean([sw[(h, "f0.34")]["pds"] for h in FOLDS]))

    arms = sorted({r["arm"] for r in rows},
                  key=lambda a: (a.split("@")[1], float(a.split("@")[0]
                                                        .replace("lib", "").split("+")[0]),
                                 a))
    have = [a for a in arms if all((h, a) in by for h in FOLDS)]

    def m(arm, key):
        return float(np.mean([by[(h, arm)][key] for h in FOLDS]))

    # ------------------------------------------------------------- known-answer rows
    gates = []
    if all((h, "lib0@a0") in by for h in FOLDS):
        d = max(abs(by[(h, "lib0@a0")]["cos_agg"] - sw[(h, "f0")]["cos_agg"])
                for h in FOLDS)
        gates.append(dict(name="KA lib0@a0 == f0 (CD4 off)", measured=f"{d:.3e}",
                          passed=bool(d < 1e-12)))
    if ("orion_hct116", "lib0.5@a0") in by:
        v = by[("orion_hct116", "lib0.5@a0")]["cos_agg"]
        gates.append(dict(name="KA 207 §5.5 lib0.5 on orion_hct116 cos_agg",
                          expected=0.024304, measured=round(v, 6),
                          passed=bool(abs(v - 0.024304) < 1e-5)))
    for g in gates:
        print(f"[{'PASS' if g['passed'] else 'FAIL'}] {g['name']}: {g.get('measured')}"
              + (f" (expected {g['expected']})" if "expected" in g else ""))

    # ------------------------------------------------------------------ the table
    print(f"\n{'arm':<24}{'alpha':>8}{'dlv/inc':>9}{'cos_agg':>10}{'vs ship%':>10}"
          f"{'vs f0.34%':>11}{'pds':>10}{'d_pds/sig':>11}{'ovr(k=1)':>10}{'class':>10}")
    print("-" * 113)
    out = []
    for arm in have:
        a = m(arm, "alpha")
        dn = m(arm, "delivered_ratio")
        c = m(arm, "cos_agg")
        p = m(arm, "pds")
        amp_ok = abs(dn - 1.0) < 0.02
        rel_ship = (c - ship_cos) / ship_cos
        o_pds = (p - ship_pds) / SPAN_PDS / 6.0 * CALIB
        if amp_ok:
            bc = OURS_COS * (1.0 + rel_ship)
            o_mse = (score_mse(mse_at(bc)) - score_mse(mse_at(OURS_COS))) / 6.0
            total = o_pds + o_mse
            cls = "shape"
        else:
            o_mse, total, cls = None, None, "AMPLITUDE"
        rec = dict(arm=arm, alpha=a, delivered_ratio=dn, cos_agg=c, pds=p,
                   rel_vs_shipped=rel_ship, rel_vs_incumbent=(c - inc_cos) / inc_cos,
                   d_pds_vs_shipped=p - ship_pds,
                   pds_sigma_vs_incumbent=(p - inc_pds) / SIGMA_PDS,
                   o_pds=o_pds, o_mse=o_mse, total_k1=total,
                   amplitude_matched=bool(amp_ok), arm_class=cls,
                   per_fold_cos={h: by[(h, arm)]["cos_agg"] for h in FOLDS})
        out.append(rec)
        tot = f"{total:>+10.5f}" if total is not None else f"{'--':>10}"
        print(f"{arm:<24}{a:>8.4f}{dn:>9.4f}{c:>10.6f}{100*rel_ship:>+10.2f}"
              f"{100*rec['rel_vs_incumbent']:>+11.2f}{p:>10.6f}"
              f"{rec['pds_sigma_vs_incumbent']:>11.2f}{tot}{cls:>10}")

    matched = [r for r in out if r["amplitude_matched"]]
    best = max(matched, key=lambda r: r["total_k1"]) if matched else None
    lib1 = [r for r in out if r["arm"] == "lib1@adel"]
    summary = dict(
        k_cos=k_cos, shipped_cos=ship_cos, incumbent_cos=inc_cos,
        best_amplitude_matched=best,
        one_over_k_arm=lib1[0] if lib1 else None,
        residual_incumbent=dict(cos_agg=inc_cos, pds=inc_pds,
                                total_k1=[r for r in price["rows"]
                                          if r["f"] == 0.34][0]["total_k1"]),
        gates=gates, gate_shape=GATE_SHAPE, gate_amplitude=GATE_AMP)
    if best:
        print(f"\nbest amplitude-matched library arm: {best['arm']}  "
              f"{best['total_k1']:+.5f} against the residual incumbent's "
              f"{summary['residual_incumbent']['total_k1']:+.5f} "
              f"(both vs the shipped f = 0.10)")
    if lib1:
        r = lib1[0]
        print(f"the 1/k route (`lib1@adel`, CD4 as one more EQUAL line): cos_agg "
              f"{r['cos_agg']:.6f} = {100*r['rel_vs_incumbent']:+.2f}% against the "
              f"incumbent residual, alpha {r['alpha']:.4f}")
    (WORK / "formanalyse216.json").write_text(json.dumps(
        dict(rows=out, **summary), indent=2))


if __name__ == "__main__":
    main()
