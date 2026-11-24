"""216 — digest the operator refit and price it against the incumbent.

Reference arms (`res0`, `res0.10`, `res0.34`) are re-scored through the SHIPPED operator
inside `op216.py` itself and carried in each row, so every refit number has its own
same-fold, same-code baseline and no cross-script comparison is needed.

Pricing uses the same two routes as everywhere else: `pds` through PLOLO x1.011, `mse`
through 096/197's `X`-only identity at matched delivered amplitude.  `nmae` is NOT priced
here -- if an arm survives this table it has to go through 188's exact identity before it
can be believed, and §6.0b measured that amplifying a library direction is exactly what
destroys `nmae`.

    uv run python opanalyse216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
SPAN_PDS, CALIB = 0.452, 1.011
SIGMA_PDS = 0.001494
OURS_COS = 0.079249
SUM_D, KC, CLAIM = 0.22914103411323902, 0.0054102686242422674, 0.24614814267726884
A_TERM, B_TERM = 0.2331583274832182, 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
ALPHA0 = 0.982


def mse_at(c):
    q, eps = CLAIM / SUM_D, A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw):
    return float(min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0))


def main() -> None:
    rows = json.loads((WORK / "op216.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in rows}
    sw = {(r["fold"], r["arm"]): r for r in json.loads((WORK / "sweep216.json").read_text())}
    have_folds = [f for f in FOLDS if any(k[0] == f for k in by)]
    arms = sorted({r["arm"] for r in rows})
    complete = [a for a in arms if all((f, a) in by for f in have_folds)]

    # the incumbent, through the SHIPPED operator, on the same folds
    inc_cos = float(np.mean([sw[(f, "f0.34")]["cos_agg"] for f in have_folds]))
    inc_pds = float(np.mean([sw[(f, "f0.34")]["pds"] for f in have_folds]))
    ship_cos = float(np.mean([sw[(f, "f0.1")]["cos_agg"] for f in have_folds]))
    ship_pds = float(np.mean([sw[(f, "f0.1")]["pds"] for f in have_folds]))

    print(f"folds present: {have_folds}")
    print(f"incumbent (shipped operator, res f=0.34) on these folds: "
          f"cos_agg {inc_cos:.6f}  pds {inc_pds:.6f}\n")
    print(f"{'arm':<30}{'beta':>9}{'alpha':>8}{'dlv':>8}{'cos_agg':>10}{'vs inc%':>9}"
          f"{'pds':>10}{'d_pds/sig':>11}{'ovr pds':>10}{'ovr mse':>10}{'TOTAL':>10}"
          f"{'class':>11}")
    print("-" * 136)
    out = []
    for a in complete:
        c = float(np.mean([by[(f, a)]["cos_agg"] for f in have_folds]))
        p = float(np.mean([by[(f, a)]["pds"] for f in have_folds]))
        al = float(np.mean([by[(f, a)]["alpha"] for f in have_folds]))
        dn = float(np.mean([by[(f, a)]["delivered_ratio"] for f in have_folds]))
        beta = float(np.mean([by[(f, a)]["spectrum"]["beta"] for f in have_folds]))
        rel_ship = (c - ship_cos) / ship_cos
        o_pds = (p - ship_pds) / SPAN_PDS / 6.0 * CALIB
        bc = OURS_COS * (1.0 + rel_ship)
        o_mse = (score_mse(mse_at(bc)) - score_mse(mse_at(OURS_COS))) / 6.0
        tot = o_pds + o_mse
        amp = abs(al / ALPHA0 - 1.0) > 0.02 or abs(dn - 1.0) > 0.02
        out.append(dict(arm=a, beta=beta, alpha=al, delivered=dn, cos_agg=c, pds=p,
                        rel_vs_incumbent=(c - inc_cos) / inc_cos,
                        pds_sigma_vs_incumbent=(p - inc_pds) / SIGMA_PDS,
                        o_pds=o_pds, o_mse=o_mse, total_k1=tot,
                        arm_class="AMPLITUDE" if amp else "shape",
                        per_fold={f: dict(cos_agg=by[(f, a)]["cos_agg"],
                                          pds=by[(f, a)]["pds"],
                                          alpha=by[(f, a)]["alpha"]) for f in have_folds},
                        ref_cos=by[(have_folds[0], a)].get("ref_cos_agg")))
        print(f"{a:<30}{beta:>9.5f}{al:>8.4f}{dn:>8.4f}{c:>10.6f}"
              f"{100*out[-1]['rel_vs_incumbent']:>+9.2f}{p:>10.6f}"
              f"{out[-1]['pds_sigma_vs_incumbent']:>+11.2f}{o_pds:>+10.5f}{o_mse:>+10.5f}"
              f"{tot:>+10.5f}{out[-1]['arm_class']:>11}")

    inc_tot = None
    price = json.loads((WORK / "price216.json").read_text())
    inc_tot = [r for r in price["rows"] if r["f"] == 0.34][0]["total_k1"]
    print(f"\nincumbent's own priced total on the full 4-fold set (pds+mse route): "
          f"{inc_tot:+.5f}")
    print(f"\n{'arm':<30}{'d vs incumbent':>16}{'gate':>9}{'x gate':>9}  verdict")
    print("-" * 78)
    verdicts = []
    for r in out:
        d = r["total_k1"] - inc_tot
        gate = GATE_AMP if r["arm_class"] == "AMPLITUDE" else GATE_SHAPE
        v = "CANDIDATE - price nmae" if d > gate else "VETO"
        verdicts.append(dict(arm=r["arm"], delta=d, gate=gate, verdict=v))
        print(f"{r['arm']:<30}{d:>+16.5f}{gate:>9.4f}{d/gate:>9.2f}  {v}")

    best = max(verdicts, key=lambda v: v["delta"])
    (WORK / "opanalyse216.json").write_text(json.dumps(
        dict(rows=out, verdicts=verdicts, folds=have_folds, incumbent_total=inc_tot,
             incumbent_cos=inc_cos, incumbent_pds=inc_pds, best=best,
             note="pds+mse route only; nmae unpriced"), indent=2))
    print(f"\nbest: {best['arm']}  {best['delta']:+.5f} against {best['gate']:+.4f}")


if __name__ == "__main__":
    main()
