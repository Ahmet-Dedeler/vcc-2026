"""243 Phase 3 — is there ONE convention under which ONE factor retrodicts all three board
direction points (a14->a15, cd4f034, h1union) within +-15%?  If not, the smallest model that
does, its falsifier, and the re-price of 237's ARM6 and 228's per-gene arm under it.

Inputs: work/matrix243.json (Phase 1), work/invert243.json (Phase 2), the board rows,
229's price229.json (ARM6 per fold), 228's analyse228.json / price228.json.

    uv run python rule243.py
"""
from __future__ import annotations

import json

import numpy as np

from lib243 import CANCER, COEFF, CTX, FOLDS, ROOT, SUMD, WORK, board_row, compose218

M = json.loads((WORK / "matrix243.json").read_text())
I = json.loads((WORK / "invert243.json").read_text())
SPAN_MSE, SPAN_PDS = 0.950162566, 0.452          # board216.py / 197 (score_mse map), 162
STEPS = ("a14->a15", "cd4f034", "h1union")
COVERAGE = {"a14->a15": 276 / 300, "cd4f034": 276 / 300, "h1union": 77 / 300}   # metas
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def cell(step, functional, tag, folds, agg):
    for c in M["cells"]:
        if (c["step"], c["functional"], c["tag"], c["folds"], c["agg"]) == (step, functional, tag, folds, agg):
            return c["rel"]
    raise KeyError((step, functional, tag, folds, agg))


def main() -> None:
    S = I["steps"]
    board = {
        "floor": {s: S[s]["rel_floor"] for s in STEPS},
        "measured": {s: S[s]["rel_meas"] for s in STEPS},
        "pds_route": {s: S[s]["rel_pds_route"] for s in STEPS},
    }
    P("=== board-side relative cosine gains, three readings ===")
    for k, v in board.items():
        P(f"  {k:10s} " + "  ".join(f"{s} {100*v[s]:+.3f}%" for s in STEPS))

    # ---- the factor table ----------------------------------------------------------------
    P("\n=== factor = board rel / PLOLO rel, every convention; 'k*' = the single factor that "
      "minimises the worst |deviation|, and that worst deviation ===")
    P(f"  {'functional':10s}{'tag':5s}{'folds':6s}{'agg':8s}{'board':10s}"
      f"{'a14->a15':>10s}{'cd4f034':>10s}{'h1union':>10s}{'k*':>8s}{'worst dev':>11s}{'<=15%':>7s}")
    table = []
    for functional in ("cos_agg", "macro"):
        for tag in ("a23", "a15"):
            for folds in ("3c", "4"):
                for agg in ("perfold", "pooled"):
                    for reading in ("floor", "measured", "pds_route"):
                        ks = {s: board[reading][s] / cell(s, functional, tag, folds, agg)
                              for s in STEPS}
                        vals = np.array([ks[s] for s in STEPS])
                        if (vals <= 0).any():
                            kstar, worst = float("nan"), float("nan")
                        else:
                            kstar = float(np.sqrt(vals.max() * vals.min()))   # minimax in log
                            worst = float(max(abs(v / kstar - 1) for v in vals))
                        row = dict(functional=functional, tag=tag, folds=folds, agg=agg,
                                   reading=reading, k=ks, kstar=kstar, worst=worst,
                                   ok=bool(worst <= 0.15))
                        table.append(row)
                        P(f"  {functional:10s}{tag:5s}{folds:6s}{agg:8s}{reading:10s}"
                          + "".join(f"{ks[s]:>10.3f}" for s in STEPS)
                          + f"{kstar:>8.3f}{100*worst:>+10.1f}%{'  YES' if row['ok'] else '   no':>7s}")
    passing = [r for r in table if r["ok"]]
    P(f"\n  conventions passing the +-15% single-factor test: {len(passing)} of {len(table)}")
    for r in passing:
        P(f"    {r['functional']} {r['tag']} {r['folds']} {r['agg']} {r['reading']}: "
          + ", ".join(f"{s} {r['k'][s]:.3f}" for s in STEPS) + f"  k* {r['kstar']:.3f}")
    best = min(table, key=lambda r: (np.nan_to_num(r["worst"], nan=9)))
    P(f"  closest: {best['functional']} {best['tag']} {best['folds']} {best['agg']} "
      f"{best['reading']}: " + ", ".join(f"{s} {best['k'][s]:.3f}" for s in STEPS)
      + f"  k* {best['kstar']:.3f}  worst {100*best['worst']:+.1f}%")

    # ---- the two-class model on the reference convention -----------------------------------
    REF = dict(functional="cos_agg", folds="3c", agg="perfold", reading="measured")
    P("\n=== two-class model on the reference convention (cos_agg, 3 cancer folds, per-fold "
      "mean, MEASURED board cosine) ===")
    classes = {}
    for tag in ("a23", "a15"):
        ks = {s: board["measured"][s] / cell(s, "cos_agg", tag, "3c", "perfold") for s in STEPS}
        k_dose = float(np.sqrt(ks["a14->a15"] * ks["cd4f034"]))
        spread_dose = max(abs(ks[s] / k_dose - 1) for s in ("a14->a15", "cd4f034"))
        k_new = ks["h1union"]
        # coverage-linear alternative through the same two group values
        slope = (k_dose - k_new) / (COVERAGE["cd4f034"] - COVERAGE["h1union"])
        classes[tag] = dict(k=ks, k_dose=k_dose, dose_within=spread_dose, k_newsource=k_new,
                            cov_slope=slope, cov_intercept=k_new - slope * COVERAGE["h1union"])
        P(f"  tag {tag}: k a14->a15 {ks['a14->a15']:.3f}, cd4f034 {ks['cd4f034']:.3f} -> "
          f"k_dose {k_dose:.3f} (both within {100*spread_dose:+.1f}%);  k_newsource {k_new:.3f}")
        P(f"           coverage-linear: k(cov) = {classes[tag]['cov_intercept']:.3f} + "
          f"{slope:.3f} * cov   (same fit, 2 params, 2 distinct coverages)")

    # ---- the pds and DE-four exchange rates, per step (for the re-price) -------------------
    P("\n=== the other legs, realised per step ===")
    rates = {}
    b = {"a14": board_row("a14"), "rc_callset": board_row("rc_callset"),
         "cd4f034": board_row("cd4f034"), "h1union": board_row("h1union")}
    pairs = {"a14->a15": ("a14", "rc_callset"), "cd4f034": ("rc_callset", "cd4f034"),
             "h1union": ("cd4f034", "h1union")}
    P(f"  {'step':10s}{'d pds board':>12s}{'PLOLO 3c':>10s}{'k_pds 3c':>10s}{'PLOLO 4':>10s}"
      f"{'k_pds 4':>9s}{'DE-four ovr':>13s}{'per meas dcos':>14s}{'mse ovr':>9s}{'per dcos':>9s}")
    for s, (x, y) in pairs.items():
        dpb = b[y]["raw"]["pds"] - b[x]["raw"]["pds"]
        # a15's scaled members are the ledger's; a14->a15 uses rc_callset's row as a15's
        # profile row for pds/mse (identical to 16 s.f.) but NOT for the DE members, which
        # rc_hybrid/rc_callset moved.  For a14->a15's DE-four use 219's value (its own ledger
        # arithmetic on a14/a15 scaled members).
        if s == "a14->a15":
            de4 = 0.0008791666666666658                  # 219 phase0_219.json propagation
            mse_ovr = 0.0015736666666666653
        else:
            de4 = sum(b[y]["scaled"][m] - b[x]["scaled"][m] for m in ("nmae", "fid", "reach", "jac")) / 6
            mse_ovr = (b[y]["scaled"]["mse"] - b[x]["scaled"]["mse"]) / 6
        p3 = cell(s, "pds", "a23", "3c", "mean")
        p4 = cell(s, "pds", "a23", "4", "mean")
        dcm = S[s]["d_cos_meas"]
        rates[s] = dict(d_pds_board=dpb, plolo_3c=p3, plolo_4=p4, k_pds_3c=dpb / p3, k_pds_4=dpb / p4,
                        de4_ovr=de4, de4_per_dcos=de4 / dcm, mse_ovr=mse_ovr, mse_per_dcos=mse_ovr / dcm)
        P(f"  {s:10s}{dpb:>+12.6f}{p3:>+10.6f}{dpb/p3:>10.3f}{p4:>+10.6f}{dpb/p4:>9.3f}"
          f"{de4:>+13.6f}{de4/dcm:>14.3f}{mse_ovr:>+9.5f}{mse_per_dcos:>9.3f}" if False else
          f"  {s:10s}{dpb:>+12.6f}{p3:>+10.6f}{dpb/p3:>10.3f}{p4:>+10.6f}{dpb/p4:>9.3f}"
          f"{de4:>+13.6f}{de4/dcm:>14.3f}{mse_ovr:>+9.5f}{mse_ovr/dcm:>9.3f}")

    # ---- the re-price ----------------------------------------------------------------------
    inc = I["arms"]["h1union"]
    cos_inc = inc["inv"]["cos"]
    mse_slope = 2.0 * inc["inv"]["mean_sqrtA"] / np.sqrt(SUMD)     # d raw / d cos, exact identity
    P(f"\n=== re-price against the incumbent h1union: measured cos {cos_inc:.6f}, "
      f"d raw mse / d cos = -{mse_slope:.4f} (floor line: -{COEFF:.4f}) ===")

    # ARM6 at g=0.35 (229 price229.json, per fold, against h1union's per-fold base in 218's
    # compose218.json -- the same object 229's KA-5 reproduces)
    p229 = json.loads((ROOT / "experiments/229-stack4-build/work/price229.json").read_text())
    c218 = compose218()
    arm6 = {}
    for a in p229["arms"]:
        if a["arm"] == "ARM6" and abs(a["g"] - 0.35) < 1e-9:
            base = c218[(a["fold"], "cd40.34+h1PAIRCARRY0.28")]
            arm6[a["fold"]] = dict(rel_agg=a["cos_agg"] / base["cos_agg"] - 1.0,
                                   rel_macro=a["rel_cos_pct_vs_h1union"] / 100.0,
                                   d_pds=a["d_pds_vs_h1union"], carrying=a["carrying_rows"],
                                   n_rows=a["n_rows"])
    arm6_agg3 = float(np.mean([arm6[f]["rel_agg"] for f in CANCER]))
    arm6_mac3 = float(np.mean([arm6[f]["rel_macro"] for f in CANCER]))
    arm6_pds3 = float(np.mean([arm6[f]["d_pds"] for f in CANCER]))
    ka229 = abs(100 * arm6_mac3 - 2.865) < 5e-3
    P(f"  ARM6 g=0.35 PLOLO (3 cancer, per-fold mean): rel cos_agg {100*arm6_agg3:+.3f}%, "
      f"rel macro {100*arm6_mac3:+.3f}% (229 §6.1: +2.865% -> {'PASS' if ka229 else 'FAIL'}), "
      f"d pds {arm6_pds3:+.6f} (229: +0.00269); coverage 176/300")
    assert ka229

    # 228 gain_reg@0.25 (analyse228.json per_variant, production alpha; price228.json @adel)
    a228 = json.loads((ROOT / "experiments/228-theta-gain-plolo/work/analyse228.json").read_text())
    v228 = a228["per_variant"]["gain_reg@0.25"]
    p228 = json.loads((ROOT / "experiments/228-theta-gain-plolo/work/price228.json").read_text())["gain_reg@0.25|adel"]
    P(f"  228 gain_reg@0.25 PLOLO (3 cancer): rel cos_agg {v228['cancer_rel_cos_agg']:+.3f}% at "
      f"production alpha (228 §3: +4.328%), {p228['rel_cos_agg_pct']:+.3f}% @adel (228 §7: +4.54%); "
      f"d pds {v228['cancer_d_pds']:+.6f} / {p228['d_pds']:+.6f}; 228's own o_nmae(wald89) "
      f"{p228['o_nmae']['wald89']:+.5f}; coverage 300/300 (occupied rows re-shaped)")

    def price(label, rel_plolo, d_pds_plolo, k_cos, k_pds, de4_rate, note=""):
        dcos = k_cos * rel_plolo * cos_inc
        d_raw = -mse_slope * dcos
        o_mse = -d_raw / SPAN_MSE / 6.0
        o_pds = k_pds * d_pds_plolo / SPAN_PDS / 6.0
        o_de4 = de4_rate * dcos
        tot = o_mse + o_pds + o_de4
        P(f"    {label:44s} k_cos {k_cos:.2f} k_pds {k_pds:.2f} DE4 {de4_rate:.2f}: "
          f"dcos {dcos:+.6f}  raw mse {inc['board']['raw']['mse']+d_raw:.5f}  "
          f"o_mse {o_mse:+.5f}  o_pds {o_pds:+.5f}  o_DE4 {o_de4:+.5f}  TOTAL {tot:+.5f} {note}")
        return dict(label=label, k_cos=k_cos, k_pds=k_pds, de4_rate=de4_rate, dcos=dcos,
                    raw_mse=inc["board"]["raw"]["mse"] + d_raw, o_mse=o_mse, o_pds=o_pds,
                    o_de4=o_de4, total=tot, note=note)

    kn = classes["a23"]["k_newsource"]
    kn15 = classes["a15"]["k_newsource"]
    kd = classes["a23"]["k_dose"]
    kcov6 = classes["a23"]["cov_intercept"] + classes["a23"]["cov_slope"] * (176 / 300)
    de4_new, de4_lo, de4_hi = rates["h1union"]["de4_per_dcos"], rates["a14->a15"]["de4_per_dcos"], rates["cd4f034"]["de4_per_dcos"]
    kp_new, kp_lo, kp_hi = rates["h1union"]["k_pds_3c"], min(r["k_pds_3c"] for r in rates.values()), max(r["k_pds_3c"] for r in rates.values())
    P("\n  ARM6 (g=0.35), the direction legs under each reading of the rule:")
    pr = {"ARM6": [
        price("new-source class, point (h1union's own rates)", arm6_agg3, arm6_pds3, kn, kp_new, de4_new),
        price("new-source class, a15-tag k", arm6_agg3, arm6_pds3, kn15, kp_new, de4_new),
        price("new-source class, LOW band", arm6_agg3, arm6_pds3, min(kn, kn15), kp_lo, de4_lo),
        price("new-source class, HIGH band", arm6_agg3, arm6_pds3, max(kn, kn15), kp_hi, de4_hi),
        price("coverage-linear model (cov 0.587)", arm6_agg3, arm6_pds3, kcov6, kp_new, de4_new),
        price("dose class (if the union point was the outlier)", arm6_agg3, arm6_pds3, kd, kp_new, de4_hi),
    ]}
    P(f"  229's six-member measured (its own nmae identity leg -0.004783): -0.00219;  "
      f"237's bins: <= -0.0015 gate validated / (-0.0015, +0.0015) inconclusive / >= +0.0015 gate wrong")

    P("\n  228 gain_reg@0.25, the direction legs (mse + pds) under each class; the DE-four leg is "
      "NOT carried at the sign-propagation rate for this arm (228 §4-6: it moves |t|-weighted "
      "sign accuracy by +0.00027 only; its nmae is a magnitude effect the rule does not price):")
    kp_dose = float(np.sqrt(rates["a14->a15"]["k_pds_3c"] * rates["cd4f034"]["k_pds_3c"]))
    pr["228"] = [
        price("occupied-rows (dose) class, production alpha", v228["cancer_rel_cos_agg"] / 100, v228["cancer_d_pds"], kd, kp_dose, 0.0),
        price("occupied-rows (dose) class, @adel", p228["rel_cos_agg_pct"] / 100, p228["d_pds"], kd, kp_dose, 0.0),
        price("new-source class, production alpha", v228["cancer_rel_cos_agg"] / 100, v228["cancer_d_pds"], kn, kp_new, 0.0),
        price("coverage-linear model (cov 1.0)", v228["cancer_rel_cos_agg"] / 100, v228["cancer_d_pds"],
              classes["a23"]["cov_intercept"] + classes["a23"]["cov_slope"] * 1.0, kp_dose, 0.0),
    ]
    P(f"  228's own price of the same two legs at k_cos 1.384: o_pds {p228['o_pds']:+.5f} + o_mse "
      f"{p228['o_mse']['k1.384']:+.5f} = {p228['o_pds']+p228['o_mse']['k1.384']:+.5f}; its nmae leg "
      f"{p228['o_nmae']['wald89']:+.5f} (identity instrument, sign unvalidated) decides the arm")

    out = dict(board=board, table=table, passing=[r for r in table if r["ok"]], closest=best,
               classes=classes, rates=rates, coverage=COVERAGE, incumbent_cos=cos_inc,
               mse_slope=mse_slope, arm6=dict(per_fold=arm6, rel_agg_3c=arm6_agg3,
                                              rel_macro_3c=arm6_mac3, d_pds_3c=arm6_pds3),
               reprice=pr)
    (WORK / "rule243.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "rule243.log").write_text("\n".join(LINES))
    P(f"\nwrote {WORK/'rule243.json'}")


if __name__ == "__main__":
    main()
