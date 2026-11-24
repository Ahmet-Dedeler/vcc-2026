"""230 Phase 2 — price 208's Hamilton rider under every toll law the data supports.

Three things this file adds to `fit230.py`:

  1. THE DIRECTION OF HAMILTON'S PROFILE CHANGE.  Every other arm in the digest-differing class
     changes the INTENDED profile.  Hamilton changes only its integer REALISATION, and 196/208
     measured that it moves the delivered profile TOWARD the intended one.  Measured here as a
     distance to `x = n L p` under each rounding rule, so "0.03% tilt correction" stops being
     a description and becomes a number with the same units as the arms it is compared to.
  2. THE WITHIN-CLASS FIT.  Rule 13's STEP is a claim about the digest-differing class.  Three
     of its members are now scored (a23, `h1union`, `cd4f034`), so the class can be fit on its
     own instead of borrowing two points that turned out not to belong to it.
  3. THE CELL-CHANNEL ALTERNATIVE.  The axis a21/a22 actually moved, and what it charges
     Hamilton.

    uv run python phase2_230.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
sys.path.insert(0, str(HERE))
import recon196                                                        # noqa: E402
from magnitude230 import hamilton                                      # noqa: E402

SCALE_REACH = 1.1292          # BOARD-LEDGER: scaled = 1.1292 * (raw - 0.0790)
COMPOSITE_219 = 0.002203      # 219 s2.2, 214's shipped renderer, KA4b
SIX_219 = 0.003904            # 219 s2.2 six-member, KA4d
FLOOR = 0.0019                # rule 13's shape-only floor


def main() -> None:
    F = json.loads((WORK / "fit230.json").read_text())
    tab = {t["arm"]: t for t in F["table"]}
    ham = F["ham"]

    # ---- 0. which CELL of rule 13's table does each scored arm actually fall in? --------
    # The table's first row has TWO conjunctive conditions: the arm must be an EMISSION-LAW
    # change AND its per-gene panel-total digest must differ.  `emission_law` below is True
    # only if the arm changes what `emit_counts` does (its mode, spread, phi or the |z|
    # re-assignment); a change to the EFFECT, or to how the totals are rounded, is not one.
    EMISSION_LAW = {"a8": True, "a10": True, "a21": True, "a22": True, "rc_hybrid": True,
                    "rc_callset": True, "a23": False, "cd4f034": False, "h1union": False}
    print("=== 0. RULE 13's TABLE, APPLIED TO EVERY SCORED ARM ===\n")
    print(f"{'arm':<12}{'emission law?':>15}{'digest differs?':>17}   {'rule-13 row':<28}"
          f"{'toll charged':>14}{'REALISED':>11}")
    cells: dict[str, list[float]] = {}
    for a, t in tab.items():
        el, dd = EMISSION_LAW[a], t["digest_same"] == 0
        rowname = ("emission, profile-CHANGING" if (el and dd) else
                   "emission, profile-PRESERVING" if el else "profile-class")
        charged = -0.0146 if (el and dd) else 0.0
        cells.setdefault(rowname, []).append(t["d_reach"])
        print(f"{a:<12}{str(el):>15}{str(dd):>17}   {rowname:<28}{charged:>+14.4f}"
              f"{t['d_reach']:>+11.6f}")
    print()
    for rowname in ("emission, profile-CHANGING", "emission, profile-PRESERVING",
                    "profile-class"):
        v = cells.get(rowname, [])
        if not v:
            print(f"  {rowname:<30} n = 0   *** THE ROW THE -0.0146 IS CHARGED TO IS "
                  f"EMPTY ***")
        else:
            print(f"  {rowname:<30} n = {len(v)}   mean {np.mean(v):+.6f}   "
                  f"range [{min(v):+.6f}, {max(v):+.6f}]")
    ys = np.array([t["d_reach"] for t in tab.values()])
    best_loo = min(F["fits"][k]["loo_rmse"] for k in F["fits"])
    print(f"\n  sd of board d_reach over all {ys.size} arms: {ys.std(ddof=1):.6f}. "
          f"Best LOO RMSE of any\n  magnitude law: {best_loo:.6f} -- "
          f"{best_loo/ys.std(ddof=1):.2f}x the sd, i.e. every fitted law predicts a held-out\n"
          f"  arm WORSE than predicting the mean does.\n")

    # ---- 1. Hamilton's change is a CORRECTION, and here is its size --------------------
    print("=== 1. DISTANCE TO THE INTENDED PROFILE `x = n L p`, per rounding rule ===\n")
    rows = {"floor": [], "ham": []}
    for ctx in ("A", "B", "C"):
        D = recon196.load_chal(ctx)
        for i in range(0, D["sums"].shape[0], 5):            # every 5th block, 180 total
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=400, seed=5500 + i)
            x = R.x
            s = float(x.sum())
            for rule, T in (("floor", np.floor(x)), ("ham", hamilton(x).astype(float))):
                rows[rule].append(float(np.abs(T - x).sum() / s))
    fl, hm_ = float(np.mean(rows["floor"])), float(np.mean(rows["ham"]))
    print(f"{'rule':<10}{'relL1 vs INTENDED':>20}{'208 tilt (log2)':>18}")
    print(f"{'floor':<10}{fl:>20.3e}{-0.012511:>+18.6f}   <- a15 / the incumbent")
    print(f"{'Hamilton':<10}{hm_:>20.3e}{+0.000339:>+18.6f}   <- 208 rider A")
    print(f"\n  Hamilton's delivered profile is {fl/hm_:.2f}x CLOSER to the intended one, and "
          f"its\n  size tilt is {abs(0.012511/0.000339):.0f}x smaller. Its 'profile change' is a "
          f"REPAIR of the\n  incumbent's rounding error, not a new perturbation. Every other arm "
          f"in the\n  digest-differing class changes the intended profile itself.")

    # ---- 2. the within-class fit -------------------------------------------------------
    cls = ["a23", "h1union", "cd4f034"]
    m = np.array([tab[a]["relL1"] for a in cls])
    y = np.array([tab[a]["d_reach"] for a in cls])
    A = np.column_stack([np.ones_like(m), m])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred_ham_lin = float(beta[0] + beta[1] * ham["relL1"])
    r2 = 1.0 - float(((y - A @ beta) ** 2).sum()) / float(((y - y.mean()) ** 2).sum())
    print(f"\n=== 2. THE DIGEST-DIFFERING CLASS, FIT ON ITS OWN MEMBERS (n = 3) ===\n")
    print(f"{'arm':<12}{'relL1':>11}{'board d_reach':>15}")
    for a in cls:
        print(f"{a:<12}{tab[a]['relL1']:>11.3e}{tab[a]['d_reach']:>+15.6f}")
    print(f"{'ham':<12}{ham['relL1']:>11.3e}{'?':>15}")
    print(f"\n  LINEAR within class: d_reach = {beta[0]:+.6f} {beta[1]:+.4f} * relL1   "
          f"R^2 = {r2:.4f}")
    print(f"  the slope is POSITIVE: inside this class, the BIGGER the profile change the "
          f"BETTER\n  `reach` did. Rule 13 assumes the opposite sign.")
    print(f"  predicted toll at Hamilton's magnitude: {pred_ham_lin:+.6f} scaled "
          f"= {pred_ham_lin/6:+.6f} overall")

    # ---- 3. the cell channel ------------------------------------------------------------
    allarms = list(tab)
    cv = np.array([tab[a]["lib_cv"] for a in allarms])
    ya = np.array([tab[a]["d_reach"] for a in allarms])
    Ac = np.column_stack([np.ones_like(cv), cv])
    bc, *_ = np.linalg.lstsq(Ac, ya, rcond=None)
    r2c = 1.0 - float(((ya - Ac @ bc) ** 2).sum()) / float(((ya - ya.mean()) ** 2).sum())
    loo = []
    for k in range(len(ya)):
        keep = np.arange(len(ya)) != k
        bk, *_ = np.linalg.lstsq(Ac[keep], ya[keep], rcond=None)
        loo.append(ya[k] - Ac[k] @ bk)
    print(f"\n=== 3. THE AXIS a21/a22 ACTUALLY MOVED — emitted per-cell library CV ===\n")
    print(f"  LINEAR on lib CV over all {len(ya)} arms: d_reach = {bc[0]:+.6f} "
          f"{bc[1]:+.4f} * CV   R^2 = {r2c:.4f}   LOO RMSE {np.sqrt(np.mean(np.square(loo))):.6f}")
    print(f"  best relL1 law R^2 was {max(F['fits'][k]['r2'] for k in F['fits']):.4f}.")
    print(f"  Hamilton's emitted library CV is the incumbent's ({tab['cd4f034']['lib_cv']:.4f}) "
          f"by construction --\n  the rounding rule is UPSTREAM of the emitter's cell law -- so "
          f"this axis charges it {bc[1]*0.0:+.6f}.")

    # ---- 4. the price -------------------------------------------------------------------
    zsd = F["falsifier2"]["zero_sd"]
    nzm, nzsd = F["falsifier2"]["nonzero_mean"], F["falsifier2"]["nonzero_sd"]
    laws = [
        ("rule 13 as written (STEP, -0.0146)", -0.0146,
         "mean(a21,a22) -- BOTH digest-IDENTICAL, so calibrated outside the class"),
        ("most punitive in-class point (a23)", -0.007933,
         "the worst realised outcome of any digest-differing arm"),
        ("digest-differing CLASS MEAN", nzm,
         f"a23 / h1union / cd4f034, sd {nzsd:.6f}"),
        ("within-class LINEAR at ham's relL1", pred_ham_lin, "phase2 section 2, n=3, R^2 "
         f"{r2:.3f}"),
        ("STEP fit on all 9 board points", F["fits"]["relL1|STEP"]["toll_at_ham"],
         f"R^2 {F['fits']['relL1|STEP']['r2']:.4f}"),
        ("LINEAR fit on all 9 board points", F["fits"]["relL1|LINEAR"]["toll_at_ham"],
         f"R^2 {F['fits']['relL1|LINEAR']['r2']:.4f}"),
        ("cell-channel law (the axis a21/a22 moved)", 0.0,
         "Hamilton does not touch the emitter's cell law"),
    ]
    print(f"\n=== 4. 208's HAMILTON RIDER, PRICED UNDER EVERY CANDIDATE LAW ===")
    print(f"    rule-13 composite on the shipped renderer = {COMPOSITE_219:+.6f} "
          f"(219 s2.2); floor {FLOOR:+.4f}\n")
    print(f"{'toll law':<44}{'toll scaled':>13}{'on overall':>12}{'composite':>12}"
          f"{'vs floor':>10}{'verdict':>9}")
    out = []
    for name, toll_scaled, why in laws:
        ov = toll_scaled / 6.0
        comp = COMPOSITE_219 + ov
        v = "SHIP" if comp > FLOOR else "VETO"
        out.append(dict(law=name, toll_scaled=toll_scaled, toll_overall=ov,
                        composite=comp, verdict=v, why=why))
        print(f"{name:<44}{toll_scaled:>+13.6f}{ov:>+12.6f}{comp:>+12.6f}"
              f"{comp/FLOOR:>10.2f}{v:>9}")
    ci = 1.96 * nzsd / np.sqrt(3)
    print(f"\n  95% interval on the class mean toll: {nzm:+.6f} +- {ci:.6f} scaled "
          f"= {(nzm-ci)/6:+.6f} .. {(nzm+ci)/6:+.6f} overall")
    print(f"  composite under that interval: {COMPOSITE_219+(nzm-ci)/6:+.6f} .. "
          f"{COMPOSITE_219+(nzm+ci)/6:+.6f}   (floor {FLOOR:+.4f})")
    print(f"\n  For scale: the UNEXPLAINED spread of board d_reach among digest-IDENTICAL "
          f"arms is\n  sd {zsd:.6f} scaled = {zsd/6:.6f} on overall -- "
          f"{zsd/6/FLOOR:.1f}x the whole floor. No\n  toll term of any shape can be estimated "
          f"to the precision this gate needs.")

    json.dump(dict(intended=dict(floor_relL1=fl, ham_relL1=hm_, ratio=fl / hm_),
                   within_class=dict(beta=beta.tolist(), r2=r2, pred_ham=pred_ham_lin),
                   cell_channel=dict(beta=bc.tolist(), r2=r2c,
                                     loo_rmse=float(np.sqrt(np.mean(np.square(loo))))),
                   prices=out, class_ci=float(ci)),
              open(WORK / "phase2_230.json", "w"), indent=1)
    print("\nwrote work/phase2_230.json")


if __name__ == "__main__":
    main()
