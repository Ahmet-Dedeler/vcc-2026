"""241 — read `work/pds241.json` into leave-one-out `pds` under each renderer base.

Two predictors are carried side by side and neither is hidden:
  4-fold   162's PRE-COMMITTED unweighted mean over all four donor folds — the one 200
           published, so it is the only apples-to-apples comparison to 200's table.
  3-cancer falsifier 3's: `kolf2` is pluripotent, gains 3-5x on this instrument and has
           already forced one band reissue (222).  It is reported SEPARATELY, never folded in.

Scaling: raw -> scaled by the board `pds` span 0.45208 (`debench.ANCHOR`); scaled -> overall
by rule 13's k = 1.011 and /6.  k = 0.899 (162's all-arm fit) is the sensitivity.

PHASE-0 GATES (abort): 200's four SHIP absolutes, the SHIP/SHIP-dup control, the a15->a23
mass-fix price -0.004014 scaled, and 200's published `pds` LOO column.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS4 = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
SPAN = 0.95228 - 0.50020
K13, K899 = 1.011, 0.899
SIGMA = 0.001494

# 200 RESULTS.md Phase 0 / Phase 1
KA_SHIP_ABS = dict(orion_hct116=0.791823, orion_hek293t=0.777375,
                   k562_gwps=0.789121, kolf2=0.610128)
KA_MASSFIX_SCALED = -0.004014
KA_LOO_OVERALL = {          # 200's `pds` column, overall units, sign = "costs us to remove"
    "noWhiten": +0.025593, "weights_4_1_1": +0.001312, "alpha_v11": +0.000022,
    "noGate": +0.001169, "noCD4": +0.002611, "noMassFix": +0.000676,
    "noPrune": +0.000392, "noPM": +0.001885,
}
LOO = ["noWhiten", "noPM", "noCD4", "alpha_v11", "weights_4_1_1", "noGate", "noPrune",
       "noMassFix", "addMassFix"]


def main() -> None:
    rows = json.loads(open(os.path.join(WORK, "pds241.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows})

    def chal(fold, arm):
        vs = [R[(fold, arm, f"chal{c}/delivered")]["pds"] for c in "ABC"
              if (fold, arm, f"chal{c}/delivered") in R]
        return float(np.mean(vs)) if len(vs) == 3 else None

    def delta(arm, base):
        per = {f: chal(f, arm) - chal(f, base) for f in FOLDS4
               if chal(f, arm) is not None and chal(f, base) is not None}
        if not per:
            return None
        c = [per[f] for f in CANCER if f in per]
        allv = np.asarray([per[f] for f in FOLDS4 if f in per])
        d = dict(per=per, n=len(per),
                 mean4=float(allv.mean()), sd4=float(allv.std()),
                 scaled4=float(allv.mean()) / SPAN,
                 signs="".join("+" if per[f] > 0 else ("-" if per[f] < 0 else "0")
                               for f in FOLDS4 if f in per))
        d["overall4"] = d["scaled4"] * K13 / 6
        d["overall4_k899"] = d["scaled4"] * K899 / 6
        if c:
            d["mean3"] = float(np.mean(c))
            d["sd3"] = float(np.std(c))
            d["scaled3"] = d["mean3"] / SPAN
            d["overall3"] = d["scaled3"] * K13 / 6
        if "kolf2" in per:
            d["kolf2"] = per["kolf2"]
            d["kolf2_overall"] = per["kolf2"] / SPAN * K13 / 6
        return d

    out = {}
    print("=" * 122)
    print("ABSOLUTE `pds`, challenge frame, per fold")
    print("=" * 122)
    print(f"  {'arm':<22}" + "".join(f"{f:>16}" for f in FOLDS4)
          + f"{'4-fold':>10}{'3-cancer':>10}")
    ka_abs_ok = True
    for tag in ("MF", "NOMF", "SHIPPED"):
        for suffix in ("BASE", "BASE-dup"):
            a = f"{tag}:{suffix}"
            if a not in arms:
                continue
            v = {f: chal(f, a) for f in FOLDS4}
            m4 = np.mean([v[f] for f in FOLDS4])
            m3 = np.mean([v[f] for f in CANCER])
            print(f"  {a:<22}" + "".join(f"{v[f]:16.6f}" for f in FOLDS4)
                  + f"{m4:10.6f}{m3:10.6f}")
            out[a] = dict(per_fold=v, mean4=float(m4), mean3=float(m3))
    # Phase-0 gate 1: MF:addMassFix does not exist; MF:BASE IS 200's SHIP
    print()
    print("PHASE-0 GATE 1 — `MF:BASE` must reproduce 200's published SHIP absolutes")
    for f in FOLDS4:
        got, exp = chal(f, "MF:BASE"), KA_SHIP_ABS[f]
        d = abs(got - exp)
        ka_abs_ok &= d < 5e-6
        print(f"  {f:<16} got {got:.6f}  200 {exp:.6f}  |d| {d:.2e}  "
              f"{'PASS' if d < 5e-6 else 'FAIL'}")
    # Phase-0 gate 2: the state-carry control
    for tag in ("MF", "NOMF", "SHIPPED"):
        w = max(abs(chal(f, f"{tag}:BASE") - chal(f, f"{tag}:BASE-dup")) for f in FOLDS4)
        print(f"  control {tag}:BASE vs BASE-dup  worst |d| = {w:.2e}  "
              f"{'PASS' if w == 0.0 else 'FAIL'}")
    # Phase-0 gate 3: the mass-fix price
    d = delta("MF:noMassFix", "MF:BASE")
    print(f"  a15 -> a23 mass fix, scaled  got {d['scaled4']:+.6f}  "
          f"200 {KA_MASSFIX_SCALED:+.6f}  |d| {abs(d['scaled4']-KA_MASSFIX_SCALED):.2e}  "
          f"{'PASS' if abs(d['scaled4']-KA_MASSFIX_SCALED) < 5e-6 else 'FAIL'}")

    def block(tag, title):
        print()
        print("=" * 122)
        print(title)
        print("=" * 122)
        print(f"  {'component removed':<20}{'raw d4':>10}{'sd4':>9}{'scaled4':>10}"
              f"{'ovr4 k1.011':>13}{'ovr3 cancer':>13}{'kolf2 ovr':>11}"
              f"{'200 pds':>10}{'signs':>8}")
        res = {}
        for a in LOO:
            name = f"{tag}:{a}"
            if name not in arms:
                continue
            r = delta(name, f"{tag}:BASE")
            if r is None:
                continue
            res[a] = r
            # LOO sign flip: positive = the component EARNS its place
            v4, v3 = -r["overall4"], -r.get("overall3", float("nan"))
            vk = -r.get("kolf2_overall", float("nan"))
            exp = KA_LOO_OVERALL.get(a)
            print(f"  {a:<20}{r['mean4']:>+10.6f}{r['sd4']:>9.6f}{r['scaled4']:>+10.6f}"
                  f"{v4:>+13.6f}{v3:>+13.6f}{vk:>+11.6f}"
                  f"{(f'{exp:+.6f}' if exp is not None else '—'):>10}{r['signs']:>8}"
                  + ("" if r["n"] == 4 else f"   ({r['n']} fold only)"))
        return res

    out["MF"] = block("MF", "LEAVE ONE OUT — base = 200's renderer (VCC_MASS_FIX=minnorm)")
    out["NOMF"] = block("NOMF", "LEAVE ONE OUT — base = THE SHIPPED RENDERER (mass fix unset)")
    out["SHIPPED"] = block("SHIPPED",
                           "LEAVE ONE OUT — base = shipped renderer AT THE SHIPPED CD4 DOSE f=0.34")

    # Phase-0 gate 4: 200's own pds LOO column
    print()
    print("PHASE-0 GATE 4 — 200's published `pds` LOO column, reproduced on the MF base")
    worst = 0.0
    for a, exp in KA_LOO_OVERALL.items():
        if a not in out["MF"]:
            continue
        got = -out["MF"][a]["overall4"]
        worst = max(worst, abs(got - exp))
        print(f"  {a:<20} got {got:+.6f}  200 {exp:+.6f}  |d| {abs(got-exp):.2e}  "
              f"{'ok' if abs(got-exp) < 1e-5 else 'MISMATCH'}")
    print(f"  worst |d| {worst:.2e}   {'PASS' if worst < 1e-5 else 'FAIL'}")

    json.dump(out, open(os.path.join(WORK, "report_pds241.json"), "w"), indent=1,
              default=float)
    print("\nwrote work/report_pds241.json")


if __name__ == "__main__":
    sys.exit(main())
