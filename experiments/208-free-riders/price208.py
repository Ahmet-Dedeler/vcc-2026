"""208 Phase 3 — assemble the six-member price of each rider and run rule 13's gate.

Nothing is computed here that was not measured by `pds208.py`, `mse208.py` or `six208.py`;
this file only reads their JSON and applies rule 13, so the composite in RESULTS.md is
reproducible from one command and cannot drift from the tables above it.

    uv run python price208.py
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

CALIB = 1.011                    # 162's fitted PLOLO -> board scale for `pds`
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
# The three readings of rule 13's `reach` toll, in OVERALL units.  AGENTS.md is internally
# inconsistent about which applies to an arm that moves the delivered per-gene TOTALS but
# neither the emission law nor the effect, so all three are carried.
TOLLS = {
    "none (196's reading; 169 measured +0.0006 for profile arms)": 0.0,
    "a23's REALISED profile-arm reach cost (-0.007933 scaled)": -0.007933 / 6,
    "the emission-law toll (-0.0146 scaled), charged anyway": -0.0146 / 6,
}
ARMS = {"ham_tier": "A  — Hamilton rounding",
        "floor_flat": "B  — flat 1.00 spread",
        "ham_flat": "A o B — composed"}
# Which of pds / mse each arm can move at all.  The spread policy cannot move the delivered
# per-gene totals (`emit_counts` computes them before `spread` is read) and `emit208.py`
# asserted the digest identity on 85/85 h1-test blocks, so B's `pds` is 0 by identity.
PDS_MOVES = {"ham_tier": True, "floor_flat": False, "ham_flat": True}


def main() -> None:
    pds = {r["rule"]: r for r in json.loads((WORK / "pds208.json").read_text())["rows"]}
    mse = json.loads((WORK / "mse208.json").read_text())["rows"]
    six = json.loads((WORK / "six208.json").read_text())["board"]
    base_mse = mse["rc_callset|floor|tier"]

    d_pds_ham = pds["ham"]["pds_overall"] * CALIB
    print(f"`pds`  PLOLO 4-fold, floor {pds['floor']['pds']:.6f} -> ham "
          f"{pds['ham']['pds']:.6f}  = {pds['ham']['d_pds']:+.6f} raw = "
          f"{pds['ham']['pds_overall']:+.6f} overall, x{CALIB} = {d_pds_ham:+.6f}")
    print(f"       macro cosine {pds['floor']['cos']:.6f} -> {pds['ham']['cos']:.6f}  "
          f"({100*pds['ham']['cos_rel']:+.2f}% relative, "
          f"{pds['ham']['cos_overall']:+.6f} overall) — INSIDE 204's +-0.20% resolution")
    print(f"`mse`  claim {base_mse['claim']:.5f} -> ", end="")
    for k in ("rc_callset|ham|tier", "rc_callset|floor|flat", "rc_callset|ham|flat"):
        print(f"{mse[k]['claim']:.5f} ", end="")
    print(f" — every arm SATURATED, d overall 0.00000")

    for pool in ("native", "wald89", "matched"):
        print(f"\n=== SIX-MEMBER PRICE, marginal to `rc_callset`  (DE pool: {pool}) ===")
        print(f"{'member':<10}{'instrument':<34}" +
              "".join(f"{ARMS[a].split(' ')[0]:>12}" for a in ARMS))
        rows = {}
        for a in ARMS:
            b = six[f"{pool}|{a}"]
            m = mse[f"rc_callset|{'ham' if a.startswith('ham') else 'floor'}|"
                    f"{'flat' if a.endswith('flat') else 'tier'}"]
            rows[a] = dict(
                pds=d_pds_ham if PDS_MOVES[a] else 0.0,
                mse=m["d_overall"] - base_mse["d_overall"],
                reach=b["reach"], fid=b["fid"], jac=b["jac"], nmae=b["nmae"],
                reach_alt=b["reach_alt"])
            rows[a]["sum"] = sum(rows[a][k] for k in
                                 ("pds", "mse", "reach", "fid", "jac", "nmae"))
            rows[a]["sum_alt"] = rows[a]["sum"] - b["reach"] + b["reach_alt"]
        for member, inst in (("pds", "PLOLO 4-fold x 1.011"),
                             ("mse", "096 `claim` route, challenge"),
                             ("reach", "h1-test x 1.290 (rc_hybrid's)"),
                             ("fid", "h1-test precision"),
                             ("jac", "Phase 1's S* convention, NP pinned"),
                             ("nmae", "h1-test raw / span")):
            print(f"{member:<10}{inst:<34}" +
                  "".join(f"{rows[a][member]:>+12.6f}" for a in ARMS))
        print(f"{'SUM':<10}{'':<34}" + "".join(f"{rows[a]['sum']:>+12.6f}" for a in ARMS))
        print(f"{'SUM@0.4277':<10}{'(reach sensitivity)':<34}" +
              "".join(f"{rows[a]['sum_alt']:>+12.6f}" for a in ARMS))
        add = rows["ham_tier"]["sum"] + rows["floor_flat"]["sum"]
        print(f"\n  COMPOSITION: A {rows['ham_tier']['sum']:+.6f} + B "
              f"{rows['floor_flat']['sum']:+.6f} = {add:+.6f}   MEASURED A o B "
              f"{rows['ham_flat']['sum']:+.6f}   super-additivity "
              f"{rows['ham_flat']['sum']-add:+.6f}")
        if pool == "native":
            keep = rows

    # ---- rule 13 --------------------------------------------------------------------
    print("\n=== RULE 13, THE COMPOSED SLOT GATE ===")
    print("  rule 13 prices ONLY `pds` (PLOLO) and `mse` (096 claim).  The +0.0019 floor IS")
    print("  the realised RMS of the four DE members it does not price, so adding them to")
    print("  the left-hand side and comparing to the same floor would double-count.  Both")
    print("  numbers are given.\n")
    print(f"{'arm':<22}{'priced (pds+mse)':>18}{'toll':>12}{'composite':>12}"
          f"{'gate':>9}{'verdict':>9}{'  +DE (best estimate)':>22}")
    out = {}
    for a in ARMS:
        priced = keep[a]["pds"] + keep[a]["mse"]
        for tname, toll in TOLLS.items():
            comp = priced + toll
            v = "SHIP" if comp > GATE_SHAPE else "VETO"
            if tname.startswith("none"):
                print(f"{ARMS[a]:<22}{priced:>18.6f}{toll:>12.6f}{comp:>12.6f}"
                      f"{GATE_SHAPE:>9.4f}{v:>9}{keep[a]['sum']:>22.6f}")
            else:
                print(f"{'  '+tname[:18]:<22}{'':>18}{toll:>12.6f}{comp:>12.6f}"
                      f"{GATE_SHAPE:>9.4f}{v:>9}")
            out.setdefault(a, {})[tname] = dict(priced=priced, toll=toll, composite=comp,
                                                verdict=v)
        out[a]["six_member_sum"] = keep[a]["sum"]
    json.dump(dict(rows={a: keep[a] for a in ARMS}, gate=out, calib=CALIB,
                   gate_shape=GATE_SHAPE), open(WORK / "price208.json", "w"), indent=1,
              default=float)


if __name__ == "__main__":
    main()
