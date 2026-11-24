"""241 — the PLOLO transfer-cosine delta of every component, in PLOLO units.

FALSIFIER 4 is binding here: the PLOLO -> board cosine factor is NOT identified (1.384 vs
1.000; 231 measured a 7.84x spread over four board pairs).  So this script quotes PLOLO
cosine in PLOLO units and NAMES its fold set, functional, gate and tag.  The rate-card fill
it prints for members h1-test cannot measure is a RANGE over k_cos in {1.000, 1.384}, marked
as an unidentified conversion, and it is never the headline for any row.

  functional : `plolo.endpoints`'s `cos` — the macro (per-row, then mean) cosine between the
               DELIVERED log-pseudobulk effect and the held library's own effect, on the
               fold's scored gene set
  gate       : the arm's own `--min-cpm` (10.0 except in the gate rows)
  frames     : challenge control frames A/B/C, mean over the three
  fold set   : 3 cancer folds (k562_gwps, orion_hct116, orion_hek293t) headline;
               `kolf2` (pluripotent) reported separately, never folded in
  tag        : 241, base = MF / NOMF / SHIPPED as labelled
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
# 218's rate card, per +0.01 of ABSOLUTE BOARD cosine, from `cd4f034`'s own board row
RATE = dict(mse=+0.003634, reach=+0.002446, nmae=+0.000980, fid=+0.000847,
            jac=+0.000100, pds=-0.000748)
LOO = ["noWhiten", "noPM", "noCD4", "alpha_v11", "weights_4_1_1", "noGate", "noPrune",
       "noMassFix", "addMassFix"]


def main() -> None:
    rows = json.loads(open(os.path.join(WORK, "pds241.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows})

    def chal(fold, arm, field):
        vs = [R[(fold, arm, f"chal{c}/delivered")][field] for c in "ABC"
              if (fold, arm, f"chal{c}/delivered") in R]
        return float(np.mean(vs)) if len(vs) == 3 else None

    out = {}
    for tag in ("MF", "NOMF", "SHIPPED"):
        base = f"{tag}:BASE"
        b3 = np.mean([chal(f, base, "cos") for f in CANCER])
        bk = chal("kolf2", base, "cos")
        print()
        print("=" * 124)
        print(f"TRANSFER COSINE, base {tag}   (3-cancer base cos {b3:.6f}, kolf2 {bk:.6f})")
        print("=" * 124)
        print(f"  {'component removed':<20}{'Δcos 3-cancer':>15}{'rel %':>9}"
              f"{'Δcos kolf2':>13}{'rel %':>9}"
              f"{'DE-4 fill @k=1.000':>21}{'@k=1.384':>11}")
        rec = {}
        for a in LOO:
            name = f"{tag}:{a}"
            if name not in arms:
                continue
            c3 = [chal(f, name, "cos") for f in CANCER]
            if any(v is None for v in c3):
                d3 = rel3 = float("nan")
            else:
                d3 = float(np.mean(c3)) - b3
                rel3 = 100.0 * d3 / b3
            ck = chal("kolf2", name, "cos")
            dk = (ck - bk) if ck is not None else float("nan")
            relk = 100.0 * dk / bk if ck is not None else float("nan")
            # the four members h1-test cannot always measure, filled from the rate card;
            # LOO sign convention: positive = the component earns its place
            fills = {}
            for k in (1.000, 1.384):
                dcos_board = -d3 * k                    # removing it loses this much cosine
                fills[k] = sum(RATE[m] * dcos_board / 0.01
                               for m in ("reach", "fid", "jac", "nmae"))
            rec[a] = dict(dcos3=d3, rel3=rel3, dcos_kolf2=dk, rel_kolf2=relk,
                          de4_fill_k1000=fills[1.000], de4_fill_k1384=fills[1.384])
            print(f"  {a:<20}{d3:>+15.6f}{rel3:>+9.2f}{dk:>+13.6f}{relk:>+9.2f}"
                  f"{fills[1.000]:>+21.6f}{fills[1.384]:>+11.6f}")
        out[tag] = dict(base_cos3=float(b3), base_cos_kolf2=float(bk), loo=rec)
    json.dump(out, open(os.path.join(WORK, "cos241.json"), "w"), indent=1, default=float)
    print("\nwrote work/cos241.json")
    print("\nNOTE: `DE-4 fill` is 218's cd4f034 rate card applied to a PLOLO cosine delta "
          "through an UNIDENTIFIED\nk_cos.  It is a sanity range, never a headline "
          "(falsifier 4).")


if __name__ == "__main__":
    sys.exit(main())
