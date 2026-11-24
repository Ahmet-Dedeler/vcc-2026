"""804p — the same known-answer grading as grade804.py, on 801's ABCp folds (A/B/C with day-recipe held-out operators).

Reads experiments/801-lab-folds/work/ka-abcp/bulk/ABCp__<held>/<arm>.npz (ka804.py, unmodified, run by 801 on EC2), prices
each A/B/C row with stats800 at the frozen settings, and writes work/grade804p.json beside 062's grade804.json for comparison.
grade804.py itself is imported read-only (its ROWS, board pairs and grade rule); its main() is not called.
    uv run python experiments/804-known-answers/grade804p.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grade804 as G  # noqa: E402

RUNP = HERE.parents[1] / "experiments/801-lab-folds/work/ka-abcp"
RUN710P = HERE / "work" / "run710p"
G.RUN = RUNP
G.surface_of = lambda name: "ABC" if name.startswith("ABCp__") else "W"
S = G.S


def main():
    old = json.loads((HERE / "work" / "grade804.json").read_text())
    out = {}
    for row, (pk, cls, sides, hf) in G.ROWS.items():
        if "ABC" not in sides:
            continue
        cand, base = sides["ABC"]
        # row 710's arrays come from price710p.py (run ka804-710p-b on 801's ABCp designs), the rest from 801's ka-abcp
        G.RUN = RUN710P if row == "710" else RUNP
        F = G.load_bulk("ABC", cand, base, hf)
        if not F:
            print(f"{row:7s} ABCp: arms absent ({cand} / {base})")
            continue
        brd = G.board_raw(pk)
        res = S.price(F, members=("pds", "mse"), n_boot=4000, seed=800)
        dec = S.decide(res["ci95"][0], res["point"], res["ci95"][1])
        m = res["members"]
        rec = dict(row=row, surface="ABCp", cls=cls, cand=cand, base=base, n_folds=len(F), point=res["point"],
                   ci95=res["ci95"], p_le0=res["p_le0"], folds_pos=res["folds_pos"], decision=dec,
                   pds_raw=m["pds"]["raw"], mse_raw=m["mse"]["raw"], r=m["r"]["raw"], board=brd,
                   grade=G.grade(brd["overall"] if brd else None, dec),
                   per_fold=[dict(fold=p["fold"], overall=p["overall"], r=p["raw"]["r"], pds=p["raw"]["pds"])
                             for p in res["per_fold"]])
        if brd and brd["overall"]:
            rec["rho_overall"] = res["point"] / brd["overall"]
        o = old.get(f"{row}|ABC")
        rec["abc062"] = dict(point=o["point"], ci95=o["ci95"], decision=o["decision"], grade=o["grade"]) if o else None
        out[row] = rec
        b = f"{brd['overall']:+.5f}" if brd else "  (prosp.)"
        o_s = f"{o['point']:+.5f} {o['decision']:5s} {o['grade']:11s}" if o else "-"
        rho = f"{rec['rho_overall']:.2f}" if rec.get("rho_overall") is not None else "-"
        print(f"{row:7s} ABCp n={len(F)} {res['point']:+.5f} [{res['ci95'][0]:+.5f}, {res['ci95'][1]:+.5f}] {dec:5s} "
              f"-> {rec['grade']:11s} rho {rho:5s} | 062: {o_s} | board {b} | r {100 * m['r']['raw']:+.2f}% "
              f"pds {m['pds']['raw']:+.5f} folds+ {res['folds_pos']}/{len(F)}", flush=True)
    (HERE / "work" / "grade804p.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
