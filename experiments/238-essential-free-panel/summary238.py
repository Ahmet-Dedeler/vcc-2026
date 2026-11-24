"""238 — the head-to-head tables RESULTS.md quotes, joined from the five artefact files.

Reads only `work/*.json`; computes nothing new.  Every number it prints is traceable to the
script that measured it:

  census238.json   the 0/34/30 row, its list provenance, per-source panel coverage
  carry238.json    per-fold carrying rows and their ESSENTIAL fraction, all five panels
  six238.json      all six members, every arm, every (panel, fold)
  price238.json    228's board-unit routes, raw and null-corrected
  retain238.json   233's retention ratios and board-unit projection

    uv run python summary238.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels238 as P  # noqa: E402

W = P.WORK
PANELS = ("abc", "real", "pess", "efr", "efp")
LABEL = {"abc": "abc (Arc's real 300)", "real": "real (199 realistic draw)",
         "pess": "pess (199 pessimistic draw)", "efr": "efr (realistic, ESSENTIAL-FREE)",
         "efp": "efp (pessimistic, ESSENTIAL-FREE)"}
OCTOBER = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
ALL4 = "nadig_hepg2+nadig_jurkat+rpe1+k562_essential"
ALL3 = "nadig_hepg2+nadig_jurkat+rpe1"


def load(name):
    p = W / name
    return json.loads(p.read_text()) if p.exists() else None


def sec(t):
    print(f"\n{'=' * 104}\n{t}\n{'=' * 104}")


def main() -> None:
    cen, carry = load("census238.json"), load("carry238.json")
    price, ret = load("price238.json"), load("retain238.json")
    six = load("six238.json") or []

    # ------------------------------------------------------------------ 1. panel properties
    sec("1.  PANEL PROPERTIES — the two Arc's real panel has, and 199's draws do not both have")
    print(f"{'panel':30}{'n':>5}{'∩ abc':>8}{'∩ essential':>13}{'k562':>7}{'hct116':>8}"
          f"{'hek293t':>9}{'CD4 union':>11}")
    cd4 = load("cd4_238.json") or {}
    cd4_233 = {}
    p233 = P.R233 / "work/cd4_233.json"
    if p233.exists():
        cd4_233 = json.loads(p233.read_text())
    for k in PANELS:
        if not P.BUNDLES[k].exists():
            continue
        e = cen["ka1_panel_x_essential"][k]
        cv = cen["coverage_and_essential_fraction"]
        u = (cd4.get(f"{k}|coverage") or cd4_233.get(f"{k}|coverage") or {}).get("union", "—")
        print(f"{LABEL[k]:30}{e['n_panel']:>5}"
              f"{(0 if k == 'abc' else cen['disjoint'].get(k, '—')):>8}"
              f"{e['n_essential']:>13}"
              f"{cv[f'{k}|k562_gwps']['carried']:>7}"
              f"{cv[f'{k}|orion_hct116']['carried']:>8}"
              f"{cv[f'{k}|orion_hek293t']['carried']:>9}{str(u):>11}")
    print(f"\nessential list: {cen['essential_list']['path']}  n={cen['essential_list']['n']}  "
          f"sha256 {cen['essential_list']['sha256'][:16]}…")
    print(f"  upstream: {cen['essential_list']['upstream']}")

    # ------------------------------------------------- 2. carrying rows + essential fraction
    sec("2.  THE OCTOBER FOUR — carrying rows per fold, and how many of them are ESSENTIAL")
    print(f"{'panel':6}{'fold':16}" + "".join(f"{s[:14]:>16}" for s in OCTOBER)
          + f"{'UNION':>14}")
    for k in PANELS:
        for held in P.CANCER_FOLDS:
            r = (carry or {}).get(f"{k}|{held}")
            if not r:
                continue
            cells = []
            for s in OCTOBER:
                v = r["sources"].get(s)
                if v is None:
                    cells.append(f"{'(leak)':>16}")
                elif v["carrying"] == 0:
                    cells.append(f"{'0  VACUOUS':>16}")
                else:
                    cells.append(f"{v['carrying']:>8} ({100 * v['frac']:3.0f}% ess)")
            u = r["sources"]["UNION"]
            print(f"{k:6}{held:16}" + "".join(cells)
                  + f"{u['carrying']:>6} ({100 * u['frac']:3.0f}%)"
                  if u["carrying"] else f"{k:6}{held:16}" + "".join(cells) + f"{0:>14}")

    # ------------------------------------------------------------------ 3. the re-priced arms
    if price:
        sec("3.  THE FOUR SOURCES RE-PRICED — null-corrected (arm minus its own derangement)")
        print(f"{'panel':6}{'arm':<48}{'carry':>7}{'rel cos%':>10}{'d_pds':>10}"
              f"{'d_nmae':>10}{'composite (no reach)':>26}")
        for key in sorted(price):
            if not key.endswith("|NULLCORR"):
                continue
            panel, arm, lab, _ = key.split("|")
            if lab != "CLEAN2":
                continue
            v = price[key]
            c = v["composite_no_reach"]
            print(f"{panel:6}{arm:<48}{str(v['carrying_rows']):>7}{v['rel_cos']:+10.2f}"
                  f"{v['d_pds']:+10.5f}{v['d_nmae']:+10.5f}"
                  f"{c['1.0']:+13.5f} .. {c['1.384']:+.5f}")

        sec("4.  WITHIN-PANEL DECOMPOSITION — the same residual on its ESSENTIAL rows vs its "
            "NON-ESSENTIAL rows")
        print(f"{'panel':6}{'arm':<48}{'rows':>6}{'rel cos%':>10}{'d_pds':>10}"
              f"{'composite (no reach)':>26}")
        for key in sorted(price):
            if not key.endswith("|NULLCORR"):
                continue
            panel, arm, lab, _ = key.split("|")
            if lab != "CLEAN2" or ("_essonly@" not in arm and "_noness@" not in arm
                                   and not arm.startswith(ALL4)
                                   and not arm.startswith(ALL3)):
                continue
            v = price[key]
            c = v["composite_no_reach"]
            print(f"{panel:6}{arm:<48}{str(v['carrying_rows']):>6}{v['rel_cos']:+10.2f}"
                  f"{v['d_pds']:+10.5f}{c['1.0']:+13.5f} .. {c['1.384']:+.5f}")

    # ------------------------------------------------------------------ 5. retention
    if ret:
        sec("5.  RETENTION — credit above the deranged null, new panel / abc")
        keys = ("cos_agg", "pds_at140", "mse", "nmae", "fid", "reach", "jac")
        print(f"{'member':12}" + "".join(f"{p:>12}" for p in PANELS[1:]))
        for k in keys:
            row = []
            for p in PANELS[1:]:
                t = ret.get(f"retention|{p}")
                row.append(f"{100 * t[k][-1]:11.1f}%" if t else f"{'—':>12}")
            print(f"{k:12}" + "".join(row))
        print(f"{'DIRECTION':12}" + "".join(
            f"{100 * ret[f'channel|{p}']['direction']:11.1f}%"
            if f"channel|{p}" in ret else f"{'—':>12}" for p in PANELS[1:]))
        print(f"{'ORDERING':12}" + "".join(
            f"{100 * ret[f'channel|{p}']['ordering']:11.1f}%"
            if f"channel|{p}" in ret else f"{'—':>12}" for p in PANELS[1:]))
        bp = ret.get("board_projection", {})
        print(f"\n{'board projection (233 §3.3)':30}{'overall':>10}{'of today':>11}")
        print(f"{'today (h1union, board)':30}{0.21183096:10.5f}{100.0:10.1f}%")
        for p in PANELS[1:]:
            if p in bp:
                print(f"{LABEL[p]:30}{bp[p]['overall']:10.5f}"
                      f"{100 * bp[p]['of_today']:10.1f}%")

    # ------------------------------------------------------------------ 6. vacuity ledger
    vac = [r for r in six if r.get("vacuous")]
    if vac:
        sec("6.  VACUITY LEDGER — sources the instrument REFUSED to price (0 carrying rows)")
        seen = set()
        for r in vac:
            k = (r["panel"], r["arm"])
            if k in seen:
                continue
            seen.add(k)
            n = sum(1 for x in vac if x["panel"] == r["panel"] and x["arm"] == r["arm"])
            print(f"  {r['panel']:6} {r['arm']:<28} on {n} fold(s); "
                  f"|panel ∩ source| = {r['panel_overlap']}")


if __name__ == "__main__":
    main()
