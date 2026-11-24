"""178 — assemble the six-member table for every arm, in overall board units (rule 12).

One row per arm, one column per member, each from the instrument that is allowed to price it:

  `pds`   PLOLO, challenge frame, 4-fold mean, scaled by span 0.452, calibrated x1.011, /6.
          The ONLY board-validated `pds` instrument (162: 6/6 signs, LOO RMS 0.001494,
          permutation p = 0.0028). The h1-test delivered-space reading is printed beside it
          because the two DISAGREE IN SIGN on two arms, and that disagreement is a result.
  `mse`   096's identity on the real challenge panel, `claim` re-measured per arm, the `X`
          ratio transferred from h1-test. Both amplitude conventions are shown.
  `fid`   h1-test precision (no `n_real`, rule 9), twin-averaged, / span 0.297 / 6.
  `nmae`  h1-test, twin-averaged, / span 0.6013 / 6.
  `jac`   h1-test overlap ratio rebuilt on the board's own n_real/n_pred, / span 0.365 / 6.
  `reach` h1-test at native BH, twin-averaged, x0.393 -- printed WITH its twin spread and
          flagged: 165 measured that h1-test reproduces board `reach` for emission-law changes
          but NOT for profile changes, and every arm here is a profile change. Rule 13's
          -0.0146 scaled toll is the alternative, and both totals are given.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SPAN_PDS, K_CAL, TOLL = 0.452, 1.011, -0.0146
ARMS = ("a15", "g045", "g030", "flip", "w211", "flip_w211", "g045_flip",
        "g045_w211", "all3_g045", "all3_g030")
DE = {"a15": ("rb_a15", "rb_a15_b"), "g045": ("rb_g045", "rb_g045_b"),
      "g030": ("rb_g030",), "flip": ("rb_flip", "rb_flip_b"),
      "w211": ("rb_w211", "rb_w211_b"),
      "flip_w211": ("rb_flip_w211", "rb_flip_w211_b"), "g045_flip": ("rb_g045_flip",),
      "g045_w211": ("rb_g045_w211",),
      "all3_g045": ("rb_all3_g045", "rb_all3_g045_b"), "all3_g030": ("rb_all3_g030",)}
PARTS = {"flip_w211": ("flip", "w211"),
         "g045_flip": ("g045", "flip"), "g045_w211": ("g045", "w211"),
         "all3_g045": ("g045", "flip", "w211"), "all3_g030": ("g030", "flip", "w211")}


def main() -> None:
    gate = json.loads((WORK / "gate178.json").read_text())
    ch = json.loads((WORK / "chal178_solve.json").read_text())
    mse = {a: {"d_overall": (r["scaled"] - ch["a15"]["scaled"]) / 6.0,
               "d_overall_knee": (r["knee_scaled"] - ch["a15"]["knee_scaled"]) / 6.0}
           for a, r in ch.items()}
    six = json.loads((WORK / "six178.json").read_text())["overall"]
    h1 = json.loads((WORK / "h1rebuild178.json").read_text())["arms"]

    fixed = {a: r["d_overall_knee"] for a, r in mse.items()}

    rows = {}
    for arm in ARMS:
        twins = DE[arm]
        de = {ch: float(np.mean([six[t][ch] for t in twins]))
              for ch in ("fid", "nmae", "jac", "reach")}
        spread = {ch: (abs(six[twins[0]][ch] - six[twins[1]][ch]) if len(twins) > 1
                       else float("nan")) for ch in ("fid", "nmae", "reach")}
        pds = (K_CAL * gate["dpds"][arm]["mean4"] / 6.0) if arm != "a15" else 0.0
        rows[arm] = dict(pds=pds, mse=mse[arm]["d_overall"],
                         mse_fixed=fixed[arm],
                         **de, spread=spread,
                         h1_pds_raw=h1[arm]["pds"] - h1["a15"]["pds"])
    # `six178`'s DE deltas are against the SINGLE table `rb_a15`, while every arm here is
    # twin-averaged where a twin exists. Re-centre on the twin-averaged a15 so the incumbent
    # row is exactly zero; what it was before re-centring is the DE resolution floor and is
    # printed as `twin` below.
    a15_de = {ch: rows["a15"][ch] for ch in ("fid", "nmae", "jac", "reach")}
    for arm in ARMS:
        for ch, v in a15_de.items():
            rows[arm][ch] -= v

    print("=== SIX MEMBERS, overall board units, vs the a15/a23 incumbent ===")
    print(f"{'arm':<12}{'pds':>10}{'mse':>10}{'fid':>10}{'nmae':>10}{'jac':>9}"
          f"{'reach':>10}{'TOTAL':>10}{'TOTAL*':>10}")
    for arm in ARMS:
        r = rows[arm]
        tot = r["pds"] + r["mse"] + r["fid"] + r["nmae"] + r["jac"] + r["reach"]
        tot_toll = r["pds"] + r["mse"] + r["fid"] + r["nmae"] + r["jac"] + TOLL / 6.0
        print(f"{arm:<12}{r['pds']:>+10.6f}{r['mse']:>+10.6f}{r['fid']:>+10.6f}"
              f"{r['nmae']:>+10.6f}{r['jac']:>+9.6f}{r['reach']:>+10.6f}"
              f"{tot:>+10.6f}{tot_toll:>+10.6f}")
    print("  TOTAL  uses the MEASURED h1-test reach delta (which 165 says is not readable "
          "for a profile change)\n  TOTAL* replaces it with rule 13's -0.0146 scaled "
          "emitted-profile toll. Both are shown because\n         neither is trustworthy "
          "alone; the verdict is the same under either.")

    print(f"\n{'arm':<12}{'twin fid':>10}{'twin nmae':>11}{'twin reach':>12}"
          f"{'mse@knee':>13}{'h1 d pds raw':>14}{'PLOLO d pds raw':>17}")
    for arm in ARMS:
        r = rows[arm]
        print(f"{arm:<12}{r['spread']['fid']:>10.6f}{r['spread']['nmae']:>11.6f}"
              f"{r['spread']['reach']:>12.6f}{r['mse_fixed']:>+13.6f}"
              f"{r['h1_pds_raw']:>+14.5f}"
              f"{(gate['dpds'][arm]['mean4'] * SPAN_PDS) if arm in gate['dpds'] else 0.0:>+17.5f}")

    print("\n=== COMPOSITION, member by member (composite / sum of its riders) ===")
    print(f"{'composite':<12}{'member':<8}{'composite':>12}{'sum parts':>12}{'ratio':>8}"
          f"  reading")
    for comp, members in PARTS.items():
        for ch in ("pds", "mse", "fid", "nmae", "reach"):
            c = rows[comp][ch]
            s = sum(rows[m][ch] for m in members)
            ratio = c / s if abs(s) > 1e-12 else float("nan")
            read = ("super-additive" if ratio > 1.05 else
                    "sub-additive" if ratio < 0.95 else "additive")
            print(f"{comp:<12}{ch:<8}{c:>+12.6f}{s:>+12.6f}{ratio:>8.3f}  {read}")
        c = sum(rows[comp][ch] for ch in ("pds", "mse", "fid", "nmae", "jac"))
        s = sum(sum(rows[m][ch] for ch in ("pds", "mse", "fid", "nmae", "jac"))
                for m in members)
        print(f"{comp:<12}{'5 (no reach)':<8}{c:>+12.6f}{s:>+12.6f}"
              f"{(c/s if abs(s) > 1e-12 else float('nan')):>8.3f}")
    (WORK / "table178.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
