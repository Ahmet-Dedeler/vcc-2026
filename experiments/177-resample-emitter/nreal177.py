"""177 — `fid` and `jac` across the challenge's UNKNOWN `n_real`, which is the whole verdict.

AGENTS.md 11.1 rule 9: `fid`, `jac` and `reach`'s normalisation all carry `n_real`, and no
harness we hold is in the challenge's regime.  For every previous emission arm that did not
matter, because `n_pred` stayed near 4,400 on the board and the yield cap `min(1, n_pred/n_real)`
was saturated at BOTH ends of the inferred range (240 from 153's enrichment meeting-point and
a21; 2,670 from 093, flagged circular).  **This experiment moves `n_pred` by up to 5.4x, so the
cap is no longer saturated and the sign of `d fid` depends on a quantity we cannot measure.**

Hard board-algebra bounds from 153 (independent of any harness): `jac <= n_real/n_pred` gives
`n_real >= 191`; precision <= 1 on v14 gives `n_real <= 4,716`.

This sweeps `n_real` across that whole admissible range and reports where each arm's `fid` and
`jac` change sign, so the verdict is stated as a boundary rather than a point estimate.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import ANCHOR  # noqa: E402

NP_BOARD = 4424.0            # a15's board-side n_pred (169's convention)
GRID = (191.0, 240.0, 400.0, 600.0, 812.0, 1200.0, 2000.0, 2670.0, 4000.0, 4716.0)
BASE = "base_mf"


def main() -> None:
    rec = json.load(open(os.path.join(HERE, "work", "six177.json")))
    b0 = rec[BASE]
    span_f = ANCHOR["fid"][1] - ANCHOR["fid"][0]
    span_j = ANCHOR["jac"][1] - ANCHOR["jac"][0]
    arms = [a for a in rec if a != BASE]

    print(f"board-side n_pred model: a15 = {NP_BOARD:.0f}, each arm scaled by its own "
          f"h1-test n_pred / {BASE}'s ({b0['n_pred']:.0f})\n")
    print("d `fid` in OVERALL board units, vs n_real")
    print(f"{'arm':<15}{'n_pred_bd':>10}" + "".join(f"{int(n):>10}" for n in GRID))
    for a in arms:
        r = rec[a]
        npb = NP_BOARD * r["n_pred"] / b0["n_pred"]
        cells = []
        for NR in GRID:
            f1 = r["precision"] * min(1.0, npb / NR)
            f0 = b0["precision"] * min(1.0, NP_BOARD / NR)
            cells.append((f1 - f0) / span_f / 6)
        print(f"{a:<15}{npb:>10.0f}" + "".join(f"{c:>+10.5f}" for c in cells))

    print("\nd `jac` in OVERALL board units, vs n_real  "
          "(J = |A n B| / (n_pred + n_real - |A n B|), recall carried from h1-test)")
    print(f"{'arm':<15}{'recall':>10}" + "".join(f"{int(n):>10}" for n in GRID))
    for a in arms:
        r = rec[a]
        npb = NP_BOARD * r["n_pred"] / b0["n_pred"]
        cells = []
        for NR in GRID:
            j1 = r["recall"] * NR / (npb + NR - r["recall"] * NR)
            j0 = b0["recall"] * NR / (NP_BOARD + NR - b0["recall"] * NR)
            cells.append((j1 - j0) / span_j / 6)
        print(f"{a:<15}{r['recall']:>10.4f}" + "".join(f"{c:>+10.5f}" for c in cells))

    print("\nSUM of the four DE members in OVERALL board units, vs n_real")
    print("(`reach` twin-averaged where a twin exists, x0.3930 board transfer; `nmae` from the")
    print(" DE table -- it reads the MEAN of per-cell CPM, so the pseudobulk pin does NOT hold it)")
    span_r = ANCHOR["reach"][1] - ANCHOR["reach"][0]
    NMAE_SPAN = 1.0013 - 0.400
    twin = {"rs_all": "rs_all_b", "rs_sched": "rs_sched_b", "hy_kt2000": "hy_kt2000_b",
            "rs_hkey": "rs_hkey_b", BASE: "base_mf_b"}

    def tw(a, key):
        t = twin.get(a)
        return (rec[a][key] + rec[t][key]) / 2 if t and t in rec else rec[a][key]

    base_reach = tw(BASE, "reach_f1.0")
    base_nmae = tw(BASE, "nmae")
    print(f"{'arm':<15}{'d reach':>9}{'d r ovr':>9}{'d nmae ovr':>11}"
          + "".join(f"{int(n):>10}" for n in GRID))
    for a in arms:
        if a.endswith("_b"):
            continue
        r = rec[a]
        dr = tw(a, "reach_f1.0") - base_reach
        dr_ovr = dr * 0.3930 / span_r / 6
        dn_ovr = -(tw(a, "nmae") - base_nmae) / NMAE_SPAN / 6
        npb = NP_BOARD * r["n_pred"] / b0["n_pred"]
        cells = []
        for NR in GRID:
            f1 = r["precision"] * min(1.0, npb / NR)
            f0 = b0["precision"] * min(1.0, NP_BOARD / NR)
            j1 = r["recall"] * NR / (npb + NR - r["recall"] * NR)
            j0 = b0["recall"] * NR / (NP_BOARD + NR - b0["recall"] * NR)
            cells.append(dr_ovr + dn_ovr + (f1 - f0) / span_f / 6 + (j1 - j0) / span_j / 6)
        print(f"{a:<15}{dr:>+9.5f}{dr_ovr:>+9.5f}{dn_ovr:>+11.5f}"
              + "".join(f"{c:>+10.5f}" for c in cells))
    print("\n`pds` = 0.000000 exactly (pseudobulk sha256 identical on 85/85 targets) and")
    print("`mse` = 0.000000 exactly (`claim177.py`: every arm above the 0.23268 knee).")


if __name__ == "__main__":
    main()
