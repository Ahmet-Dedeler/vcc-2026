"""313 -- how big is D5 (the panel columns inside production's norm) at PRODUCTION's panel share?

The mandated fixture is 20 targets on a 60-gene axis, so its panel is 33 % of the axis; production's panel is
300 of 18,533 = 1.6 %.  D5 is the only divergence whose SIZE is a function of that ratio, so it is swept here on
the same construction at five axis widths, reporting the per-target weight difference between the two weight
rules with everything else held fixed (no PerturbMatch, no ordering term).  Weights only; nothing else changes.

  PYTHONDONTWRITEBYTECODE=1 uv run python scale313.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import fixture313 as F
from parity313 import normrel_W, normrel_W_production_rule_rows

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
EXP = -1.5


def one(n_genes: int) -> dict:
    """The fixture's construction at a different axis width: 20 targets, panel = columns 0..19, three donors
    with the same coverage FRACTIONS, the same self-knockdown magnitude, the same seed."""
    rng = np.random.default_rng(F.SEED)
    T = F.N_TARGETS
    tg = [f"T{i:03d}" for i in range(T)]
    shared = rng.normal(0.0, 0.35, size=(T, n_genes))
    cover = {
        "k562_gwps": np.arange(0, int(round(n_genes * 50 / 60))),
        "orion_hct116": np.arange(0, n_genes),
        "orion_hek293t": np.arange(int(round(n_genes * 10 / 60)), n_genes),
    }
    rows, avail = [], []
    for name in F.DONORS:
        cv = cover[name]
        R = np.zeros((T, n_genes))
        a = np.zeros(T, bool)
        for ti, t in enumerate(tg):
            if t in F.PERT_MISSING[name]:
                continue
            row = np.zeros(n_genes)
            row[cv] = 0.55 * shared[ti, cv] + rng.normal(0.0, 0.30, size=cv.size)
            if ti in cv:
                row[ti] = -2.4 + rng.normal(0.0, 0.4)      # the donor's own knockdown, on panel column ti
            R[ti] = row
            a[ti] = True
        if name == F.BIG_ROW[0]:
            R[tg.index(F.BIG_ROW[1])] *= 8.0
        if name == F.SMALL_ROW[0]:
            R[tg.index(F.SMALL_ROW[1])] *= 0.05
        rows.append(R)
        avail.append(a)
    avail = np.stack(avail)
    panel = np.arange(T)
    non_panel = np.arange(T, n_genes)

    rows_h = [R.copy() for R in rows]                       # the 308 harness substrate: panel zeroed
    for R in rows_h:
        R[:, panel] = 0.0
    W_308 = normrel_W(rows_h, avail, non_panel, EXP)                       # norm on non-panel columns
    W_prod = normrel_W_production_rule_rows(rows, avail, np.arange(n_genes), EXP)   # norm on the full axis
    live = np.isfinite(W_308).all(1) & np.isfinite(W_prod).all(1)
    d = np.abs(W_308[live] - W_prod[live])
    share = []
    for j in range(len(rows)):
        e_all = np.linalg.norm(rows[j][avail[j]], axis=1) ** 2
        e_pan = np.linalg.norm(rows[j][avail[j]][:, panel], axis=1) ** 2
        share.append(float(np.median(e_pan / np.maximum(e_all, 1e-300))))
    return dict(n_genes=int(n_genes), panel_share_of_columns=float(T / n_genes),
                median_panel_energy_share=float(np.median(share)),
                max_abs_weight_diff=float(d.max()), median_abs_weight_diff=float(np.median(d)),
                p90_abs_weight_diff=float(np.percentile(d, 90)), n_weighted=int(live.sum()))


def main() -> None:
    out = [one(g) for g in (60, 120, 300, 600, 1200, 18533 * 20 // 300)]
    print(f"{'genes':>7s} {'panel%cols':>10s} {'panelE%':>9s} {'max|dW|':>10s} {'p90|dW|':>10s} {'med|dW|':>10s}")
    for r in out:
        print(f"{r['n_genes']:7d} {100 * r['panel_share_of_columns']:10.2f} "
              f"{100 * r['median_panel_energy_share']:9.2f} {r['max_abs_weight_diff']:10.4e} "
              f"{r['p90_abs_weight_diff']:10.4e} {r['median_abs_weight_diff']:10.4e}")
    (WORK / "scale313.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK / 'scale313.json'}")


if __name__ == "__main__":
    main()
