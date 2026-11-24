#!/usr/bin/env python3
"""528 — follow-up to fold528: the SCALE of excess dropout differs by line (median x on expressed genes: K562 0.005, H1
0.009, HEK293T 0.012, HCT116 0.016, A/B/C 0.015-0.020), and fold528's capture-optimal beta moves inversely with it. Test two
scale-free forms of the key on every fold (and H1 as the reference), both with H1's beta 22:
  norm : x * (median_H1(x) / median_line(x))     (a line-level rescale to H1's median)
  qm   : x quantile-mapped onto H1's x distribution over the 5-CPM universe
against the raw form and the shuffled falsifier. Uses fold528's cached features (work/feat528_<fold>.npz) and 480's gates.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/fold528b.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from abc528 import qmap  # noqa: E402
from fold528 import GATES, block  # noqa: E402

BETA = 22.0


def main():
    H = dict(np.load(HERE / "work/feat528_H1.npz", allow_pickle=True))
    h1x = H["x"][H["keep"]]
    med_h1 = float(np.median(h1x))
    out = {"median_x_H1": med_h1}
    for fold in ("orion_hct116", "orion_hek293t", "k562_gwps"):
        F = dict(np.load(HERE / f"work/feat528_{fold}.npz", allow_pickle=True))
        z = np.load(GATES / f"gate-{fold}.npz", allow_pickle=True)
        gg = z["genes"].astype(str)
        pos = {g: i for i, g in enumerate(F["genes"].astype(str))}
        col = np.array([pos[g] for g in gg])
        x = F["x"][col]
        lc = np.log10(np.maximum(F["mean_cpm"][col], 1e-3))
        med = float(np.median(F["x"][F["keep"]]))
        xn = x * (med_h1 / med)
        xq = qmap(x, h1x)
        targets = z["targets"].astype(str)
        q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
        gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gg)
        excl = np.isin(gg, targets)
        b0 = block(lc, excl)
        rows = {}
        for nm, key in (("cpm", lc), ("raw22", lc + BETA * x), ("norm22", lc + BETA * xn), ("qm22", lc + BETA * xq)):
            b = block(key, excl)
            rows[nm] = dict(capture=float(gate[:, b].sum() / gate.sum()), swapped=int((b & ~b0).sum()))
        # beta grid on the normalised form: where is the optimum now?
        grid = [5.0, 10.0, 15.0, 22.0, 30.0, 40.0]
        rows["norm_grid"] = {f"{b_:g}": float(gate[:, block(lc + b_ * xn, excl)].sum() / gate.sum()) for b_ in grid}
        out[fold] = dict(median_x=med, rows=rows)
        c0 = rows["cpm"]["capture"]
        print(f"{fold}: median x {med:.4f} (H1 {med_h1:.4f}); " + "  ".join(
            f"{k} {v['capture']:.4f} ({v['capture'] - c0:+.4f}, sw {v['swapped']})" for k, v in rows.items() if k != "norm_grid")
            + "  | norm grid: " + " ".join(f"b{k} {v - c0:+.4f}" for k, v in rows["norm_grid"].items()), flush=True)
    (HERE / "work/fold528b.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
