"""Same arms, but every SHAPE arm re-solved to a12's own amplitude (`A` held exactly), which is
what `--effect-scale` does in the builder.  Shape and amplitude are then separate columns:
`lfc-scale` is the amplitude question, everything else is the shape question at constant A.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae")
from arm_table import Bench, arms, DEPTHS
import nmaekit as K
import numpy as np

for W_ in ("h1", "orion"):
    Bn = Bench(W_); b = Bn.base; H = Bn.H
    print(f"\n===== {W_}  a12 base: nmae {b['n0.05']:.4f}/{b['n0.35']:.4f}/{b['n0.6']:.4f}"
          f"  pds {b['pds']:.4f}  A {b['A']:.5f}  X {b['X']:.6f} =====")
    print(f"  {'arm (A held at a12)':<26}{'X/X0':>7}{'B/B0':>7}{'pds':>8}{'Δpds':>8}"
          + "".join(f"{'Δnmae@'+str(d):>12}" for d in DEPTHS))
    for nm, kw in arms(Bn).items():
        if nm.startswith("de-offset") or nm.startswith("lfc shift"): continue
        r = Bn.metrics(kw["gain"], kw.get("shift") or 0.0, hold_A=True)
        print(f"  {nm:<26}{r['X']/b['X']:>7.3f}{r['B']/b['B']:>7.3f}{r['pds']:>8.4f}"
              f"{r['pds']-b['pds']:>+8.4f}"
              + "".join(f"{r['n'+str(d)]-b['n'+str(d)]:>+12.4f}" for d in DEPTHS))
