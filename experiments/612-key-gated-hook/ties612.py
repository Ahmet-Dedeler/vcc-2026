#!/usr/bin/env python3
"""612 (607 item 3, 650's 10-06 procedure): how many priced H1 targets have their design key TIED at the KEYGATE cut?
For each of the 85 h1-test design keys (lib521's S["key"], the gated axis run612 priced on) at q 0.5, it counts:
  ties   = targets with > 1 nonzero |key| exactly equal to np.quantile(|key[nz]|, q);
  filled = targets where patch612's tie fill changes the weak set (need = floor(q·n) − #(|key| < cut) > 0 with ties).
    uv run python experiments/612-key-gated-hook/ties612.py"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
S = L.load_state(0)
K = np.abs(np.asarray(S["key"], float))
q = 0.5
ties = filled = 0
for i in range(K.shape[0]):
    k = K[i]
    nz = k > 0
    if not nz.any():
        continue
    cut = float(np.quantile(k[nz], q))
    t = int((nz & (k == cut)).sum())
    need = int(np.floor(q * int(nz.sum()) + 1e-9)) - int((nz & (k < cut)).sum())
    ties += t > 1
    filled += (t > 0 and need > 0)
print(f"85 H1 targets priced: {ties} with a tie at the q {q} cut; {filled} where the tie fill changes the weak set")
