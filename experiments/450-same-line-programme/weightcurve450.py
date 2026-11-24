"""450g (2026-09-26, after 446bja's board read): how much of a same-line programme's gain is already in a consensus that
contains the same-line bank at weight ~1/4? H1 analogue of context A, where CLEAN4 = K562-GWPS, K562-ess, Nadig Jurkat,
Nadig HepG2 (Jurkat at ~1/4) and 446bja moved A's prior to Nadig alone (weight 1).
Predictor: sign of (1-w) z(P_down CLEAN4) + w z(P_down Pan), w in [0, 1]; score gate-block sign accuracy on h1-test
(sameline450's truth, block and gate). Light: loads what sameline450 loads."""
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
c4p, pp, gate, block, truth, adj = (g[k] for k in ("c4p", "pp", "gate", "block", "truth", "adj"))
z = g["z"]
zc, zp = z(c4p), z(pp)
both = np.isfinite(zc) & np.isfinite(zp)
print(f"genes with both priors: {both.sum()} of {len(zc)}")
rows = []
for w in (0.0, 0.1, 0.2, 0.25, 0.33, 0.5, 0.67, 0.8, 1.0):
    s = np.where(both, -np.sign((1 - w) * zc + w * zp), np.nan)   # P_down high -> predicted down (-1)
    for name, m in (("gate_block", gate & block), ("gate_all", gate), ("block_all_adj", adj & block)):
        mm = m & np.isfinite(s)[None, :]
        acc = float(np.mean((np.broadcast_to(s, truth.shape) == truth)[mm]))
        rows.append((w, name, acc, int(mm.sum())))
for name in ("gate_block", "gate_all", "block_all_adj"):
    r = [x for x in rows if x[1] == name]
    base, top = r[0][2], r[-1][2]
    print(f"{name:14s} n {r[0][3]:6d}: " + "  ".join(f"w{x[0]:g} {x[2]:.4f}" for x in r)
          + f"   | share of the w0->w1 gain already at w 0.25: {(r[3][2] - base) / (top - base + 1e-12):.2f}")
