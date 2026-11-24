#!/usr/bin/env python3
"""613d3: DESCRIPTIVE depth sensitivity of 613d (it does not change 613d's registered 496a prediction). 613d computes each
fold's add from its OWN depth (Orion 12-14k, K562 11k median library), so genes under 0.3 counts/cell get no add. A/B/C
run at ~20k (context A's median library is 20,109). This re-computes the cancer-fold rows with the add at lib 20,109 on
every surface. Output: work/block613d3.out."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613 as B  # noqa: E402
import block613c as BC  # noqa: E402
import block613d as BD  # noqa: E402
import decomp639 as DC  # noqa: E402
import ds612b as DB  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
LIB = 20109.0
I = Inst()
prof = DC.flex_profiles()
lines = []
vals = {}
for s in ("h1",) + B.CANCER:
    D = I.S[s]
    z = np.load(G480 / SURF[s], allow_pickle=True)
    add = BD.F.mde_add(z["ctrl_cpm"].astype(float) * LIB / 1e6)[None, :]
    head, block = DC.callset_for(I, s)
    unc = D["tested"] & ~head & ~block
    Ttrue = np.where(D["gate"], D["lfc"], np.nan)
    K = BC.signed_key(I, s)
    key = np.sign(K)
    C4, _ = clean4_lolo(I, s)
    c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
    rows = []
    for f in (1.0, 0.5, 0.35):
        Tp = DB.pool_T(Ttrue, z["p"].astype(float), f)
        a = BD.pred_dnmae(c4h, key, add, Tp, unc)[0]
        rows.append(a)
        vals[(s, f)] = a
    lines.append(f"== {s} (add at lib {LIB:.0f}): c4h − base " + " / ".join(f"{v:+.5f}" for v in rows) + " (f 1 / 0.5 / 0.35)")
m = float(np.mean([vals[(s, 0.35)] for s in B.CANCER]))
lines.append(f"cancer-fold mean at f 0.35 with A/B/C depth: {m:+.5f} (613d's registered value at own depth: +0.00017)")
print("\n".join(lines))
(HERE / "work/block613d3.out").write_text("\n".join(lines) + "\n")
