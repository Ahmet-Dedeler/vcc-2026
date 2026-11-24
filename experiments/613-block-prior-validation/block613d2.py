#!/usr/bin/env python3
"""613d2: the mechanism row (A1's criterion d) for 613d's add-rule instrument on h1: the key gate on CLEAN4 (c4kg = c4h on
each target's weak-key half only) vs c4h. The real emitter has c4kg − c4h nmae +0.00031 / +0.00124 / +0.00107 at f 1 /
0.5 / 0.35 (3 seeds), i.e. the gate COSTS nmae on H1. Known-answer row only; no idkg3 number is read here."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613d as BD  # noqa: E402
import block613c as BC  # noqa: E402
import decomp639 as DC  # noqa: E402
import ds612b as DB  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
I = Inst()
s = "h1"
D = I.S[s]
z = np.load(G480 / SURF[s], allow_pickle=True)
add = BD.F.mde_add(z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6)[None, :]
head, block = DC.callset_for(I, s)
unc = D["tested"] & ~head & ~block
Ttrue = np.where(D["gate"], D["lfc"], np.nan)
K = BC.signed_key(I, s)
key, ak = np.sign(K), np.abs(K)
C4, _ = clean4_lolo(I, s)
c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
weak = np.zeros(ak.shape, bool)
for i in range(ak.shape[0]):
    nz = ak[i] > 0
    weak[i] = (~nz) | (ak[i] < (np.quantile(ak[i][nz], 0.5) if nz.any() else 0.0))
c4kg = np.where(weak, np.broadcast_to(c4h, weak.shape), 0.0)
p = z["p"].astype(float)
for f in (1.0, 0.5, 0.35):
    Tp = DB.pool_T(Ttrue, p, f)
    a = BD.pred_dnmae(c4h, key, add, Tp, unc)[0]
    b = BD.pred_dnmae(c4kg, key, add, Tp, unc)[0]
    print(f"h1 f {f}: c4kg − c4h {b - a:+.5f} (real emitter {dict([(1.0, '+0.00031'), (0.5, '+0.00124'), (0.35, '+0.00107')])[f]})")
