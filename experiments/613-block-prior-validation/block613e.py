#!/usr/bin/env python3
"""613e: a held-out PREMISE check for 630, rank-neutral UP magnitude on uncalled UP genes ([d4cd5f]). It is not a price.
630 raises the delivered |lfc| of UNCALLED genes whose delivered sign is UP. On a gate gene that helps only if the UP
sign is right more often than wrong (and |T| leaves headroom above the add). This measures, on 636's surfaces at 521's
pools, among uncalled gate genes whose delivered sign is UP:
  acc_up    = P(real lfc > 0);
  headroom  = the median |T| / add, and the share of cells with |T| > 2·add.
The delivered sign is the key's (no hook, 446bj2z's form) and, beside it, c4h's (the plain hook where CLEAN4 is
confident). Committed before its first run. Reading: acc_up > 0.5 on the cancer folds supports 630's premise off H1;
< 0.5 means it would hurt there. No rule is attached; it is a second reading for d4cd5f's October case (11.1c.2).
Output: work/block613e.out.
"""
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
I = Inst()
lines = []
for s in ("h1",) + B.CANCER:
    D = I.S[s]
    z = np.load(G480 / SURF[s], allow_pickle=True)
    add = BD.F.mde_add(z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6)[None, :]
    head, block = DC.callset_for(I, s)
    unc = D["tested"] & ~head & ~block
    Ttrue = np.where(D["gate"], D["lfc"], np.nan)
    key = np.sign(BC.signed_key(I, s))
    C4, _ = clean4_lolo(I, s)
    c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
    hooked = np.where(np.broadcast_to(c4h, key.shape) != 0, np.broadcast_to(c4h, key.shape), key)
    p = z["p"].astype(float)
    parts = []
    for f in (1.0, 0.35):
        Tp = DB.pool_T(Ttrue, p, f)
        G = np.isfinite(Tp) & unc
        for lab, sg in (("key", key), ("c4h", hooked)):
            m = G & (sg > 0) & (add > 0)
            acc = float((Tp[m] > 0).mean()) if m.any() else float("nan")
            ratio = np.abs(Tp[m]) / np.broadcast_to(add, Tp.shape)[m]
            parts.append(f"f {f:g} {lab}: acc_up {acc:.3f} n/t {m.sum() / len(D['T']):.0f} |T|/add med {np.median(ratio):.1f} "
                         f">2add {float((ratio > 2).mean()):.2f}")
    lines.append(f"== {s}: " + " | ".join(parts))
print("\n".join(lines))
(HERE / "work/block613e.out").write_text("\n".join(lines) + "\n")
