#!/usr/bin/env python3
"""613e2: 613e's footprint restated in nmae's own units ([d4cd5f]'s point 3: cells per target conflate n_real). For each
scored target (>= 10 pool genes), the share of its Σ|T| that sits on UNCALLED gate genes whose delivered sign is UP
(key; c4h beside it), then the mean over targets. nmae normalises per target, so this share bounds a lever on those genes.
Also the head / block / uncalled split of the gated mass for reference. Output: work/block613e2.out."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613 as B  # noqa: E402
import block613c as BC  # noqa: E402
import decomp639 as DC  # noqa: E402
import ds612b as DB  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
I = Inst()
lines = []
for s in ("h1",) + B.CANCER:
    D = I.S[s]
    z = np.load(G480 / SURF[s], allow_pickle=True)
    head, block = DC.callset_for(I, s)
    unc = D["tested"] & ~head & ~block
    Ttrue = np.where(D["gate"], D["lfc"], np.nan)
    key = np.sign(BC.signed_key(I, s))
    C4, _ = clean4_lolo(I, s)
    c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
    hooked = np.where(np.broadcast_to(c4h, key.shape) != 0, np.broadcast_to(c4h, key.shape), key)
    parts = []
    for f in (1.0, 0.35):
        Tp = DB.pool_T(Ttrue, z["p"].astype(float), f)
        A = np.where(np.isfinite(Tp), np.abs(Tp), 0.0)
        tot = A.sum(1)
        ok = (np.isfinite(Tp).sum(1) >= 10) & (tot > 0)
        sh = lambda M: float(np.mean((A * M).sum(1)[ok] / tot[ok]))  # noqa: E731
        parts.append(f"f {f:g}: head {sh(head):.2f} block {sh(block):.2f} uncalled {sh(unc):.2f} | uncalled-UP key {sh(unc & (key > 0)):.3f} "
                     f"c4h {sh(unc & (hooked > 0)):.3f} (targets {int(ok.sum())})")
    lines.append(f"== {s}: " + " || ".join(parts))
print("\n".join(lines))
(HERE / "work/block613e2.out").write_text("\n".join(lines) + "\n")
