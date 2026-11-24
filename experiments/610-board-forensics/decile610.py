import sys
from pathlib import Path
import numpy as np
R = Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0, str(R / "experiments/613-block-prior-validation"))
import block613 as B  # noqa
from inst636 import G480, SURF, Inst, clean4_lolo
I = Inst()
for s in ("h1",) + B.CANCER:
    D = I.S[s]; lfc, gate, tested = D["lfc"], D["gate"], D["tested"]
    z = np.load(G480 / SURF[s], allow_pickle=True); q = z["q"].astype(float); cpm = z["ctrl_cpm"].astype(float)
    null = tested & np.isfinite(q) & (q >= 0.05) & (lfc != 0)
    C4, _ = clean4_lolo(I, s); c4dn = C4 > 0.5; cdef = np.isfinite(C4)
    order = np.argsort(-cpm, kind="stable"); rank = np.empty(len(cpm), int); rank[order] = np.arange(len(cpm))
    edges = [0, 2200, 4400, 6600, 8800, 11000, len(cpm)]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        gm = (rank >= a) & (rank < b) & cdef
        m = null & gm[None, :]
        acc = float(((c4dn[None, :] == (lfc < 0)) & m).sum() / max(m.sum(), 1))
        dn = float(((lfc < 0) & m).sum() / max(m.sum(), 1))
        out.append(f"{a}-{b}: {acc:.4f} (down {dn:.3f}, n/t {m.sum() / len(D['T']):.0f})")
    print(f"{s:8s} CLEAN4 null precision by expression rank: " + " | ".join(out))
