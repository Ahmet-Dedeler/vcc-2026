#!/usr/bin/env python3
"""613d6 (descriptive): ungated CLEAN4 at tau 0 (c4all0: CLEAN4 sign on every covered uncalled gene) vs the ungated v3 identity hook and
c4h by 613d's add rule: is the key gate needed once the source is v3? Output:
work/block613d6.out."""
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
prof = DC.flex_profiles()
lines = []
for s in ("h1",) + B.CANCER:
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
    if s == "h1":
        own, refs = prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
    else:
        own, refs = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float))), [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
    Iv = DC.identity(own, refs, list(D["G"]))
    zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    sc = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan)
    v3 = np.where(np.isfinite(sc) & (sc != 0), np.where(sc > 0, -1.0, 1.0), 0.0)
    weak = np.zeros(ak.shape, bool)
    for i in range(ak.shape[0]):
        nz = ak[i] > 0
        weak[i] = (~nz) | (ak[i] < (np.quantile(ak[i][nz], 0.5) if nz.any() else 0.0))
    c4kg = np.where(weak, np.broadcast_to(c4h, weak.shape), 0.0)
    kgid = np.where(weak, np.broadcast_to(v3, weak.shape), 0.0)
    parts = []
    for f in (1.0, 0.5, 0.35):
        Tp = DB.pool_T(Ttrue, z["p"].astype(float), f)
        a = BD.pred_dnmae(c4h, key, add, Tp, unc)[0]
        b = BD.pred_dnmae(c4kg, key, add, Tp, unc)[0]
        c = BD.pred_dnmae(kgid, key, add, Tp, unc)[0]
        u = BD.pred_dnmae(np.broadcast_to(v3, weak.shape), key, add, Tp, unc)[0]
        c40 = np.where(np.isfinite(C4) & (C4 != 0.5), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
        w = BD.pred_dnmae(np.broadcast_to(c40, weak.shape), key, add, Tp, unc)[0]
        parts.append(f"f {f:g}: c4all0−c4h {w - a:+.5f}  ungated v3−c4h {u - a:+.5f}  (identity over tau-0 CLEAN4 {u - w:+.5f})")
    lines.append(f"== {s}: " + " | ".join(parts))
print("\n".join(lines))
(HERE / "work/block613d6.out").write_text("\n".join(lines) + "\n")
