#!/usr/bin/env python3
"""613d7: a THIRD known-answer row for 613d. KGID in its v1 form (idsign612: identity alone where CLEAN4 is undefined, the
identity-only extension; key-gated) vs c4h on h1. The real emitter has idkg − c4h −0.00315 / −0.00019 / +0.00067 (3 seeds).
Known-answer row only. Output: work/block613d7.out."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613c as BC  # noqa: E402
import block613d as BD  # noqa: E402
import decomp639 as DC  # noqa: E402
import ds612b as DB  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
I = Inst()
prof = DC.flex_profiles()
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
Iv = DC.identity(prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")], list(D["G"]))
zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
sc1 = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan_to_num(zi))      # v1: identity fills CLEAN4's gaps
v1 = np.where(np.isfinite(sc1) & (sc1 != 0), np.where(sc1 > 0, -1.0, 1.0), 0.0)
weak = np.zeros(ak.shape, bool)
for i in range(ak.shape[0]):
    nz = ak[i] > 0
    weak[i] = (~nz) | (ak[i] < (np.quantile(ak[i][nz], 0.5) if nz.any() else 0.0))
kg1 = np.where(weak, np.broadcast_to(v1, weak.shape), 0.0)
real = {1.0: -0.00315, 0.5: -0.00019, 0.35: 0.00067}
out = []
for f in (1.0, 0.5, 0.35):
    Tp = DB.pool_T(Ttrue, z["p"].astype(float), f)
    a = BD.pred_dnmae(c4h, key, add, Tp, unc)[0]
    b = BD.pred_dnmae(kg1, key, add, Tp, unc)[0]
    out.append(f"h1 f {f:g}: KGID-v1 − c4h {b - a:+.5f} (real emitter idkg − c4h {real[f]:+.5f}; |err| {abs(b - a - real[f]):.5f})")
print("\n".join(out))
(HERE / "work/block613d7.out").write_text("\n".join(out) + "\n")
