"""810 dev: why decomp810's KA-706 differs from dirswap706-abc.json on ABC/k562_gwps (local, through gate602, ~2 GB)."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/810-strong-direction"))
sys.path.insert(0, str(ROOT / "experiments/706-strong-tier-budget"))
import numpy as np
import decomp810 as D
import engine801 as E
import fold522 as FD
import dirswap706 as DS
held = sys.argv[1] if len(sys.argv) > 1 else "k562_gwps"
F = E.load(f"ABC/{held}", check=False)
S = F.S
z706 = FD.zsrc(S)
rows = [DS.donor_rows(s, S.targets)[0] for s in S.sources]
U706 = np.zeros_like(S.F.fST)
for R in rows:
    nr = np.linalg.norm(R, axis=1, keepdims=True)
    U706 += np.where(nr > 0, R / np.maximum(nr, 1e-300), 0.0)
U706[:, ~S.F.nonpanel] = 0.0
U810 = D.unit_mean([D.raw_rows(b, F.targets)[0] for b in F.sources], F.nonpanel)
print("sources", S.sources, F.sources, "U diff", np.abs(U706 - U810).max())
have = np.linalg.norm(U706, axis=1) > 0
tz = np.full(S.n, -1); tz[have] = DS.FD_tiers(z706[have])
print("n_sel", int((tz == 2).sum()))
wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
s0, a0 = S.score(wp)
fin0 = S.F.fST.copy()
nf, nu = np.linalg.norm(fin0, axis=1), np.linalg.norm(U706, axis=1)
Fn = fin0.copy()
for i in np.flatnonzero((tz == 2) & (nf > 0) & (nu > 0)):
    m = 0.75 * fin0[i] / nf[i] + 0.25 * U706[i] / nu[i]
    Fn[i] = nf[i] * m / np.linalg.norm(m)
S.F.fST = Fn
s1, a1 = S.score(wp)
S.F.fST = fin0
live = np.flatnonzero(a0["d"] > 0)
print("706 path dX", 100 * (DS.X_of(np.asarray(a1["x"]), np.asarray(a1["a"]), a0["d"], live) / DS.X_of(np.asarray(a0["x"]), np.asarray(a0["a"]), a0["d"], live) - 1), "dpds", s1["pds"] - s0["pds"])
pct = E.named(F, "pct")
e0, b0 = E.price(F, pct)
e1, b1 = E.price(F, E.Arm(w=pct.w, q=1.75, fin=Fn))
print("w equal", np.abs(pct.w - wp).max(), "base cos", s0["cos_agg"], e0["cos_agg"])
print("engine dX", 100 * (E.cos_agg_arr(b1["x"], b1["a"], b1["d"], live) / E.cos_agg_arr(b0["x"], b0["a"], b0["d"], live) - 1), "dpds", e1["pds"] - e0["pds"])
