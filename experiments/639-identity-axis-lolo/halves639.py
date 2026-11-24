#!/usr/bin/env python3
"""639c - is the identity gain on h1 a shared-control-noise artifact? Identity from control half A, truth vs half B
(independent per-gene halves, 635's construction); CLEAN4 + lambda*identity vs CLEAN4 on gate genes, paired."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo
import ident639 as M
z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(z["keep"], bool); G = z["genes"].astype(str)[keep]; T = z["targets"].astype(str)
q, lfc, ref_lin, cs = z["q"], z["lfc"], np.asarray(z["ref_lin"], float), z["cs"]
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (G[None, :] != T[:, None])
mp = np.where(np.isfinite(lfc), ref_lin[None, :] * np.exp2(lfc), np.nan)
rng = np.random.default_rng(639); n = cs.shape[1]; h = n // 2
mA = np.empty(cs.shape[0]); mB = np.empty(cs.shape[0])
for i in range(0, cs.shape[0], 1000):
    X = np.expm1(cs[i:i+1000]); idx = np.argsort(rng.random(X.shape), axis=1); Xs = np.take_along_axis(X, idx, axis=1)
    mA[i:i+1000] = Xs[:, :h].mean(1); mB[i:i+1000] = Xs[:, h:].mean(1)
truthB = np.log2((mp + 1e-9) / (mB[None, :] + 1e-9))
# other profiles (never h1): the three fold surfaces + A/B/C
prof = {}
for s, f in SURF.items():
    if s == "h1": continue
    zz = np.load(G480 / f, allow_pickle=True); prof[s] = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
prof.update(M.abc_profiles())
gi = {g: i for i, g in enumerate(G)}
cg = sorted(set(G) & set.intersection(*(set(p) for p in prof.values())))
cols = np.array([gi[g] for g in cg])
others = np.mean([np.log1p(np.array([prof[k][g] for g in cg])) for k in prof], 0)
I = Inst(); C4, _ = clean4_lolo(I, "h1")
assert (I.S["h1"]["G"] == G).all()
zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
for name, mctrl, truth in (("identity from FULL controls, truth vs full (shared noise)", ref_lin, lfc),
                           ("identity from half A, truth vs half B (clean)", mA, truthB),
                           ("identity from half B, truth vs half B (shared half)", mB, truthB)):
    iv = np.full(len(G), np.nan); iv[cols] = np.log1p(mctrl[cols]) - others
    zi = (iv - np.nanmean(iv)) / np.nanstd(iv)
    both = np.isfinite(zc) & np.isfinite(zi); m = gate & both[None, :] & np.isfinite(truth) & (truth != 0)
    dn = truth < 0
    acc = lambda sc: float((((sc > 0)[None, :] == dn) & m).sum() / m.sum())
    print(f"{name:58s} C4 {acc(zc):.4f}  C4+0.5I {acc(zc + 0.5*zi):.4f}  C4+1I {acc(zc + zi):.4f}  I alone {acc(zi):.4f}  (pairs {int(m.sum())})")
