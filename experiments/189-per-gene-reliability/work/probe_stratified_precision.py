"""The one arm that survived: theta reordering INSIDE |H| strata, priced on call-set precision.

`fid_p = k / max(n_pred_p, n_real_p)` reads sign precision directly, so precision at a fixed
`n_pred` is the transferable quantity (rule 9 forbids reading absolute `fid` off this surface).
Measured with a per-TARGET bootstrap and a permutation null on theta, at four call sizes.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness189 import Surface189, check
from members189 import topn_call

T = Surface189(); check(T)
absH = np.abs(T.H); base = np.log(np.maximum(absH, 1e-30))

def strat_key(theta, frac=0.05, lam=2.0):
    key = absH.copy()
    for p in range(len(T.perts)):
        sel = absH[p] >= np.quantile(absH[p], 1-frac)
        key[p, sel] = np.exp(base[p, sel] + lam*theta[sel]) + 1e3
    return key

def prec_per_target(call):
    S = T.S; out = []
    for p in range(len(T.perts)):
        idx = np.flatnonzero(call[p] & S.adj[p])
        if idx.size == 0: continue
        out.append(((np.sign(S.pred_lfc[p, idx]) == np.sign(S.real_lfc[p, idx])).sum(),
                    idx.size))
    return np.array(out, float)

def boot(a, b, n=20000, seed=7):
    rng = np.random.default_rng(seed); d = a-b
    idx = rng.integers(0, len(d), size=(n, len(d)))
    bt = d[idx].mean(1)
    return float(d.mean()), float(bt.std(ddof=1)), float((bt>0).mean())

print(f"\n{'n_pred':>7}{'|H| prec':>10}{'strat prec':>12}{'delta':>9}{'boot SE':>10}{'P>0':>7}"
      f"{'null mean d':>13}{'null sd':>9}{'z':>7}")
for n in (25, 50, 100, 200, 400):
    ph = prec_per_target(topn_call(T, absH, n))
    base_acc = ph[:,0]/ph[:,1]
    ps = prec_per_target(topn_call(T, strat_key(T.theta), n))
    d, se, pb = boot(ps[:,0]/ps[:,1], base_acc)
    nul = []
    for s in range(30):
        pn = prec_per_target(topn_call(
            T, strat_key(np.random.default_rng(9000+s).permutation(T.theta)), n))
        nul.append(float((pn[:,0]/pn[:,1] - base_acc).mean()))
    print(f"{n:>7}{ph[:,0].sum()/ph[:,1].sum():>10.4f}{ps[:,0].sum()/ps[:,1].sum():>12.4f}"
          f"{d:>+9.4f}{se:>10.4f}{pb:>7.3f}{np.mean(nul):>+13.4f}{np.std(nul,ddof=1):>9.4f}"
          f"{(d-np.mean(nul))/np.std(nul,ddof=1):>+7.2f}")
