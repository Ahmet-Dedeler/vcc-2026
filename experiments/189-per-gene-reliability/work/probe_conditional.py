"""The within-|H|-decile theta gap: is it real, is it expression, and can it be exploited?

`probe_mechanism.py` found that inside the FIRST |H| decile of the reference DE pool -- exactly
where `reach` reads -- the top theta tercile signs 0.6623 against 0.5911 for the bottom, a
+0.0712 gap, while spearman(theta, gene-mean |H|) is +0.031.  Three questions:

  1. Is the gap real?  Bootstrap over TARGETS + a permutation null on theta.
  2. Is it expression in disguise?  The same statistic keyed on control expression.
  3. Can it be exploited?  A key that reorders ONLY inside the top-X% of |H| per target, and a
     q_w weight applied only there.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness189 import Surface189, check, paired_bootstrap
from members189 import precision, topn_call

T = Surface189(); check(T)
absH = np.abs(T.H)
ze = (T.b_ctrl - T.b_ctrl.mean())/T.b_ctrl.std()
base = np.log(np.maximum(absH, 1e-30))
ref = T.reach_vec(absH)

def decile_gap(key, decile=0, seed=None):
    """Per-target (lo-tercile acc, hi-tercile acc, n) inside one |H| decile."""
    k = key if seed is None else np.random.default_rng(seed).permutation(key)
    out = []
    for p in range(len(T.perts)):
        idx = np.flatnonzero(T.pool[p])
        if idx.size < 30: continue
        o = idx[np.argsort(-absH[p, idx])]
        e = (np.arange(11)*o.size)//10
        blk = o[e[decile]:e[decile+1]]
        t = blk.size//3
        if t == 0: continue
        q = blk[np.argsort(k[blk])]
        out.append((T.correct[p, q[:t]].mean(), T.correct[p, q[-t:]].mean(), t))
    a = np.array(out)
    return a

print(f"\n1/2 — the within-|H|-decile tercile gap, decile 1 (the front of the reach walk)")
print(f"{'key':<24}{'lo acc':>9}{'hi acc':>9}{'gap':>9}{'boot SE':>10}{'P>0':>7}"
      f"{'null mean':>11}{'null sd':>9}{'z':>7}")
for name, key in (("theta (the prior)", T.theta), ("control expression", ze),
                  ("gene-mean |H| itself", absH.mean(0))):
    a = decile_gap(key)
    d = a[:,1]-a[:,0]
    g, se, lo, hi, pb = paired_bootstrap(a[:,1], a[:,0])
    nul = [ (lambda x: x[:,1].mean()-x[:,0].mean())(decile_gap(key, seed=4000+s))
            for s in range(30) ]
    print(f"{name:<24}{a[:,0].mean():>9.4f}{a[:,1].mean():>9.4f}{g:>+9.4f}{se:>10.4f}"
          f"{pb:>7.3f}{np.mean(nul):>+11.4f}{np.std(nul,ddof=1):>9.4f}"
          f"{(g-np.mean(nul))/np.std(nul,ddof=1):>+7.2f}")

print(f"\n3 — exploit: reorder by theta only inside the top-X% of |H| per target")
print(f"{'arm':<44}{'reach':>9}{'d vs |H|':>10}{'SE':>8}{'P(bet)':>8}{'prec@50':>9}{'q_w':>9}")
K0, q0 = T.mass_ratio(), T.qw()
print(f"{'|H| (incumbent)':<44}{0.23533:>9.5f}{0.0:>+10.5f}{0.0:>8.5f}{'--':>8}"
      f"{precision(T, topn_call(T, absH, 50))[0]:>9.4f}{q0:>9.5f}")
for frac in (0.02, 0.05, 0.10, 0.25):
    for lam in (0.5, 1.0, 2.0):
        # rank |H| within each target; apply the theta tilt only to the top `frac`
        key = absH.copy()
        for p in range(len(T.perts)):
            n = key.shape[1]
            cut = np.quantile(absH[p], 1-frac)
            sel = absH[p] >= cut
            key[p, sel] = np.exp(base[p, sel] + lam*T.theta[sel])
            # keep the tilted block strictly above the rest
            key[p, sel] += 1e3
        r = T.reach(key)
        d, se, lo, hi, pb = paired_bootstrap(T.reach_vec(key), ref)
        w = np.ones(len(T.genes))
        pr,_ = precision(T, topn_call(T, key, 50))
        print(f"{f'top {frac:.0%} of |H| reordered by exp({lam:+.1f} theta)':<44}"
              f"{r:>9.5f}{d:>+10.5f}{se:>8.5f}{pb:>8.3f}{pr:>9.4f}{'--':>9}")
