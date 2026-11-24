"""Why does an arm with WORSE signs post 2x the reach?  Decompose R_p = k*/n_real."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bench import load_frame
s, inc, gdst, pok = load_frame()
def frame(M):
    P = np.zeros_like(s.pred_lfc); P[:, gdst] = M[:, gdst]; return P
def dbl(M):
    M = M - M.mean(0, keepdims=True); return M - M.mean(1, keepdims=True)
B0 = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "coef_k0.npz"), allow_pickle=True)["beta"].astype(np.float64)
arms = {"incumbent": inc, "+beta k=0": frame(dbl(B0))}
for name, P in arms.items():
    s.pred_lfc = P; ks=[]; nr=[]; Rs=[]
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0: continue
        o = np.lexsort((s.genes[pool], -np.abs(P[p, pool]))); q = pool[o]
        q = q[s.adj[p, q]]
        if q.size == 0: ks.append(0); nr.append(int(s.real_sig[p].sum())); Rs.append(0.0); continue
        ok = np.sign(P[p, q]) == np.sign(s.real_lfc[p, q])
        pur = np.cumsum(ok)/np.arange(1, q.size+1); d = np.flatnonzero(pur >= 0.9)
        k = int(d[-1]+1) if d.size else 0
        ks.append(k); nr.append(int(s.real_sig[p].sum())); Rs.append(k/int(s.real_sig[p].sum()))
    ks=np.array(ks); nr=np.array(nr); Rs=np.array(Rs)
    small = nr <= np.median(nr)
    print(f"{name:12} mean R {Rs.mean():.4f}  median k* {np.median(ks):5.1f}  mean k* {ks.mean():6.2f}  "
          f"frac k*=0 {(ks==0).mean():.3f}  R on small-n_real {Rs[small].mean():.4f} vs large {Rs[~small].mean():.4f}")
print(f"n_real: median {np.median(nr):.0f}, min {nr.min()}, max {nr.max()}")
