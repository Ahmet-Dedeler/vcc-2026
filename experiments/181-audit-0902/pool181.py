"""181 — does `pds` delta scale with the DECOY-POOL SIZE?  178 concluded 'h1-test understates
pds costs of un-whitening arms ~6x; PLOLO is the pds instrument' from PLOLO -0.0245 vs 173's
h1-test -0.0042.  But h1-test ranks 85 decoys and PLOLO ranks 269-294 (the board ranks 299).
162 already noted the same effect ('magnitudes on the Orion pair are 4-5x larger than 156's,
consistent with the decoy pool going from 132 to ~290').  If pds deltas grow with the pool,
the '6x' is a frame difference and must NOT be applied as a general instrument multiplier.

Test: take a fixed, realistic perturbation of the prediction (flip the sign of the smallest-|lfc|
entries, i.e. 173's sign-repair shape) and measure |delta pds| against the number of
perturbations in the ranking, by subsampling h1-test's 85.
What this reads if the mechanism is absent: a flat line.
"""
import numpy as np, sys
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
S = Surface('experiments/177-resample-emitter/de177.KA_a15_raw.npz')
P0 = S.pred_lfc.copy()

# 173's shape: flip the sign of the 15.7% smallest-|lfc| gated entries
def make_flip(frac=0.157, seed=0):
    P = P0.copy()
    for p in range(P.shape[0]):
        v = np.abs(P[p]); k = int(frac*np.isfinite(v).sum())
        idx = np.argsort(v)[:k]
        P[p, idx] *= -1
    return P
PF = make_flip()

rng = np.random.default_rng(0)
print(f"{'n perts':>8} {'pds base':>10} {'pds flip':>10} {'d pds':>10}  (mean of 25 random subsets)")
for n in (10, 20, 40, 60, 85):
    d=[]; b=[]
    for r in range(25 if n < 85 else 1):
        idx = np.sort(rng.choice(85, n, replace=False)) if n < 85 else np.arange(85)
        s = Surface('experiments/177-resample-emitter/de177.KA_a15_raw.npz')
        s.perts=s.perts[idx]; s.real_lfc=s.real_lfc[idx]; s.pred_lfc=s.pred_lfc[idx]
        s.real_sig=s.real_sig[idx]; s.pred_sig=s.pred_sig[idx]; s.adj=s.adj[idx]
        s.self_idx=s.self_idx[idx]
        p0=s.pds(pred_lfc=P0[idx]); p1=s.pds(pred_lfc=PF[idx])
        b.append(p0); d.append(p1-p0)
    print(f'{n:8d} {np.mean(b):10.5f} {np.mean(b)+np.mean(d):10.5f} {np.mean(d):+10.5f}')
print()
print('If |d pds| grows with the number of ranked perturbations, then 173 (85 decoys),')
print('178/PLOLO (269-294) and the board (299) are three different frames, and the ratio')
print('between them is a frame effect, not an instrument error.')
