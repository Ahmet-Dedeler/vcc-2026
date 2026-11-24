"""Q3 — asymmetric shrinkage of the anti-informative UP half.

Q2 measured, on the repaired harness, over the reference's DE pool:
    our DOWN calls  0.6666 accurate   (n=5794)
    our UP   calls  0.4742 accurate   (n=5515)   <- worse than chance
    reference DOWN-rate 0.5979, ours 0.5123

`nmae` in our regime is ~ 1 - <sign(lfc_ref), lfc_hat> / <|lfc_ref|>, so magnitude carried on
a wrong-sign gene is charged twice: it fails to cancel |r| AND adds |h|.  Damping the half of
our predictions that is anti-informative should therefore move `nmae` directly, and demote
those genes in any |lfc|-keyed ordering at the same time.

lam = 1.0 is production.  lam = 0 zeroes every up-prediction.  Deployable: the transform
reads only our own sign, never the reference.

COST NOT PRICED HERE: this changes the mean profile, so `pds` and `mse` both move and must be
measured separately before any arm is built.  DE-bench prices the three DE members only.
"""
import numpy as np
from debench import Surface, overall_delta

s = Surface("../095-clean-harness/de_thinned.npz")
base_nmae, base_reach = s.nmae(), s.reach()
print(f"production: nmae {base_nmae:.5f}   reach {base_reach:.5f}\n")
print(f"{'lam':>6} {'nmae':>9} {'Δnmae':>9} {'reach':>9} {'Δreach':>9} {'sign acc':>9}")
for lam in (1.0, 0.9, 0.75, 0.5, 0.25, 0.0):
    L = np.where(s.pred_lfc < 0, s.pred_lfc, s.pred_lfc * lam)
    key = [np.abs(L[p]) for p in range(len(s.perts))]
    n, r = s.nmae(pred_lfc=L), s.reach(key=key)
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        num += int((np.sign(L[p, i]) == np.sign(s.real_lfc[p, i])).sum()); den += i.size
    print(f"{lam:6.2f} {n:9.5f} {n-base_nmae:+9.5f} {r:9.5f} {r-base_nmae*0:+9.5f} {num/den:9.4f}")

print("\nsame family, but scaling the informative DOWN half UP instead (mu), lam fixed at 1")
for mu in (1.0, 1.25, 1.5, 2.0, 3.0):
    L = np.where(s.pred_lfc < 0, s.pred_lfc * mu, s.pred_lfc)
    print(f"  mu={mu:4.2f}  nmae {s.nmae(pred_lfc=L):9.5f}  Δ {s.nmae(pred_lfc=L)-base_nmae:+9.5f}")

print("\nreach for the |lfc| key under asymmetric damping (the ordering effect alone)")
for lam in (1.0, 0.5, 0.0):
    L = np.where(s.pred_lfc < 0, s.pred_lfc, s.pred_lfc * lam)
    r = s.reach(key=[np.abs(L[p]) for p in range(len(s.perts))])
    print(f"  lam={lam:4.2f}  reach {r:9.5f}  Δoverall vs production order {overall_delta(d_reach=r-base_reach):+9.5f}")
