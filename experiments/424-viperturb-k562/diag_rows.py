"""424 diagnostic: VIPerturb's prepared panel rows vs the production donors' (E062 space, panel columns zeroed)."""
import sys
sys.dont_write_bytecode = True
import numpy as np
from lib424 import MAIN, Donors, PANEL, OFFICIAL, vip_full_rows, pre_normrel, arm_view
pos = {g: i for i, g in enumerate(OFFICIAL)}
pp = np.asarray([pos[g] for g in PANEL])
dn = Donors(list(PANEL), MAIN)
Fv, Av, au = vip_full_rows()
F = {s: dn.F[s].copy() for s in MAIN}; F["vip"] = Fv.copy()
A = {s: dn.A[s] for s in MAIN}; A["vip"] = Av
for s in F: F[s][:, pp] = 0.0
u = lambda M: M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-300)
allc = A["k562_gwps"] & A["orion_hct116"] & A["orion_hek293t"] & A["vip"]
print("targets carried by all four:", allc.sum())
names = list(F)
for s in names:
    print(f"{s:14s} median row norm (panel, carried) {np.median(np.linalg.norm(F[s][A[s]], axis=1)):.3f}  normrel med {dn.med.get(s, float('nan')):.3f}")
print("per-target cosine (median over targets carried by all four), same target / different-target mean:")
for i, a in enumerate(names):
    for b in names[i + 1:]:
        S = u(F[a][allc]) @ u(F[b][allc]).T
        n = S.shape[0]
        print(f"  {a:14s} x {b:14s} same {np.median(np.diag(S)):+.4f} (mean {np.diag(S).mean():+.4f})  diff {S[~np.eye(n, dtype=bool)].mean():+.4f}")
# mean norms: base mean of 3 vs mean of 4 on the all-four targets
m3 = np.mean([F[s][allc] for s in MAIN], axis=0); m4 = np.mean([F[s][allc] for s in names], axis=0)
print("median |mean of 3| %.3f  |mean of 4| %.3f  ratio^2 %.3f" % (np.median(np.linalg.norm(m3, axis=1)), np.median(np.linalg.norm(m4, axis=1)),
      (np.median(np.linalg.norm(m4, axis=1)) / np.median(np.linalg.norm(m3, axis=1))) ** 2))
# the common component: mean over targets of each donor's rows
for s in names:
    c = F[s][A[s]].mean(0)
    print(f"{s:14s} |mean over panel rows| {np.linalg.norm(c):.3f}  share of row energy {np.linalg.norm(c)**2 / np.mean(np.linalg.norm(F[s][A[s]], axis=1)**2):.3f}")
cm = {s: F[s][A[s]].mean(0) for s in names}
for i, a in enumerate(names):
    for b in names[i + 1:]:
        print(f"  cos(common {a}, common {b}) {float(u(cm[a][None])[0] @ u(cm[b][None])[0]):+.3f}")
