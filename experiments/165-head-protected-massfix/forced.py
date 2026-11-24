"""165 — the part of the profile move that NO placement can remove.

Once `sum_g clip(expm1(b + h), 0) == TS` holds, a gene the correction does not touch is
delivered at `expm1(b+D)/TS` where a15 delivered it at `expm1(b+D)/S`.  Its delivered log
effect therefore moves by exactly `+log(S/TS)`, with no free parameter, and because `log1p`
saturates that lift is LARGEST on the highest-expressed genes -- which is where the reach head
lives (81st expression percentile, `place.py`).

Minimum-L2 puts `eps = lam*w`, `lam < 0`, `w = exp(b+h)` on those same genes, so it partly
CANCELS the forced lift.  Head protection removes the cancellation and leaves the bare lift.
This script measures both, deterministically -- no emitter, no seed, no cell-eval2.
"""
import numpy as np
TS = 5.0e4
z = np.load("work/profiles_all.npz", allow_pickle=True)
b = z["b_ctrl"]; keep = z["keep"]
e = np.load("work/effects.npz", allow_pickle=True)
D, gate, mr = e["effect"], e["gate"], e["mass_ratio"]
h = np.load("work/head.npz")["head"]                 # full-axis order, ascending a15 pred p
P0 = np.log1p(TS * z["p_a15"].astype(np.float64)) - b[None, :]
for arm in ("mn", "hp_p2000", "hp_p9000", "eps_tail", "cst"):
    P = np.log1p(TS * z[f"p_{arm}"].astype(np.float64)) - b[None, :]
    d = P - P0
    head = np.array([d[i][h[i][:500]] for i in range(d.shape[0])])
    print(f"{arm:<10} mean delivered-effect move on the top-500 by p: {head.mean():+.6f}"
          f"   |move| {np.abs(head).mean():.6f}")
print(f"\nforced lift log(S/TS) = {np.log(mr).mean():+.6f}  (mass ratio {mr.mean():.6f})")
lin = np.clip(np.expm1(b[None, :] + D), 0.0, None)
forced = D - (np.log1p(lin / mr[:, None]) - b[None, :])
fh = np.array([forced[i][h[i][:500]] for i in range(forced.shape[0])])
print(f"analytic forced move on the top-500 by p (eps = 0 there): {fh.mean():+.6f}")
print(f"analytic forced move over ALL gated genes:                "
      f"{forced[:, gate].mean():+.6f}")

print("\nmean DELIVERED-EFFECT move by control-expression decile (10 = highest):")
expr = b.copy(); g = np.flatnonzero(gate)
q = np.digitize(expr[g], np.percentile(expr[g], np.arange(10, 100, 10)))
print(f"{'decile':<8}" + "".join(f"{i+1:>9}" for i in range(10)))
for arm in ("mn", "hp_p2000", "eps_tail"):
    P = np.log1p(TS * z[f"p_{arm}"].astype(np.float64)) - b[None, :]
    d = (P - P0)[:, g].mean(0)
    print(f"{arm:<8}" + "".join(f"{d[q == i].mean():>+9.5f}" for i in range(10)))
print(f"{'forced':<8}" + "".join(f"{forced[:, g].mean(0)[q == i].mean():>+9.5f}"
                                 for i in range(10)))
