"""Q2 — the sign confusion matrix, and whether an asymmetric key beats a magnitude key.

094's arm a1 was built on an HCT116 down-asymmetry and withdrawn when h1-test appeared to
refute it.  That was measured on the DEPTH-CONFOUNDED harness (90.8% down).  This is the
repaired one.  Measured exactly, not through the emitter.
"""
import numpy as np
from debench import Surface, overall_delta

s = Surface("../095-clean-harness/de_thinned.npz")
P = len(s.perts)

rs, ps = [], []
for p in range(P):
    i = np.flatnonzero(s.real_sig[p] & s.adj[p])
    rs.append(np.sign(s.real_lfc[p, i])); ps.append(np.sign(s.pred_lfc[p, i]))
r = np.concatenate(rs); q = np.concatenate(ps)

print(f"reference DE pool: {r.size} (pert, gene) pairs")
print(f"  reference is DOWN   {(r<0).mean():.4f}")
print(f"  we predict  DOWN    {(q<0).mean():.4f}")
print("\nconfusion — of the genes WE call in each direction, how often the reference agrees")
for lab, m in (("we say DOWN", q < 0), ("we say UP  ", q > 0)):
    acc = (r[m] == q[m]).mean()
    print(f"  {lab}: n={m.sum():6d}   accuracy {acc:.4f}")

print("\nreach under asymmetric keys (deterministic)")
keys = {
    "production":                   None,
    "|predicted lfc|":              [np.abs(s.pred_lfc[p]) for p in range(P)],
    "predicted-DOWN first":         [(s.pred_lfc[p] < 0).astype(float) for p in range(P)],
    "DOWN first, then |lfc|":       [(s.pred_lfc[p] < 0) * 1e6 + np.abs(s.pred_lfc[p]) for p in range(P)],
    "|lfc| but DOWN-only eligible": [np.where(s.pred_lfc[p] < 0, np.abs(s.pred_lfc[p]), -1.0) for p in range(P)],
}
base = None
print(f"{'key':32} {'raw reach':>10} {'Δoverall':>10}")
for name, k in keys.items():
    v = s.reach(key=k)
    if base is None: base = v
    print(f"{name:32} {v:10.5f} {overall_delta(d_reach=v-base):+10.5f}")

# per-perturbation down-rate spread: is the asymmetry global or per-target?
rates = []
for p in range(P):
    i = np.flatnonzero(s.real_sig[p] & s.adj[p])
    if i.size >= 10: rates.append((np.sign(s.real_lfc[p, i]) < 0).mean())
rates = np.array(rates)
print(f"\nper-perturbation reference DOWN-rate: median {np.median(rates):.3f}  "
      f"p10 {np.percentile(rates,10):.3f}  p90 {np.percentile(rates,90):.3f}  "
      f"frac of perts >0.5: {(rates>0.5).mean():.3f}")
