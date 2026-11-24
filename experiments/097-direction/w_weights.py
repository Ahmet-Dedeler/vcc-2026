"""Rule 3: print the intervention itself, not just its effect."""
import sys
import numpy as np
from w_arms import key_normrel, soft
from wkit import Ctx, uniform

for name in ("h1-test", "orion"):
    c = Ctx(name)
    keys = key_normrel(c)
    dn = list(c.model.names)
    f = soft(c, keys, -1.5)
    Wt = np.array([[f(p)[n] for n in dn] for p in c.perts])
    Wt = Wt / Wt.sum(1, keepdims=True)
    print(f"\n##### {name}: normrel^-1.5 normalised donor weights over {len(c.perts)} rows")
    for i, n in enumerate(dn):
        q = np.percentile(Wt[:, i], [0, 10, 25, 50, 75, 90, 100])
        print(f"  {n:<16} min {q[0]:.3f} p10 {q[1]:.3f} p25 {q[2]:.3f} med {q[3]:.3f} "
              f"p75 {q[4]:.3f} p90 {q[5]:.3f} max {q[6]:.3f}")
    k = keys[dn[0]] / np.maximum(keys[dn[1]], 1e-12)
    print(f"  donor effect-norm RATIO k562/hek per row: median {np.median(k):.3f}  "
          f"p10 {np.percentile(k,10):.3f}  p90 {np.percentile(k,90):.3f}")
    print(f"  rows where the weight moves >10pp from 0.5: "
          f"{int(np.sum(np.abs(Wt[:,0]-0.5) > 0.10))}/{len(c.perts)}")
    b = c.run(uniform); E1 = c.effects(f)
    nb = np.linalg.norm(b['E'], axis=1)
    Ed = E1 * (nb/np.maximum(np.linalg.norm(E1,axis=1),1e-30))[:,None]
    cs = np.einsum('ij,ij->i', Ed, b['E'])/np.maximum(np.linalg.norm(Ed,axis=1)*nb,1e-30)
    print(f"  cos(E_arm_row, E_base_row): median {np.median(cs):.4f}  min {cs.min():.4f}")
    print(f"  median ||E|| base {np.median(nb):.4f} -> arm(dir-only, by construction) "
          f"{np.median(np.linalg.norm(Ed,axis=1)):.4f}")
    xb = np.einsum('ij,ij->i', c.dh_exact(b['E']), c.d)*c.w
    xa = np.einsum('ij,ij->i', c.dh_exact(Ed), c.d)*c.w
    print(f"  per-row X: rows improved {int(np.sum(xa>xb))}/{len(c.perts)}; "
          f"top-quartile-|X| rows carry {100*np.sort(np.abs(xb))[-len(c.perts)//4:].sum()/np.abs(xb).sum():.1f}% of |X|")
    del c
