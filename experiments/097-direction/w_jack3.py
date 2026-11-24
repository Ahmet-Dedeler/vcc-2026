"""Jackknife error bars for the three-donor arms, against each harness's two-donor base."""
import numpy as np
from w_arms import key_normrel, soft
from wkit import Ctx, uniform

def rows(c, E):
    dh = c.dh_exact(E)
    return np.einsum("ij,ij->i", dh, c.d) * c.w, np.einsum("ij,ij->i", dh, dh) * c.w

for harness, extra in (("h1-test", "orion_hct116"), ("orion", "h1_2025")):
    c2 = Ctx(harness); b2 = c2.run(uniform); xb, ab = rows(c2, b2["E"])
    nb = np.linalg.norm(b2["E"], axis=1); n = len(c2.perts); del c2
    c = Ctx(harness, extra_donors=(extra,))
    keys = key_normrel(c)
    E3 = c.effects(soft(c, keys, -1.5))
    Eu = c.effects(uniform)
    # dir-only must preserve the row norms of the base WITHIN this donor set, i.e. the
    # 3-donor uniform -- using the 2-donor norms measures a different arm entirely.
    n3 = np.linalg.norm(Eu, axis=1)
    arms = {"a7 3-donor uniform": Eu,
            "a7 3-donor + normrel^-1.5": E3,
            "a7 3-donor + normrel^-1.5 dir":
                E3 * (n3/np.maximum(np.linalg.norm(E3,axis=1),1e-30))[:,None]}
    print(f"\n##### {harness} + {extra}  (n={n}; {sum(p in c.model.pert_index[extra] for p in c.perts)} rows carried by the new donor)")
    Wu = None
    for lbl, Ea in arms.items():
        xa, aa = rows(c, Ea)
        full = np.sqrt((xa.sum()**2/aa.sum()) / (xb.sum()**2/ab.sum()))
        js = np.array([np.sqrt(((xa.sum()-xa[i])**2/(aa.sum()-aa[i]))
                               / ((xb.sum()-xb[i])**2/(ab.sum()-ab[i]))) for i in range(n)])
        se = np.sqrt((n-1)/n*((js-js.mean())**2).sum())
        if Wu is None: Wu = (xa.sum()**2/aa.sum())
        rel = np.sqrt((xa.sum()**2/aa.sum())/Wu)
        print(f"  {lbl:<32} vs 2-donor {full:.4f} +/- {1.96*se:.4f}   "
              f"vs 3-donor uniform {rel:.4f}")
    del c
