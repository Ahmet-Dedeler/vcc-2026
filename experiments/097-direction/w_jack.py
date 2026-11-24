"""Sampling error bar on sqrt(W_arm/W_base) over PERTURBATIONS -- the error bar that matters.

The emitter noise floor (w_extra) says how reproducible a number is; this says whether 85 or
132 perturbations are enough to see it.  X and A are plain sums of per-row terms, so a
delete-one-perturbation jackknife is exact and free.
"""
import numpy as np
from w_arms import key_normrel, soft
from wkit import Ctx, fixed, uniform

def rows(c, E):
    dh = c.dh_exact(E)
    return np.einsum("ij,ij->i", dh, c.d) * c.w, np.einsum("ij,ij->i", dh, dh) * c.w

for name in ("h1-test", "orion"):
    c = Ctx(name)
    b = c.run(uniform)
    nb = np.linalg.norm(b["E"], axis=1)
    keys = key_normrel(c)
    xb, ab = rows(c, b["E"])
    n = len(c.perts)
    print(f"\n##### {name}  (n={n} perturbations)")
    arms = {"a3 share 0.35": c.effects(fixed({"k562_gwps":0.35,"orion_hek293t":0.65})),
            "a3 share 0.60": c.effects(fixed({"k562_gwps":0.60,"orion_hek293t":0.40})),
            "a4 basal": c.effects(fixed({k: float(np.exp(10*v)) for k,v in
                                         c.model.basal_weights(c.b_ctrl).items()})),
            "a6 normrel^-1.5": c.effects(soft(c, keys, -1.5))}
    E = c.effects(soft(c, keys, -1.5))
    arms["a6 normrel^-1.5 dir"] = E * (nb/np.maximum(np.linalg.norm(E,axis=1),1e-30))[:,None]
    for lbl, Ea in arms.items():
        xa, aa = rows(c, Ea)
        full = np.sqrt((xa.sum()**2/aa.sum()) / (xb.sum()**2/ab.sum()))
        js = np.array([np.sqrt(((xa.sum()-xa[i])**2/(aa.sum()-aa[i]))
                               / ((xb.sum()-xb[i])**2/(ab.sum()-ab[i]))) for i in range(n)])
        se = np.sqrt((n-1)/n * ((js - js.mean())**2).sum())
        print(f"  {lbl:<24} mult {full:.4f}  jackknife SE {se:.4f}  "
              f"=> {full:.3f} +/- {1.96*se:.3f} (95%)  z={(full-1)/se:+.1f}")
    del c
