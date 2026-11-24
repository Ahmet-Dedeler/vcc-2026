"""Rule 3: print the intermediate quantities before believing the arm table."""
import sys, numpy as np
sys.path.insert(0,"~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, BASE_DONORS, ALPHA, MIN_CPM
C = Ctx("orion"); U=("k562_gwps","orion_hek293t","kolf2","cd4_rest")
m = C.set_donors(U)
print("per-donor _scale (median ||s||) and _renorm:")
for n in U: print(f"  {n:16s} scale={m._scale[n]:.5g} renorm={m._renorm[n]:.5g} "
                  f"ncols={m.tcols[n].size}")
gate = m.cpm_gate_mask(C.b_ctrl, MIN_CPM) & (C.b_ctrl > 0)
print(f"CPM>5 & expressed gate keeps {int(gate.sum())} / {C.G} genes")
for lbl, dons in (("base2",("k562_gwps","orion_hek293t")),
                  ("base2+kolf2",("k562_gwps","orion_hek293t","kolf2")),
                  ("base2+cd4_rest",("k562_gwps","orion_hek293t","cd4_rest"))):
    w = {n:(1.0 if n in dons else 0.0) for n in U}
    raw = np.vstack([m.effect(p,1.0,0.0,weights=dict(w)) for p in C.perts])
    print(f"\n{lbl}")
    print(f"  raw  median||E||={np.median(np.linalg.norm(raw,axis=1)):.5g} "
          f"median||E[gate]||={np.median(np.linalg.norm(raw[:,gate],axis=1)):.5g} "
          f"gate share of energy={float(np.median((raw[:,gate]**2).sum(1)/(raw**2).sum(1))):.4f}")
    E = C.effects(lambda _p, w=w: dict(w))
    print(f"  post median||E||={np.median(np.linalg.norm(E,axis=1)):.5g} "
          f"median||E[gate]||={np.median(np.linalg.norm(E[:,gate],axis=1)):.5g}")
    dh = C.dh_exact(E)
    print(f"  dh   median||dh||={np.median(np.linalg.norm(dh,axis=1)):.5g}")
    # per-donor contribution norms inside the mixture, before whitening
    for n in dons:
        v = np.zeros(C.G); 
        vals=[]
        for p in C.perts[:40]:
            i = m.pert_index[n].get(p)
            if i is None: continue
            v[:] = 0; v[m.tcols[n]] = m.s[n][i][m.dcols[n]]
            vals.append(np.linalg.norm(m._renorm[n]*v))
        print(f"    {n:16s} median||renorm*s_p|| over 40 rows = {np.median(vals):.5g}")
