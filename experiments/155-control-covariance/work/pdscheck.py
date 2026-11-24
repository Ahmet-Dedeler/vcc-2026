import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bench import Surface, SURFACE, FROZEN104
s = Surface(SURFACE)
z = np.load(FROZEN104, allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z["genes"])}
t104 = {str(t): i for i, t in enumerate(z["targets"])}
gcol = np.array([g104.get(g, -1) for g in s.genes]); gok = gcol >= 0
prow = np.array([t104.get(p, -1) for p in s.perts]); pok = prow >= 0
P = np.zeros_like(s.pred_lfc)
P[np.ix_(pok, gok)] = z["baseline"][np.ix_(prow[pok], gcol[gok])]
s.pred_lfc = P
print("pds on FULL gated axis (all 10261 genes):", round(s.pds(), 4))
print("pds exclude_panel=False:", round(s.pds(exclude_panel=False), 4))
for k in ("agreement", "full"):
    Q = np.zeros_like(s.pred_lfc); Q[np.ix_(pok, gok)] = z[k][np.ix_(prow[pok], gcol[gok])]
    s.pred_lfc = Q; print(f"pds {k}:", round(s.pds(), 4))
