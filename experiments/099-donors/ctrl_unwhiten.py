"""099 arm 5c -- THE control the post-whitening CD4 arm needs.

Adding ANY unwhitened vector back onto a whitened matrix partially undoes the operator, and 097
measured that undoing the operator raises the direction by 1.50-2.14x. So before believing that
`+ f * CD4` is about CD4, add back the BASE MIXTURE'S OWN unwhitened effect at the same f and
see how much of the gain that reproduces. If it reproduces it, the arm is a whitening-strength
knob wearing a donor's name -- and 097 already priced that knob as net negative on both surfaces.
"""
from __future__ import annotations
import sys, argparse
import numpy as np, pandas as pd
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT, ALPHA

CD4 = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
FRACS = (0.0, 0.05, 0.10, 0.25, 0.50, 1.0)

ap = argparse.ArgumentParser(); ap.add_argument("harness")
ap.add_argument("--gamma", type=float, default=0.65)
a = ap.parse_args()
C = Ctx(a.harness)
base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116") if d not in C.banned)[:2]
m = C.set_donors(base + CD4)
w0 = {n: (1.0 if n in base else 0.0) for n in m.names}
RAW = np.vstack([m.effect(p, 1.0, 0.0, weights=dict(w0)) for p in C.perts]).astype(np.float64)
E0 = RAW.copy()
b = float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1)))
E0[:, C.wcols] = C.W.apply(E0[:, C.wcols], a.gamma, 0.0, renorm=False)
E0[:, C.wcols] *= b / float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1)))
nE = float(np.median(np.linalg.norm(E0, axis=1)))

def finish(E):
    return C.dh_exact(E * (ALPHA * m.effect_scale /
                           (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)))

def raw_of(x):
    V = np.zeros_like(E0)
    for i, p in enumerate(C.perts):
        j = m.pert_index[x].get(p)
        if j is None or j >= m.s[x].shape[0]:
            continue
        V[i, m.tcols[x]] = m.s[x][j][m.dcols[x]]
    return V * (nE / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0))

def cd4_prod():
    V = sum(raw_of(x) for x in CD4) / 3.0
    have = np.ones(len(C.perts), bool)
    for x in CD4:
        have &= np.array([(m.pert_index[x].get(p) is not None
                           and m.pert_index[x][p] < m.s[x].shape[0]) for p in C.perts])
    V[~have] = 0.0
    sc = np.abs(V).mean(0); live = sc[sc > 0]
    sc = np.maximum(sc, np.percentile(live, 25.0))
    bb = float(np.median(np.linalg.norm(V, axis=1))); V = V / sc
    V *= bb / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0)
    panel = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str))
    V[:, np.array([i for i, g in enumerate(m.target_genes) if str(g) in panel])] = 0.0
    return V * (nE / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0))

SELF = RAW * (nE / (float(np.median(np.linalg.norm(RAW, axis=1))) or 1.0))
dh0 = finish(E0); X0, A0, W0 = C.XAW(dh0)
p0, n0 = C.pds(dh0), C.nmae(dh0)
print(f"# {a.harness}, gamma={a.gamma}: base kappa={C.kappa(dh0):.5f} pds={p0:.5f} nmae={n0:.5f}")
print(f"  {'added vector':<26}" + "".join(f"{'f='+str(f):>9}" for f in FRACS))
for lbl, V in (("SELF (unwhitened base mix)", SELF), ("cd4_mean3_gn_prod", cd4_prod())):
    ms, ps, ns = [], [], []
    for f in FRACS:
        dh = finish(E0 + f * V); X, A, W = C.XAW(dh)
        ms.append(float(np.sign(X) * np.sqrt(W / W0)))
        ps.append(C.pds(dh) - p0); ns.append(C.nmae(dh) - n0)
    print(f"  {lbl:<26}" + "".join(f"{v:>9.3f}" for v in ms))
    print(f"    {'d_pds':<24}" + "".join(f"{v:>+9.5f}" for v in ps))
    print(f"    {'d_nmae':<24}" + "".join(f"{v:>+9.5f}" for v in ns))
# and the equivalent whitening gammas, for the same comparison in 097's own units
print(f"\n  {'gamma':<26}" + "".join(f"{g:>9}" for g in (0.0, 0.3, 0.45, 0.65, 0.85)))
ms, ps, ns = [], [], []
for g in (0.0, 0.3, 0.45, 0.65, 0.85):
    E = RAW.copy()
    if g > 0:
        bb = float(np.median(np.linalg.norm(E[:, C.wcols], axis=1)))
        E[:, C.wcols] = C.W.apply(E[:, C.wcols], g, 0.0, renorm=False)
        E[:, C.wcols] *= bb / float(np.median(np.linalg.norm(E[:, C.wcols], axis=1)))
    dh = finish(E); X, A, W = C.XAW(dh)
    ms.append(float(np.sign(X) * np.sqrt(W / W0)))
    ps.append(C.pds(dh) - p0); ns.append(C.nmae(dh) - n0)
print(f"  {'M':<26}" + "".join(f"{v:>9.3f}" for v in ms))
print(f"  {'d_pds':<26}" + "".join(f"{v:>+9.5f}" for v in ps))
print(f"  {'d_nmae':<26}" + "".join(f"{v:>+9.5f}" for v in ns))
