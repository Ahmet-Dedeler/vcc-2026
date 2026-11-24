"""Can the whitener's `nmae` gain be bought without its `pds` cost?

`nmae` reads only genes with control CPM > 5 that the reference calls DE; `pds` reads the whole
gene vector as a cosine.  If the two live on separable gene sets, a blend that de-whitens only
where `nmae` looks and keeps production's operator elsewhere would be a Pareto move.  This tests
that directly:  dh = dh(0.65) + t * (dh(gam) - dh(0.65)) on a CPM-selected gene set only.

Amplitude held at the shipped a12 `A` throughout, so every row is a pure direction comparison.
"""
from __future__ import annotations
import hashlib, json, os, sys
import numpy as np
import pandas as pd
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from msekit import TS
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.whiten import Whitener

WHICH = sys.argv[1]
OUTDIR = {"h1": f"{ROOT}/experiments/095-clean-harness/out-h1test",
          "orion": f"{ROOT}/experiments/095-clean-harness/out"}[WHICH]
DEPTHS = (0.05, 0.35, 0.60)
ALPHA_V, EFFSCALE, MIN_CPM = 1.08, 6.817, 5.0
H = K.Harness(WHICH)
z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{WHICH}.npz", allow_pickle=True)
QQ, RL = z["q"], z["lfc"].astype(np.float64)
donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in ("k562_gwps", "orion_hek293t")}
model = MultiDonorModel(donors, H.genes, gene_norm=False, eb=False)
shared = set()
for n in model.names: shared |= set(model.target_genes[model.tcols[n]])
chal = set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str))
real_t = set(str(x) for x in H.rperts)
shared -= (chal if real_t <= chal else real_t)
ti = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
assert hashlib.sha256(wcols.tobytes()).hexdigest() == \
       json.load(open(f"{OUTDIR}/whitener_union.json"))["wcols_sha256"]
E0 = np.vstack([model.effect(str(x), 1.0, 0.0, weights={n: 1.0 for n in model.names})
                for x in H.tg]).astype(np.float64)
W = Whitener.load(f"{OUTDIR}/whitener_union.npz")


def dh_at(gam):
    E = E0.copy()
    if gam > 0:
        before = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
        E[:, wcols] = W.apply(E[:, wcols], gam, 0.0, renorm=False)
        E[:, wcols] *= before / (float(np.median(np.linalg.norm(E[:, wcols], axis=1))) or 1.0)
    E *= ALPHA_V * EFFSCALE / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
    dh = np.vstack([np.log1p(TS * model.profile_from_effect(E[i], H.bctrl, min_cpm=MIN_CPM,
                                                            gate_unexpressed=True)) - H.bctrl
                    for i in range(H.P)])
    for k, c in H.cols.items(): dh[k, c] = 0.0
    return dh


GR = {}
for a in DEPTHS:
    out = []
    for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
        s = (QQ[j] < a) & np.isfinite(RL[j]) & (RL[j] != 0)
        gs = np.flatnonzero(s)
        if gs.size >= 10: out.append((k, j, gs))
    GR[a] = out
A0, _ = K.AXR(H, H.dh0)


def score(dh, label):
    A1, _ = K.AXR(H, dh); dh = dh * float(np.sqrt(A0 / A1))
    A, X = K.AXR(H, dh); B = float(_across_pert_budget(dh, H.cols, H.w))
    L = H.lfc_hat(dh[H.rows])
    v = [float(np.mean([np.abs(L[k][gs] - RL[j][gs]).sum() / np.abs(RL[j][gs]).sum()
                        for k, j, gs in GR[a]])) for a in DEPTHS]
    return dict(lab=label, kap=X / np.sqrt(B), pds=K.pds(H, dh), n=v)


d65, d00 = dh_at(0.65), dh_at(0.0)
base = score(d65, "gamma 0.65 (production)")
cpm = np.zeros(H.G); cpm[H.gate] = H.ref_lin
rows = [base]
for T in (0.0, 5.0, 20.0, 50.0, 200.0):
    m = cpm > T
    for t in (0.5, 1.0):
        dh = d65.copy(); dh[:, m] = d65[:, m] + t * (d00[:, m] - d65[:, m])
        rows.append(score(dh, f"de-whiten t={t} on CPM>{T:g} ({int(m.sum())}g)"))
# the reverse: de-whiten only the LOW-expression genes
for T in (20.0, 200.0):
    m = (cpm > 0) & (cpm <= T)
    dh = d65.copy(); dh[:, m] = d00[:, m]
    rows.append(score(dh, f"de-whiten only 5<CPM<={T:g} ({int(m.sum())}g)"))
print(f"\n===== {WHICH}: blends, amplitude held at a12 =====")
print(f"  {'arm':<40}{'kappa':>9}{'pds':>8}{'Δpds':>9}"
      + "".join(f"{'Δnmae@'+str(a):>12}" for a in DEPTHS))
for r in rows:
    print(f"  {r['lab']:<40}{r['kap']:>9.5f}{r['pds']:>8.4f}{r['pds']-base['pds']:>+9.4f}"
          + "".join(f"{r['n'][i]-base['n'][i]:>+12.4f}" for i in range(len(DEPTHS))))
