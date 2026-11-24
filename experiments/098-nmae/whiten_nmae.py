"""Re-price the WHITENING operator on `nmae`, at CHALLENGE-LIKE gate depth.

097's whitening gate priced `nmae` at BH 0.05 and read it as "disagrees in sign, ~zero".
`gateladder.py` shows the BH-0.05 cohort ranks amplitude arms with the OPPOSITE sign to the
board, so that reading was taken on the wrong cohort.  `nmae` is, in our regime, a DIRECTION
member (NMAE = 1 - <sign(lfc_ref), lfc_hat>_w), and the whitener is the project's largest
measured directional lever -- so this is the one member whose whitening verdict could flip.

Build path is `whiten_gate.py`'s verbatim, including its gamma=0.65 reconstruction check.

    GAMMAS=0,0.15,0.30,0.45,0.65,0.85 uv run python whiten_nmae.py <h1|orion>
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
GAMMAS = [float(x) for x in os.environ.get("GAMMAS", "0,0.15,0.30,0.45,0.65,0.85").split(",")]
DEPTHS = (0.05, 0.15, 0.35, 0.60)
ALPHA_V, EFFSCALE, MIN_CPM = 1.08, 6.817, 5.0

H = K.Harness(WHICH)
z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{WHICH}.npz", allow_pickle=True)
QQ, RL = z["q"], z["lfc"].astype(np.float64)
donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in ("k562_gwps", "orion_hek293t")}
model = MultiDonorModel(donors, H.genes, gene_norm=False, eb=False)
shared = set()
for n in model.names:
    shared |= set(model.target_genes[model.tcols[n]])
chal = set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str))
real_t = set(str(x) for x in H.rperts)
shared -= (chal if real_t <= chal else real_t)
ti = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
meta = json.load(open(f"{OUTDIR}/whitener_union.json"))
sha = hashlib.sha256(wcols.tobytes()).hexdigest()
assert sha == meta["wcols_sha256"], "wcols does not reproduce the cached operator"
print(f"# {WHICH}: wcols {wcols.size} sha MATCH", flush=True)
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


def gate_rows(alpha):
    out = []
    for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
        s = (QQ[j] < alpha) & np.isfinite(RL[j]) & (RL[j] != 0)
        gs = np.flatnonzero(s)
        if gs.size >= 10: out.append((k, j, gs))
    return out


GR = {a: gate_rows(a) for a in DEPTHS}
print(f"# mean |S_p| per depth: " + "  ".join(
    f"BH{a}={np.mean([g.size for _,_,g in GR[a]]):.0f}" for a in DEPTHS), flush=True)

# reconstruction check at 0.65
dh65 = dh_at(0.65)
cs = np.einsum("ij,ij->i", dh65, H.dh0) / np.maximum(
    np.linalg.norm(dh65, axis=1) * np.linalg.norm(H.dh0, axis=1), 1e-30)
print(f"# reconstruction at gamma=0.65: per-row cosine vs shipped dh median {np.median(cs):.5f} "
      f"min {cs.min():.5f}", flush=True)

base = None
print(f"\n  {'gamma':>6}{'kappa':>9}{'A':>10}{'pds':>8}"
      + "".join(f"{'nmae@BH'+str(a):>13}" for a in DEPTHS)
      + "".join(f"{'Δ@'+str(a):>10}" for a in DEPTHS))
rows = []
for gam in GAMMAS:
    dh = dh_at(gam)
    # put every gamma on the SAME amplitude as the shipped a12 (A held), so this is a pure
    # direction comparison and not an amplitude one
    A1, _ = K.AXR(H, dh)
    A0, _ = K.AXR(H, H.dh0)
    dh = dh * float(np.sqrt(A0 / A1))
    A, X = K.AXR(H, dh); B = float(_across_pert_budget(dh, H.cols, H.w))
    L = H.lfc_hat(dh[H.rows])
    v = {}
    for a in DEPTHS:
        v[a] = float(np.mean([np.abs(L[k][gs] - RL[j][gs]).sum() / np.abs(RL[j][gs]).sum()
                              for k, j, gs in GR[a]]))
    kap = X / np.sqrt(B * 1.0)
    p = K.pds(H, dh)
    if base is None or gam == 0.65: pass
    rows.append((gam, kap, A, p, v))
b = [r for r in rows if r[0] == 0.65]
b = b[0] if b else rows[-1]
for gam, kap, A, p, v in rows:
    print(f"  {gam:>6.2f}{kap:>9.5f}{A:>10.5f}{p:>8.4f}"
          + "".join(f"{v[a]:>13.4f}" for a in DEPTHS)
          + "".join(f"{v[a]-b[4][a]:>+10.4f}" for a in DEPTHS))
print(f"  (deltas are against gamma = 0.65, production; pds at 0.65 = {b[3]:.4f}, "
      f"kappa {b[1]:.5f})")
