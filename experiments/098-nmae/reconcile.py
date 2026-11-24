"""Two checks the whitening reading depends on.

(1) `pds` -- my analytic rank-cosine against `cell_eval2.metrics.discrimination`'s own function.
(2) The amplitude convention.  097's whitening gate re-solved each gamma onto its own rho cliff;
    098 holds `A` at a12's.  On h1 the two conventions report OPPOSITE signs for `nmae` at
    gamma = 0, so which convention is in force has to be stated with the number.
"""
from __future__ import annotations
import hashlib, json, sys
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

DEPTHS = (0.05, 0.35, 0.60)
for WHICH in ("h1", "orion"):
    OUTDIR = {"h1": f"{ROOT}/experiments/095-clean-harness/out-h1test",
              "orion": f"{ROOT}/experiments/095-clean-harness/out"}[WHICH]
    H = K.Harness(WHICH)
    z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{WHICH}.npz", allow_pickle=True)
    QQ, RL = z["q"], z["lfc"].astype(np.float64)
    # ---- (1) pds cross-check ----
    try:
        from cell_eval2.metrics.discrimination import discrimination_score
        import inspect
        print(f"[{WHICH}] discrimination_score signature: "
              f"{inspect.signature(discrimination_score)}")
    except Exception as e:
        print(f"[{WHICH}] could not import discrimination_score: {e}")
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
            b0 = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
            E[:, wcols] = W.apply(E[:, wcols], gam, 0.0, renorm=False)
            E[:, wcols] *= b0 / (float(np.median(np.linalg.norm(E[:, wcols], axis=1))) or 1.0)
        E *= 1.08 * 6.817 / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
        dh = np.vstack([np.log1p(TS * model.profile_from_effect(E[i], H.bctrl, min_cpm=5.0,
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
    print(f"\n===== {WHICH}: amplitude convention =====")
    print(f"  {'gamma':>6}{'convention':>14}{'s':>8}{'A/A_a12':>9}{'pds':>8}"
          + "".join(f"{'nmae@'+str(a):>11}" for a in DEPTHS))
    for gam in (0.0, 0.30, 0.45, 0.65, 0.85):
        raw = dh_at(gam)
        A1, _ = K.AXR(H, raw); B1 = float(_across_pert_budget(raw, H.cols, H.w))
        for conv, s in (("A held", float(np.sqrt(A0 / A1))),
                        ("rho cliff", float(np.sqrt(H.claim / B1)))):
            dh = raw * s
            A, _ = K.AXR(H, dh)
            L = H.lfc_hat(dh[H.rows])
            v = [float(np.mean([np.abs(L[k][gs] - RL[j][gs]).sum() / np.abs(RL[j][gs]).sum()
                                for k, j, gs in GR[a]])) for a in DEPTHS]
            print(f"  {gam:>6.2f}{conv:>14}{s:>8.4f}{A/A0:>9.4f}{K.pds(H, dh):>8.4f}"
                  + "".join(f"{x:>11.4f}" for x in v))
