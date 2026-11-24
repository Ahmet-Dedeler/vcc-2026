"""436 (d) side question, held out from H1: which table set should the October donor sign average?

`mode="all"` (the mandate's rule: every covering donor table, i.e. the essential screens plus the three production donors)
against `mode="ess"` (the essential screens only). H1 prefers `ess` by -0.0038 +- 0.0027 nmae on its 20 covered targets
(ana436.json), which is 1.4 se and read on the instrument itself. This asks the same question on truth that is not H1.

Truth: one essential screen L at a time (its own delta row, pert minus control, log1p-CP50k). Predictors use only donors
from a DIFFERENT study than L (the board's situation: no donor shares Arc's study). Studies: Replogle = K562-essential,
RPE1, K562-GWPS; Nadig = Jurkat, HepG2; Orion = HCT116, HEK293T.
  E = sign of the mean over the cross-study ESSENTIAL screens covering the target;
  A = sign of the mean over the cross-study essential screens AND cross-study production donors covering the target.
Targets: in L and in >= 1 cross-study essential screen. Genes: L's z1000 (top 1000 by |delta_L| * sqrt(control CP50k),
430's best gate proxy, Spearman 0.94 with H1's real gate), target gene excluded. Scored where the predictor is defined;
coverage reported. Split: genes E defines (both measure) vs genes only a production donor measures (A defines, E not).
Writes work/xcheck436.json."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import signs436 as S6

HERE = Path(__file__).resolve().parent
STUDY = {"k562_essential": "replogle", "rpe1": "replogle", "k562_gwps": "replogle",
         "nadig_jurkat": "nadig", "nadig_hepg2": "nadig", "orion_hct116": "orion", "orion_hek293t": "orion"}
DIRS = {**S6.ESS_DIRS, **S6.PANEL_DIRS}
K = 1000


def load(nm):
    d = DIRS[nm]
    P = S6._lines(d / "perts.txt"); G = S6._lines(d / "genes.txt")
    M = np.load(d / "delta.npy", mmap_mode="r")
    b = np.load(d / "b_ctrl.npy")
    assert M.shape == (len(P), len(G)) and b.shape[-1] == len(G), (nm, M.shape, b.shape)
    return {p: i for i, p in enumerate(P)}, {g: i for i, g in enumerate(G)}, M, np.asarray(b, dtype=np.float64).ravel()


T = {nm: load(nm) for nm in DIRS}
ESS_LIST = set(S6._lines(S6.ESS_DIRS["k562_essential"] / "perts.txt"))
out = {"K": K, "folds": {}}
for L in S6.ESS_DIRS:
    pL, gL, ML, bL = T[L]
    genesL = list(gL)
    ess_x = [nm for nm in S6.ESS_DIRS if STUDY[nm] != STUDY[L]]
    pan_x = [nm for nm in S6.PANEL_DIRS if STUDY[nm] != STUDY[L]]
    tg = [t for t in pL if any(t in T[nm][0] for nm in ess_x)]
    colmaps = {nm: np.array([T[nm][1].get(g, -1) for g in genesL]) for nm in ess_x + pan_x}
    zw = np.sqrt(np.maximum(np.expm1(bL), 0.0))
    accE, accA, accA_on_E, accP_only, covE, covA, nP, onl = [], [], [], [], [], [], [], []
    for t in tg:
        y = np.asarray(ML[pL[t]], dtype=np.float64)
        sc = np.abs(y) * zw
        if t in gL:
            sc[gL[t]] = -1.0
        top = np.argsort(-sc, kind="stable")[:K]
        top = top[(sc[top] > 0) & (y[top] != 0)]

        def stack(names):
            R = np.full((len(names), top.size), np.nan)
            for i, nm in enumerate(names):
                if t in T[nm][0]:
                    c = colmaps[nm][top]; ok = c >= 0
                    R[i, ok] = np.asarray(T[nm][2][T[nm][0][t]], dtype=np.float64)[c[ok]]
            return R
        RE, RP = stack(ess_x), stack(pan_x)
        with np.errstate(all="ignore"):
            mE = np.nanmean(RE, 0); mA = np.nanmean(np.vstack([RE, RP]), 0)
        sE, sA = np.sign(np.nan_to_num(mE)), np.sign(np.nan_to_num(mA))
        yt = np.sign(y[top])
        dE, dA = sE != 0, sA != 0
        if dE.sum() < 20:
            continue
        accE.append(float(np.mean(sE[dE] == yt[dE]))); covE.append(float(dE.mean())); onl.append(t in ESS_LIST)
        accA.append(float(np.mean(sA[dA] == yt[dA]))); covA.append(float(dA.mean()))
        accA_on_E.append(float(np.mean(sA[dE] == yt[dE])))
        po = dA & ~dE
        nP.append(int(po.sum()))
        if po.sum() >= 20:
            accP_only.append(float(np.mean(sA[po] == yt[po])))
    d = np.array(accA_on_E) - np.array(accE)
    rec = {"cross_study_ess": ess_x, "cross_study_panel": pan_x, "n_targets": len(accE),
           "acc_E": float(np.mean(accE)), "cov_E": float(np.mean(covE)),
           "acc_A": float(np.mean(accA)), "cov_A": float(np.mean(covA)),
           "acc_A_on_E_genes": float(np.mean(accA_on_E)),
           "A_minus_E_on_E_genes": float(d.mean()), "A_minus_E_se": float(d.std(ddof=1) / np.sqrt(d.size)),
           "A_better_frac": float(np.mean(d > 0)),
           "panel_only_genes_median": float(np.median(nP)) if nP else 0.0,
           "acc_panel_only_genes": float(np.mean(accP_only)) if accP_only else None, "n_targets_panel_only": len(accP_only)}
    onl_ = np.array(onl, dtype=bool)
    for nm_, m_ in (("on_list", onl_), ("off_list", ~onl_)):
        rec[nm_] = {"n": int(m_.sum()), "acc_E": float(np.mean(np.array(accE)[m_])) if m_.any() else None,
                    "acc_A": float(np.mean(np.array(accA)[m_])) if m_.any() else None}
    out["folds"][L] = rec
    print(f"    on-list {rec['on_list']}  off-list {rec['off_list']}", flush=True)
    print(f"L={L:15s} n {rec['n_targets']:4d} | E acc {rec['acc_E']:.4f} (cov {rec['cov_E']:.2f}) | A acc on E's genes "
          f"{rec['acc_A_on_E_genes']:.4f}: A-E {rec['A_minus_E_on_E_genes']:+.4f} +- {rec['A_minus_E_se']:.4f} "
          f"(A better {rec['A_better_frac']:.2f}) | A acc all {rec['acc_A']:.4f} (cov {rec['cov_A']:.2f}) | "
          f"panel-only genes median {rec['panel_only_genes_median']:.0f}, acc {rec['acc_panel_only_genes']}", flush=True)
(HERE / "work/xcheck436.json").write_text(json.dumps(out, indent=1))
