#!/usr/bin/env python3
"""635 — how much of the "context's own mean-response programme" sign prize is reachable, once the reference control
pool's shared sampling error is removed? (council 600; challenge from [7772dc]/[4d0207] to d4cd5f.md §3)

Every target's reference lfc is taken against ONE control pool (18,400 cells), so its sampling error s_g is shared by all
targets. The leave-one-target-out (LOTO) mean response used as an "oracle" in 441 / 610 contains -s_g, and on a null gene
knowing s_g alone is worth up to 0.5 + arctan(sqrt(400/18400))/pi = 0.547 sign accuracy. No predictor can know s_g.

Here, on h1-test (323's reference: the real per-cell control CPM `cs`, and each target's per-cell CPM arithmetic mean
rebuilt exactly as m_p = ref_lin * 2^lfc), the controls are split into disjoint halves A and B (seed 635):
  truth_B   = log2(m_p / m_B)                       (the evaluated target vs control half B)
  oracle_A  = LOTO mean over q != p of log2(m_q / m_A)   (clean: shares no control cells with truth_B)
  oracle_B  = LOTO mean over q != p of log2(m_q / m_B)   (contaminated: shares half B, like 441's oracle)
  CLEAN4, all-down                                        (target-independent public signs)
Scored as sign accuracy on (i) the fid filler: top-2,200-CPM block genes NOT in the reference gate, (ii) gate genes.
Analysis only (< 4 GB). Writes work/halves635.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent / "work"
WORK.mkdir(exist_ok=True)
EPS = 1e-9

z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(z["keep"], bool)
G = z["genes"].astype(str)[keep]
T = z["targets"].astype(str)
q, lfc, ref_lin = z["q"], z["lfc"], np.asarray(z["ref_lin"], float)
cs = z["cs"]  # (genes, 18,400 control cells): per-cell log1p(CPM), SORTED along cells per gene (323 diagnose_de323)
assert cs.shape[0] == len(G)
# a random set of positions per gene is a random half of that gene's control values (cell identity across genes is not
# needed for per-gene sign accuracy); means are taken on expm1 = per-cell CPM, in gene chunks
chk = np.concatenate([np.expm1(cs[i:i + 2000]).mean(1) for i in range(0, cs.shape[0], 2000)])
print(f"KA: control mean from cs vs ref_lin max rel err {np.nanmax(np.abs(chk - ref_lin) / np.maximum(ref_lin, 1e-12)):.2e}", flush=True)

gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (G[None, :] != T[:, None])
mp = np.where(np.isfinite(lfc), ref_lin[None, :] * np.exp2(lfc), np.nan)  # per-cell CPM arithmetic mean of each target
# KA: rebuilt lfc vs stored
lfc_re = np.log2((mp + EPS) / (ref_lin[None, :] + EPS))
print(f"KA: rebuilt lfc max abs err {np.nanmax(np.abs(lfc_re - lfc)):.2e}", flush=True)

rng = np.random.default_rng(635)
# cs is sorted along cells PER GENE, so one shared position split would put the same rank positions in half A for every
# gene (a fake common shift across genes). Split each gene independently: an independent random half of its values.
n = cs.shape[1]
h = n // 2
mA = np.empty(cs.shape[0])
mB = np.empty(cs.shape[0])
for i in range(0, cs.shape[0], 1000):
    X = np.expm1(cs[i:i + 1000])
    idx = np.argsort(rng.random(X.shape), axis=1)
    Xs = np.take_along_axis(X, idx, axis=1)
    mA[i:i + 1000] = Xs[:, :h].mean(1)
    mB[i:i + 1000] = Xs[:, h:].mean(1)
    del X, idx, Xs
print(f"KA: (mA + mB)/2 vs ref_lin max rel err {np.nanmax(np.abs((mA + mB) / 2 - ref_lin) / np.maximum(ref_lin, 1e-12)):.2e}", flush=True)
del cs

LA = np.log2((mp + EPS) / (mA[None, :] + EPS))
LB = np.log2((mp + EPS) / (mB[None, :] + EPS))
LF = lfc
nT = len(T)


def loto(L):
    s = np.nansum(L, 0)
    n = np.isfinite(L).sum(0)
    return (s[None, :] - np.nan_to_num(L)) / np.maximum(n[None, :] - np.isfinite(L), 1)


OA, OB, OF = loto(LA), loto(LB), loto(LF)

cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
names = list(cp["bank_names"].astype(str))
P4 = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
c = np.array([pos.get(g, -1) for g in G])
C4 = np.full(len(G), np.nan)
C4[c >= 0] = P4[c[c >= 0]]

block = np.zeros(len(G), bool)
block[np.argsort(-ref_lin)[:2200]] = True
own = G[None, :] == T[:, None]
filler = block[None, :] & ~gate & ~own & np.isfinite(lfc)
gatem = gate


def acc(pred_dn, truthL, mask):
    m = mask & np.isfinite(truthL) & (truthL != 0)
    if pred_dn.ndim == 1:
        pred_dn = np.broadcast_to(pred_dn[None, :], truthL.shape)
    ok = (pred_dn == (truthL < 0)) & m
    w = np.abs(np.nan_to_num(truthL)) * m
    return dict(n=int(m.sum()), acc=float(ok.sum() / max(m.sum(), 1)), acc_w=float((w * ok).sum() / max(w.sum(), 1e-12)))


rows = []
for truth_name, TL in (("truth_full (ref lfc)", LF), ("truth_B (control half B)", LB)):
    for mask_name, M in (("fid filler: block, non-DE", filler), ("gate genes", gatem)):
        cov4 = np.isfinite(C4)[None, :]
        cand = [("all-down", np.ones(len(G), bool), M),
                ("CLEAN4 (its coverage)", C4 > 0.5, M & cov4),
                ("all-down on CLEAN4 coverage", np.ones(len(G), bool), M & cov4),
                ("LOTO oracle vs FULL controls (441-style)", OF < 0, M),
                ("LOTO oracle vs half A (CLEAN, no shared control)", OA < 0, M),
                ("LOTO oracle vs half B (shares half B)", OB < 0, M)]
        for lab, pdn, mm in cand:
            r = acc(pdn, TL, mm)
            rows.append(dict(truth=truth_name, mask=mask_name, source=lab, **r))
print(f"\n{'truth':26s} {'mask':26s} {'source':50s} {'n':>8s} {'acc':>7s} {'acc_w':>7s}")
for r in rows:
    print(f"{r['truth']:26s} {r['mask']:26s} {r['source']:50s} {r['n']:8d} {r['acc']:7.4f} {r['acc_w']:7.4f}")
(WORK / "halves635.json").write_text(json.dumps(rows, indent=1))
