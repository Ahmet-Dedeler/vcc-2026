#!/usr/bin/env python3
"""390b -- sign rules for a decoupled DE channel (exact sign delivery, bulk untouched), H1, analytic.
R0 delivered signs | R1 every call down | R2 |H|-head (top 400 by |E1|) keeps R0, rest of calls down |
R3 cross-fitted rule: within each |E0| decile of the call set and each delivered sign, flip to the other
sign if the training-half precision of that (decile, sign) cell is < 0.5 (2-fold over targets, 20 splits).
Outside the call set every rule keeps E0. Magnitudes |E0| everywhere. KA inherited from price390.py."""
import numpy as np, polars as pl, json
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); QHAT = s["q"].astype(float); PHAT = s["p"].astype(float)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
qrm = r["q"].astype(float); lrm = r["lfc"].astype(float)
nT, nG = E0.shape; own = np.equal.outer(targets, genes); feat = np.tile(genes, (nT, 1))
sig = np.isfinite(QHAT) & (QHAT < .05); adjud = np.isfinite(lrm) & (lrm != 0)
callad = sig & adjud & ~own; pool = (qrm < .05) & ~own; n_conf = pool.sum(1)
yld = np.minimum(1, callad.sum(1) / np.maximum(n_conf, 1))
head = np.zeros_like(sig)
for i in range(nT): head[i, np.argsort(-np.abs(E1[i]))[:400]] = True
qm = np.where(np.isfinite(QHAT), QHAT, np.inf); pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
cur = np.empty((nT, nG), dtype=int)
for i in range(nT): cur[i] = np.lexsort((feat[i], -np.abs(E0[i]), pm[i], qm[i], ~sig[i]))
refdown = lrm < 0
base = {"P_ref_down_calls": float(refdown[callad].mean()),
        "P_ref_down_head_calls": float(refdown[callad & head].mean()), "P_ref_down_block_calls": float(refdown[callad & ~head].mean()),
        "share_calls_in_head": float((callad & head).sum() / callad.sum()),
        "prec_E0_head": float(((np.sign(E0) == np.sign(lrm)))[callad & head].mean()), "prec_E0_block": float(((np.sign(E0) == np.sign(lrm)))[callad & ~head].mean()),
        "P_E0_down_calls": float((E0 < 0)[callad].mean()),
        "prec_E0_up_calls": float((lrm > 0)[callad & (E0 > 0)].mean()), "prec_E0_down_calls": float((lrm < 0)[callad & (E0 < 0)].mean()),
        "P_ref_down_gate": float(refdown[pool & adjud].mean())}
print("BASE", json.dumps({k: round(v, 4) for k, v in base.items()}))
SG0 = np.where(E0 == 0, 1.0, np.sign(E0))
R1 = np.where(sig, -1.0, SG0); R2 = np.where(sig & ~head, -1.0, SG0)
# R3 cross-fitted
dec = np.zeros_like(E0, dtype=int)
for i in range(nT):
    c = sig[i]
    if c.sum() > 10:
        qs = np.quantile(np.abs(E0[i][c]), np.linspace(0, 1, 11)[1:-1]); dec[i] = np.searchsorted(qs, np.abs(E0[i]))
agree0 = SG0 == np.sign(lrm)
R3s = []
for split in range(20):
    rng = np.random.default_rng(390 + split); half = rng.permutation(nT) < nT // 2
    R3 = SG0.copy()
    for trn, tst in ((half, ~half), (~half, half)):
        for d in range(10):
            for sg in (-1.0, 1.0):
                cell = callad & (dec == d) & (SG0 == sg)
                tr = cell & trn[:, None]
                if tr.sum() >= 30 and agree0[tr].mean() < 0.5:
                    R3[(sig & (dec == d) & (SG0 == sg)) & tst[:, None]] = -sg
    R3s.append(R3)
ref = pl.read_parquet(W / "reference_de.parquet"); tcol = ref["target"].to_numpy().astype(str); gcol = ref["feature"].to_numpy().astype(str)
gi = {g: i for i, g in enumerate(genes)}; qa = ref["p_adj"].to_numpy().astype(float); la = ref["log2_fold_change"].to_numpy().astype(float)
elig, GIDX, LREF = [], [], []
for k, t in enumerate(targets):
    m = tcol == t; g = gcol[m]; gate = (qa[m] < .05) & np.isfinite(la[m]) & (la[m] != 0) & (g != t)
    if gate.sum() < 10: continue
    elig.append(k); GIDX.append(np.array([gi[x] for x in g[gate]])); LREF.append(la[m][gate])
def walk_reach(order_idx, match, adj, nconf):
    R = np.full(len(match), np.nan)
    for i in range(len(match)):
        if nconf[i] == 0: continue
        ma = match[i][order_idx[i]][adj[i][order_idx[i]]]
        if ma.size == 0: R[i] = 0.0; continue
        okk = np.flatnonzero(np.cumsum(ma) / np.arange(1, ma.size + 1) >= 0.9)
        R[i] = (okk.max() + 1 if okk.size else 0) / nconf[i]
    return R
def metrics(SG):
    agree = SG == np.sign(lrm)
    prec = np.array([agree[i][callad[i]].mean() if callad[i].any() else np.nan for i in range(nT)])
    nm = np.array([np.abs(SG[k, GIDX[j]] * np.abs(E0[k, GIDX[j]]) - LREF[j]).sum() / np.abs(LREF[j]).sum() for j, k in enumerate(elig)])
    R = walk_reach(cur, adjud & agree & pool, adjud & pool, n_conf)
    return dict(fid=prec * yld, nmae=nm, reach=R)
M0 = metrics(SG0); rng = np.random.default_rng(20260913)
def summ(M, name):
    row = {k: float(np.nanmean(v)) for k, v in M.items()}
    for k in M:
        dv = M[k] - M0[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk)))
        bs = dv[okk][D].mean(1); row[f"d_{k}"] = float(np.nanmean(dv)); row[f"d_{k}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
    print(name, json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()})); return row
out = {"base": base, "R0": summ(M0, "R0"), "R1": summ(metrics(R1), "R1_all_down"), "R2": summ(metrics(R2), "R2_block_down")}
m3 = [metrics(R) for R in R3s]; M3 = {k: np.nanmean([m[k] for m in m3], axis=0) for k in M0}
out["R3"] = summ(M3, "R3_crossfit_mean20"); out["R3_flip_frac_calls"] = float(np.mean([(R != SG0)[sig].mean() for R in R3s]))
print("R3 flip fraction of calls", round(out["R3_flip_frac_calls"], 4))
(OUT / "signrule390.json").write_text(json.dumps(out, indent=1))
