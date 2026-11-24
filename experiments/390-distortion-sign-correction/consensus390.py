#!/usr/bin/env python3
"""390d -- per-gene cross-line sign prior on the CPM block, exact delivery, H1 analytic (scorer reference truth).
Head = top-400 |design|; block = top-2200 ref_lin among non-head (361); only CALLED block genes change sign.
Priors from 371 consensus_prior.npz P_banks (7 x 6346): LOLO6 = mean over banks except h1_2025 (leaks h1-test);
K562ESS = k562_essential alone; CLEAN4 = k562_gwps, k562_essential, nadig_jurkat, nadig_hepg2 (no debiased Orion).
Rule: covered block gene -> down if P_down > 0.5 else up; uncovered block gene -> down. Magnitudes |E0| or eps=0.002."""
import numpy as np, polars as pl, json
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
cp = np.load(REPO / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); QHAT = s["q"].astype(float); PHAT = s["p"].astype(float)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
qrm = r["q"].astype(float); lrm = r["lfc"].astype(float); ref_lin = r["ref_lin"].astype(float)
nT, nG = E0.shape; own = np.equal.outer(targets, genes); feat = np.tile(genes, (nT, 1))
sig = np.isfinite(QHAT) & (QHAT < .05); adjud = np.isfinite(lrm) & (lrm != 0)
callad = sig & adjud & ~own; pool = (qrm < .05) & ~own; n_conf = pool.sum(1)
yld = np.minimum(1, callad.sum(1) / np.maximum(n_conf, 1))
head = np.zeros((nT, nG), bool); block = np.zeros((nT, nG), bool)
for i in range(nT):
    head[i, np.argsort(-np.abs(E1[i]), kind="stable")[:400]] = True
    block[i, np.argsort(-np.where(head[i], -np.inf, ref_lin), kind="stable")[:2200]] = True
cb = block & sig
qm = np.where(np.isfinite(QHAT), QHAT, np.inf); pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
cur = np.empty((nT, nG), dtype=int)
for i in range(nT): cur[i] = np.lexsort((feat[i], -np.abs(E0[i]), pm[i], qm[i], ~sig[i]))
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float); cg = {g: j for j, g in enumerate(cp["genes"].astype(str))}
cmap = np.array([cg.get(g, -1) for g in genes]); cov = cmap >= 0
def prior(idx):
    P = np.full(nG, np.nan); P[cov] = PB[idx][:, cmap[cov]].mean(0); return P
PRI = {"LOLO6": prior([k for k, n in enumerate(names) if n != "h1_2025"]),
       "K562ESS": prior([names.index("k562_essential")]),
       "CLEAN4": prior([names.index(n) for n in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")])}
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
SG0 = np.where(E0 == 0, 1.0, np.sign(E0))
def metrics(SG, MAG):
    agree = SG == np.sign(lrm); X = SG * MAG
    prec = np.array([agree[i][callad[i]].mean() if callad[i].any() else np.nan for i in range(nT)])
    nm = np.array([np.abs(X[k, GIDX[j]] - LREF[j]).sum() / np.abs(LREF[j]).sum() for j, k in enumerate(elig)])
    return dict(fid=prec * yld, nmae=nm, reach=walk_reach(cur, adjud & agree & pool, adjud & pool, n_conf))
M0 = metrics(SG0, np.abs(E0)); rng = np.random.default_rng(20260913)
out = {"coverage_called_block": float(cov[None, :].repeat(nT, 0)[cb].mean()),
       "block_precision_E0": float((SG0 == np.sign(lrm))[cb & adjud].mean()),
       "block_P_ref_down": float((lrm < 0)[cb & adjud].mean()), "arms": {}}
print("BASE", json.dumps({k: round(v, 4) for k, v in out.items() if k != "arms"}))
EPSMAG = np.where(cb, 0.002, np.abs(E0))
arms = {"R2_down_absE0": (np.where(cb, -1.0, SG0), np.abs(E0)), "R2_down_eps": (np.where(cb, -1.0, SG0), EPSMAG)}
for pn, P in PRI.items():
    sg = np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), -1.0)
    SG = np.where(cb, sg[None, :], SG0)
    out[f"block_precision_{pn}"] = float((SG == np.sign(lrm))[cb & adjud].mean())
    arms[f"{pn}_absE0"] = (SG, np.abs(E0)); arms[f"{pn}_eps"] = (SG, EPSMAG)
for name, (SG, MAG) in arms.items():
    M = metrics(SG, MAG); row = {}
    for k in M:
        dv = M[k] - M0[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk)))
        bs = dv[okk][D].mean(1); row[k] = round(float(np.nanmean(dv)), 5); row[f"{k}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
    row["overall_est"] = round(0.561 * row["fid"] - 0.255 * row["nmae"] + 0.188 * row["reach"], 5)
    out["arms"][name] = row; print(name, json.dumps(row))
print("block precisions", json.dumps({k: round(v, 4) for k, v in out.items() if k.startswith("block_precision")}))
(OUT / "consensus390.json").write_text(json.dumps(out, indent=1))
