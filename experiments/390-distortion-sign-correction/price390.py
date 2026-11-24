#!/usr/bin/env python3
"""390 -- price a distortion-free DE channel on H1 (analytic, saved arrays only).

Mechanism (PI, 2026-09-13): the welded emitter makes each predicted cell's library equal, so the
prediction's arithmetic-mean CPM == its bulk CPM, while the reference's controls have arithmetic != bulk.
Delivered log2FC = bulk LFC + d_g, d_g = log2(bulk_ctrl / arith_ctrl) (verified: corr 1.00, slope 1.00).
Estimators (signed 85 x 10779):
  E0 delivered lfc | E7 prediction bulk LFC (= E0 - d_g) | E8 sign(E1 designed) x |E0|
  E9 E7 with E0's sign kept on the first 50 genes of the delivered p-order (reach front protected)
Metrics: |l|-weighted gate sign accuracy (343), nmae with E0 magnitudes, fid-like directional precision
on called & adjudicable & not-own (KA vs 323's recorded fid), reach in the delivered order (348 walk).
"""
import numpy as np, polars as pl, json
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
h = np.load(OUT / "h1_distortion.npz")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); QHAT = s["q"].astype(float); PHAT = s["p"].astype(float)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
BL = h["bulk_lfc"]; E7 = np.where(np.isfinite(BL), BL, E0)
E8 = np.where(E1 != 0, np.sign(E1) * np.abs(E0), E0)
qrm = r["q"].astype(float); lrm = r["lfc"].astype(float)
nT, nG = E0.shape; own = np.equal.outer(targets, genes); feat = np.tile(genes, (nT, 1))
sig = np.isfinite(QHAT) & (QHAT < .05)
qm = np.where(np.isfinite(QHAT), QHAT, np.inf); pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
cur = np.empty((nT, nG), dtype=int); rank = np.empty((nT, nG), dtype=int)
for i in range(nT):
    cur[i] = np.lexsort((feat[i], -np.abs(E0[i]), pm[i], qm[i], ~sig[i])); rank[i, cur[i]] = np.arange(nG)
E9 = np.where(rank < 50, E0, E7)
EST = {"E0": E0, "E7": E7, "E8": E8, "E9": E9}
# gate (343)
ref = pl.read_parquet(W / "reference_de.parquet"); tcol = ref["target"].to_numpy().astype(str); gcol = ref["feature"].to_numpy().astype(str)
gi = {g: i for i, g in enumerate(genes)}; qa = ref["p_adj"].to_numpy().astype(float); la = ref["log2_fold_change"].to_numpy().astype(float)
elig, GIDX, LREF = [], [], []
for k, t in enumerate(targets):
    m = tcol == t; g = gcol[m]; q = qa[m]; l = la[m]
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10: continue
    elig.append(k); GIDX.append(np.array([gi[x] for x in g[gate]])); LREF.append(l[gate])
elig = np.array(elig)
adjud = np.isfinite(lrm) & (lrm != 0); pool = (qrm < .05) & ~own; n_conf = pool.sum(1); scored = n_conf > 0
def walk_reach(order_idx, match, adj, nconf):
    R = np.full(len(match), np.nan)
    for i in range(len(match)):
        if nconf[i] == 0: continue
        ma = match[i][order_idx[i]][adj[i][order_idx[i]]]
        if ma.size == 0: R[i] = 0.0; continue
        okk = np.flatnonzero(np.cumsum(ma) / np.arange(1, ma.size + 1) >= 0.9)
        R[i] = (okk.max() + 1 if okk.size else 0) / nconf[i]
    return R
res = {}
for name, X in EST.items():
    pw, nm, prec = [], [], np.full(nT, np.nan)
    for j, k in enumerate(elig):
        idx, l = GIDX[j], LREF[j]; D = np.abs(l).sum(); x = X[k, idx]
        pw.append((np.abs(l) * (np.sign(x) == np.sign(l))).sum() / D)
        nm.append(np.abs(np.sign(x) * np.abs(E0[k, idx]) - l).sum() / D)
    callad = sig & adjud & ~own
    agree = (np.sign(X) == np.sign(lrm)) & (X != 0)
    for i in range(nT):
        c = callad[i]
        if c.any(): prec[i] = agree[i][c].mean()
    match = adjud & (X != 0) & (np.sign(X) == np.sign(lrm))
    R = walk_reach(cur, match & pool, adjud & pool, n_conf)
    res[name] = dict(pw=np.array(pw), nm=np.array(nm), prec=prec, reach=R)
# KA: fid = precision x yield?  try yield = min(1, n_called_adj / n_real_pool)
npred = (sig & adjud & ~own).sum(1); yld = np.minimum(1, npred / np.maximum(n_conf, 1))
ka = {"prec_mean": float(np.nanmean(res["E0"]["prec"])), "prec_x_yield_mean": float(np.nanmean(res["E0"]["prec"] * yld)),
      "recorded_fid": 0.49755452, "recorded_nmae": 0.95083646, "nmae_E0": float(res["E0"]["nm"].mean()),
      "recorded_reach": 0.25434726, "reach_E0": float(np.nanmean(res["E0"]["reach"]))}
print("KA", json.dumps(ka))
rng = np.random.default_rng(20260913); D58 = rng.integers(0, len(elig), (10000, len(elig)))
idx85 = np.flatnonzero(np.isfinite(res["E0"]["prec"])); D85 = rng.integers(0, len(idx85), (10000, len(idx85)))
sc = np.flatnonzero(scored); DSC = rng.integers(0, len(sc), (10000, len(sc)))
out = {"ka": ka, "arms": {}}
for name in EST:
    v = res[name]; row = {"pw_gate": float(v["pw"].mean()), "nmae_e0mag": float(v["nm"].mean()),
                          "prec": float(np.nanmean(v["prec"])), "reach": float(np.nanmean(v["reach"]))}
    if name != "E0":
        b = res["E0"]
        for key, dv, DR in (("pw_gate", v["pw"] - b["pw"], D58), ("nmae_e0mag", b["nm"] - v["nm"], D58),
                            ("prec", (v["prec"] - b["prec"])[idx85], D85), ("reach", (v["reach"] - b["reach"])[sc], DSC)):
            bs = dv[DR].mean(1); row[f"d_{key}"] = float(dv.mean()); row[f"d_{key}_ci"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
    out["arms"][name] = row
    print(name, json.dumps(row))
(OUT / "price390.json").write_text(json.dumps(out, indent=1))
