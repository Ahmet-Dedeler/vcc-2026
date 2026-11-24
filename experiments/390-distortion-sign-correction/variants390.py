#!/usr/bin/env python3
"""390g -- next variants on the DE channel, all in 390 C5b's sign-preserving mass-balanced form (H1, analytic, exact
known answers inherited). Steered set S (called genes); flip set F within S takes the prior sign at |E0|; down-signed
targets in S are shrunk by one per-target factor k in [0, 1] (floor 0.002) so S's arithmetic mass is unchanged.
Arms: base (block, CLEAN4) | head+block CLEAN4 | block, flip only |E0| below the target's block q50 / q80 |
LOLO6 / K562ESS priors | CLEAN4 flips only where |P - 0.5| >= 0.02 / 0.05."""
import numpy as np, polars as pl, json
from pathlib import Path
from scipy.optimize import brentq
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
cp = np.load(REPO / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
keep = r["keep"].astype(bool); genes = r["genes"].astype(str)[keep]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); QHAT = s["q"].astype(float); PHAT = s["p"].astype(float)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
qrm = r["q"].astype(float); lrm = r["lfc"].astype(float); ref_lin = r["ref_lin"].astype(float)
sums = s["sums"].astype(float); A0 = (sums / sums.sum(1, keepdims=True) * 1e6)[:, keep]
nT, nG = E0.shape; own = np.equal.outer(targets, genes); feat = np.tile(genes, (nT, 1))
sig = np.isfinite(QHAT) & (QHAT < .05); adjud = np.isfinite(lrm) & (lrm != 0)
callad = sig & adjud & ~own; pool = (qrm < .05) & ~own; n_conf = pool.sum(1)
yld = np.minimum(1, callad.sum(1) / np.maximum(n_conf, 1))
head = np.zeros((nT, nG), bool); block = np.zeros((nT, nG), bool)
for i in range(nT):
    head[i, np.argsort(-np.abs(E1[i]), kind="stable")[:400]] = True
    block[i, np.argsort(-np.where(head[i], -np.inf, ref_lin), kind="stable")[:2200]] = True
qm = np.where(np.isfinite(QHAT), QHAT, np.inf); pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
cur = np.empty((nT, nG), dtype=int)
for i in range(nT): cur[i] = np.lexsort((feat[i], -np.abs(E0[i]), pm[i], qm[i], ~sig[i]))
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float); cg = {g: j for j, g in enumerate(cp["genes"].astype(str))}
cmap = np.array([cg.get(g, -1) for g in genes]); cov = cmap >= 0
def prior(banks):
    P = np.full(nG, np.nan); P[cov] = PB[[names.index(b) for b in banks]][:, cmap[cov]].mean(0); return P
PRI = {"CLEAN4": prior(("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")),
       "LOLO6": prior([n for n in names if n != "h1_2025"]), "K562ESS": prior(("k562_essential",))}
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
def metrics(L):
    SG = np.where(L == 0, 1.0, np.sign(L)); agree = SG == np.sign(lrm)
    prec = np.array([agree[i][callad[i]].mean() if callad[i].any() else np.nan for i in range(nT)])
    nm = np.array([np.abs(L[k, GIDX[j]] - LREF[j]).sum() / np.abs(LREF[j]).sum() for j, k in enumerate(elig)])
    return dict(fid=prec * yld, nmae=nm, reach=walk_reach(cur, adjud & agree & pool, adjud & pool, n_conf))
def c5b(steer_mask, flip_fn, P):
    psign = np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), -1.0)
    L = E0.copy(); ks = []; flips = []
    for i in range(nT):
        S = steer_mask[i] & sig[i]; F = S & flip_fn(i, P)
        tgt = np.where(F, psign * np.abs(E0[i]), E0[i]); side = S & (tgt < 0)
        t2 = lambda k: np.where(side, np.minimum(tgt * k, -0.002), tgt)
        R = lambda k: (A0[i, S] * (2 ** (t2(k)[S] - E0[i, S]) - 1)).sum()
        k = brentq(R, 0.0, 1.0) if R(0.0) * R(1.0) < 0 else (0.0 if abs(R(0.0)) < abs(R(1.0)) else 1.0)
        L[i, S] = t2(k)[S]; ks.append(k); flips.append(float((np.sign(L[i, S]) != np.sign(E0[i, S])).mean()))
        Rres = (A0[i, S] * (2 ** (L[i, S] - E0[i, S]) - 1)).sum(); O = 1e6 - A0[i, S].sum()
        L[i, ~S] = E0[i, ~S] + np.log2((O - Rres) / O)
    return L, np.array(ks), np.array(flips)
allF = lambda i, P: np.ones(nG, bool)
def weakF(q):
    def f(i, P):
        b = block[i] & sig[i]; thr = np.quantile(np.abs(E0[i, b]), q) if b.any() else 0
        return np.abs(E0[i]) <= thr
    return f
def confF(c): return lambda i, P: ~np.isfinite(P) | (np.abs(P - 0.5) >= c)
notcalled = ~sig
ARMS = [("base_block_CLEAN4", block, allF, "CLEAN4"), ("headblock_CLEAN4", head | block, allF, "CLEAN4"),
        ("block_weak50_CLEAN4", block, weakF(0.5), "CLEAN4"), ("block_weak80_CLEAN4", block, weakF(0.8), "CLEAN4"),
        ("block_LOLO6", block, allF, "LOLO6"), ("block_K562ESS", block, allF, "K562ESS"),
        ("block_conf02_CLEAN4", block, confF(0.02), "CLEAN4"), ("block_conf05_CLEAN4", block, confF(0.05), "CLEAN4")]
import sys
if len(sys.argv) > 1 and sys.argv[1] == "order":
    # 390 §8: reach order by design |H| (production key) with a boost where the design's sign agrees with the prior.
    # Signs are 392's (C5b block CLEAN4); only the order moves, so fid / nmae / jac are unchanged by construction.
    L, _, _ = c5b(block, allF, PRI["CLEAN4"])
    SG = np.where(L == 0, 1.0, np.sign(L)); agree_ref = (SG == np.sign(lrm))
    P = PRI["CLEAN4"]; psign = np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), 0.0)
    agree = (np.sign(E1) == psign[None, :]) & (psign[None, :] != 0)
    match = adjud & agree_ref & pool
    def order_by(key):
        o = np.empty((nT, nG), dtype=int)
        for i in range(nT): o[i] = np.lexsort((feat[i], -key[i], ~sig[i]))
        return o
    base = walk_reach(order_by(np.abs(E1)), match, adjud & pool, n_conf)
    base_E0 = walk_reach(order_by(np.abs(E0)), adjud & (np.where(E0 == 0, 1.0, np.sign(E0)) == np.sign(lrm)) & pool, adjud & pool, n_conf)
    deliv = walk_reach(cur, match, adjud & pool, n_conf)
    print("reach: 392 signs, delivered order", round(float(np.nanmean(deliv)), 5), "| |H| order", round(float(np.nanmean(base)), 5),
          "| cap15 signs, |E0| order", round(float(np.nanmean(base_E0)), 5), flush=True)
    rng = np.random.default_rng(20260913); res = {"delivered_order": float(np.nanmean(deliv)), "H_order": float(np.nanmean(base))}
    for beta in (0.25, 1.0, 4.0, 1e3):
        R = walk_reach(order_by(np.abs(E1) * (1 + beta * agree)), match, adjud & pool, n_conf)
        dv = R - base; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk))); bs = dv[okk][D].mean(1)
        res[f"beta_{beta}"] = {"reach": round(float(np.nanmean(R)), 5), "d_vs_H_order": round(float(np.nanmean(dv)), 5),
                               "ci": [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]}
        print(f"beta {beta}", json.dumps(res[f"beta_{beta}"]), flush=True)
    sh = rng.permutation(nG)
    Rs = walk_reach(order_by(np.abs(E1) * (1 + 1.0 * agree[:, sh])), match, adjud & pool, n_conf)
    res["beta_1_shuffled_agreement"] = round(float(np.nanmean(Rs - base)), 5); print("shuffled agreement beta 1 d:", res["beta_1_shuffled_agreement"])
    (OUT / "variants390_order.json").write_text(json.dumps(res, indent=1)); sys.exit(0)
if len(sys.argv) > 1 and sys.argv[1] == "noncalled":
    def c5b_S(Smask, flip_fn, P):
        psign = np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), -1.0)
        L = E0.copy()
        for i in range(nT):
            S = Smask[i] & np.isfinite(P); F = S & flip_fn(i, P)
            tgt = np.where(F, psign * np.abs(E0[i]), E0[i]); side = S & (tgt < 0)
            t2 = lambda k: np.where(side, np.minimum(tgt * k, -0.002), tgt)
            R = lambda k: (A0[i, S] * (2 ** (t2(k)[S] - E0[i, S]) - 1)).sum()
            k = brentq(R, 0.0, 1.0) if R(0.0) * R(1.0) < 0 else (0.0 if abs(R(0.0)) < abs(R(1.0)) else 1.0)
            L[i, S] = t2(k)[S]
            Rres = (A0[i, S] * (2 ** (L[i, S] - E0[i, S]) - 1)).sum(); O = 1e6 - A0[i, S].sum()
            L[i, ~S] = E0[i, ~S] + np.log2((O - Rres) / O)
        return L
    nc = np.broadcast_to(ref_lin > 0, (nT, nG)) & ~sig      # universe genes we do not call (scored only by nmae / reach tail)
    M0 = metrics(E0); rng = np.random.default_rng(20260913); res = {}
    for name, (Smask, fn) in {"noncalled_CLEAN4": (nc, allF), "noncalled_conf05_CLEAN4": (nc, confF(0.05)),
                              "block_plus_noncalled_CLEAN4": ((block & sig) | nc, allF)}.items():
        L = c5b_S(Smask, fn, PRI["CLEAN4"]); M = metrics(L); row = {}
        for kk in M:
            dv = M[kk] - M0[kk]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk)))
            bs = dv[okk][D].mean(1); row[kk] = round(float(np.nanmean(dv)), 5); row[f"{kk}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
        row["overall_est"] = round(0.561 * row["fid"] - 0.255 * row["nmae"] + 0.188 * row["reach"], 5)
        res[name] = row; print(name, json.dumps(row), flush=True)
    (OUT / "variants390_noncalled.json").write_text(json.dumps(res, indent=1)); sys.exit(0)
M0 = metrics(E0); rng = np.random.default_rng(20260913); out = {}
for name, mask, fn, pn in ARMS:
    L, ks, fl = c5b(mask, fn, PRI[pn]); M = metrics(L); row = {}
    for kk in M:
        dv = M[kk] - M0[kk]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk)))
        bs = dv[okk][D].mean(1); row[kk] = round(float(np.nanmean(dv)), 5); row[f"{kk}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
    row["overall_est"] = round(0.561 * row["fid"] - 0.255 * row["nmae"] + 0.188 * row["reach"], 5)
    row["k_median"] = round(float(np.median(ks)), 3); row["flip_frac_steered_median"] = round(float(np.median(fl)), 3)
    out[name] = row; print(name, json.dumps(row), flush=True)
(OUT / "variants390.json").write_text(json.dumps(out, indent=1))
