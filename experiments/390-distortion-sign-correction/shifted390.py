#!/usr/bin/env python3
"""390f -- CLEAN4 consensus signs on the called CPM block, priced WITH the compositional shift (H1, analytic).
Every per-cell CPM sums to 1e6, so the block's arithmetic-mass change R is absorbed by all other genes: their delivered
LFC moves by sh = log2((O - R)/O), O = 1e6 - block mass (uniform-absorption model). Arms:
  C0  CLEAN4 signs at |E0|, no shift (consensus390's optimistic number, for reference)
  C1  CLEAN4 signs at |E0| + shift on every non-block gene (what an unbalanced emitter would deliver)
  C2  balanced: add one per-target constant c to the block's target LFCs so R = 0 (no shift); signs re-read after c
  C3  confident-only: flip only where |P_down - 0.5| >= 0.02, others keep delivered sign; + shift
  C4  C3 balanced like C2
fid / nmae / reach exact as signrule390; paired bootstrap vs delivered."""
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
cb = block & sig
qm = np.where(np.isfinite(QHAT), QHAT, np.inf); pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
cur = np.empty((nT, nG), dtype=int)
for i in range(nT): cur[i] = np.lexsort((feat[i], -np.abs(E0[i]), pm[i], qm[i], ~sig[i]))
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float); cg = {g: j for j, g in enumerate(cp["genes"].astype(str))}
cmap = np.array([cg.get(g, -1) for g in genes]); cov = cmap >= 0
idx = [names.index(n) for n in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]
P = np.full(nG, np.nan); P[cov] = PB[idx][:, cmap[cov]].mean(0)
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
def build(flipmask_fn, balance, shift):
    L = E0.copy(); shifts, cs = [], []
    for i in range(nT):
        b = cb[i]; fm = flipmask_fn(i) & b
        psign = np.where(P > 0.5, -1.0, 1.0)
        tgt = np.where(fm, psign * np.abs(E0[i]), E0[i])          # block target LFC (non-flipped keep delivered)
        c = 0.0
        if balance:
            f = lambda c_: (A0[i, b] * (2 ** (tgt[b] + c_ - E0[i, b]) - 1)).sum()
            c = brentq(f, -0.5, 0.5)
        L[i, b] = tgt[b] + c; cs.append(c)
        R = (A0[i, b] * (2 ** (L[i, b] - E0[i, b]) - 1)).sum(); O = 1e6 - A0[i, b].sum()
        sh = np.log2((O - R) / O) if shift else 0.0; shifts.append(sh)
        L[i, ~b] = E0[i, ~b] + sh
    return L, np.array(shifts), np.array(cs)
allcov = lambda i: np.isfinite(P) | ~np.isfinite(P)             # uncovered genes: P nan -> P>0.5 False -> up; handle below
def cover_all(i): return np.ones(nG, bool)
Pfill = np.where(np.isfinite(P), P, 1.0)                          # uncovered -> down (371 fallback)
P = Pfill
conf = lambda i: np.abs(P - 0.5) >= 0.02
M0 = metrics(E0); rng = np.random.default_rng(20260913); out = {}
def build_signkeep(mode):
    """C5: zero the block's arithmetic-mass change WITHOUT crossing zero. After CLEAN4 flips at |E0|, solve one per-target
    factor k on magnitudes: mode 'up' scales up-signed block targets (k >= 0), mode 'down' scales down-signed ones; if the
    chosen side cannot reach R = 0 inside [0, 4] the residual goes to the uniform shift. Signs never change."""
    L = E0.copy(); shifts, cs, clips = [], [], []
    psign = np.where(P > 0.5, -1.0, 1.0)
    lo, hi = -0.3528, 0.3029
    for i in range(nT):
        b = cb[i]; tgt = psign * np.abs(E0[i]); side = (tgt > 0) if mode == "up" else (tgt < 0)
        def Rk(k):
            t2 = np.where(side, tgt * k, tgt)
            return (A0[i, b] * (2 ** (t2[b] - E0[i, b]) - 1)).sum()
        if Rk(0.0) * Rk(4.0) < 0: k = brentq(Rk, 0.0, 4.0)
        else: k = 0.0 if abs(Rk(0.0)) < abs(Rk(4.0)) else 4.0
        t2 = np.where(side, tgt * k, tgt); L[i, b] = t2[b]; cs.append(k)
        mv = t2[b] - E0[i, b]; clips.append(float(((mv < lo) | (mv > hi)).mean()))
        R = (A0[i, b] * (2 ** (L[i, b] - E0[i, b]) - 1)).sum(); O = 1e6 - A0[i, b].sum()
        sh = np.log2((O - R) / O); shifts.append(sh); L[i, ~b] = E0[i, ~b] + sh
    return L, np.array(shifts), np.array(cs), np.array(clips)
ARMS = [("C0_noshift", cover_all, False, False), ("C1_shift", cover_all, False, True),
        ("C2_balanced", cover_all, True, True), ("C3_conf_shift", conf, False, True), ("C4_conf_balanced", conf, True, True),
        ("C5a_scale_up", "up", None, None), ("C5b_scale_down", "down", None, None)]
for name, fn, bal, sh in ARMS:
    if isinstance(fn, str):
        L, shifts, cs, clips = build_signkeep(fn)
    else:
        L, shifts, cs = build(fn, bal, sh); clips = np.zeros(nT)
    M = metrics(L); row = {}
    for k in M:
        dv = M[k] - M0[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (10000, len(okk)))
        bs = dv[okk][D].mean(1); row[k] = round(float(np.nanmean(dv)), 5); row[f"{k}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
    row["overall_est"] = round(0.561 * row["fid"] - 0.255 * row["nmae"] + 0.188 * row["reach"], 5)
    row["shift_median"] = round(float(np.median(shifts)), 5); row["c_or_k_q10_50_90"] = [round(float(x), 4) for x in np.quantile(cs, [.1, .5, .9])]
    row["clip_median"] = round(float(np.median(clips)), 5)
    row["block_prec"] = round(float(((np.sign(L) == np.sign(lrm)))[cb & adjud].mean()), 4)
    out[name] = row; print(name, json.dumps(row))
(OUT / "shifted390.json").write_text(json.dumps(out, indent=1))
