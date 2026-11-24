#!/usr/bin/env python3
"""397 -- neighbour-conditioned block sign prior, H1 analytic (frame copied from 390/consensus390.py; scorer reference truth).

Question: does averaging block signs over a target's NEAREST GWPS perturbations (a target-specific, denoised tilt) beat the
marginal per-gene tilt (CLEAN4, block precision 0.5223 on H1) and the delivered signs (0.5083)?
Neighbours: Pearson correlation of the target's own GWPS response row with every other GWPS row over the GWPS gene axis
(the target itself excluded). Votes: the neighbours' GWPS rows (+ optional K562-essential / Jurkat rows of the same
perturbations). Controls: random neighbour sets of the same size (5 draws), the marginal GWPS-only tilt, CLEAN4.
Only CALLED block genes change sign; head keeps delivered signs; magnitudes |E0| (392's construction)."""
import numpy as np, polars as pl, json, sys, time
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/397-knn-tilt/work"; OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
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
M0 = metrics(SG0, np.abs(E0)); rng = np.random.default_rng(20260916)
# ---- CLEAN4 marginal (390's construction) for the known answer
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float); cg = {g: j for j, g in enumerate(cp["genes"].astype(str))}
cmap = np.array([cg.get(g, -1) for g in genes]); cov4 = cmap >= 0
P4 = np.full(nG, np.nan); P4[cov4] = PB[[names.index(n) for n in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]][:, cmap[cov4]].mean(0)
sg4 = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), -1.0)
def block_prec(SG): return float((SG == np.sign(lrm))[cb & adjud].mean())
KA = {"block_precision_E0": block_prec(SG0), "block_precision_CLEAN4": block_prec(np.where(cb, sg4[None, :], SG0))}
print("KNOWN ANSWER (390: E0 0.5083, CLEAN4 0.5223):", json.dumps({k: round(v, 4) for k, v in KA.items()}), flush=True)
# ---- banks
def load_bank(name):
    d = REPO / "data/donors" / name
    D = np.load(d / "delta.npy", mmap_mode="r"); P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    return D, P, G
GW, GWP, GWG = load_bank("k562_gwps"); GW = np.asarray(GW, dtype=np.float32)
gwp = {p: i for i, p in enumerate(GWP)}; gwg = {g: i for i, g in enumerate(GWG)}
gmap = np.array([gwg.get(g, -1) for g in genes]); covg = gmap >= 0
print(f"GWPS {GW.shape}; H1 targets in GWPS {sum(t in gwp for t in targets)}/{nT}; ref genes covered {covg.sum()}/{nG}", flush=True)
extra = {}
for b in ("k562_essential", "nadig_jurkat"):
    D, P, G = load_bank(b); D = np.asarray(D, dtype=np.float32); pi = {p: i for i, p in enumerate(P)}; gi_b = {g: i for i, g in enumerate(G)}
    extra[b] = (D, pi, np.array([gi_b.get(g, -1) for g in genes]))
# neighbour similarity: Pearson over the GWPS gene axis (rows centred)
GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
def votes(rows_idx, banks):
    """P_down over the given GWPS perturbation indices (+ same perturbations' rows in extra banks); NaN where no vote."""
    dn = np.zeros(nG); nn = np.zeros(nG)
    sub = GW[rows_idx][:, np.where(covg, gmap, 0)]; nzm = sub != 0
    dn[covg] += ((sub < 0) & nzm).sum(0)[covg]; nn[covg] += nzm.sum(0)[covg]
    for b in banks:
        D, pi, gm = extra[b]; ridx = [pi[GWP[j]] for j in rows_idx if GWP[j] in pi]
        if not ridx: continue
        cvb = gm >= 0; subb = D[ridx][:, np.where(cvb, gm, 0)]; nzb = subb != 0
        dn[cvb] += ((subb < 0) & nzb).sum(0)[cvb]; nn[cvb] += nzb.sum(0)[cvb]
    P = np.where(nn > 0, dn / np.maximum(nn, 1), np.nan); return P, nn
P_marg_gw, _ = votes(np.arange(len(GWP)), [])
P_marg_3, _ = votes(np.arange(len(GWP)), ["k562_essential", "nadig_jurkat"])
def signs_from(P, fallback):
    return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fallback)
def eval_arm(name, SG, store):
    M = metrics(SG, np.abs(E0)); row = {"block_precision": round(block_prec(SG), 5)}
    for k in M:
        dv = M[k] - M0[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (4000, len(okk)))
        bs = dv[okk][D].mean(1); row[k] = round(float(np.nanmean(dv)), 5); row[f"{k}_ci"] = [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)]
    row["overall_est"] = round(0.561 * row["fid"] - 0.255 * row["nmae"] + 0.188 * row["reach"], 5)
    store[name] = row; print(name, json.dumps(row), flush=True)
out = {"known_answer": KA, "arms": {}}
eval_arm("CLEAN4_marginal", np.where(cb, sg4[None, :], SG0), out["arms"])
eval_arm("GWPS_marginal", np.where(cb, signs_from(P_marg_gw, sg4)[None, :], SG0), out["arms"])
eval_arm("GWPS+ess+jurkat_marginal", np.where(cb, signs_from(P_marg_3, sg4)[None, :], SG0), out["arms"])
# ---- kNN tilts
has = np.array([t in gwp for t in targets])
sims = np.zeros((nT, len(GWP)))
for i, t in enumerate(targets):
    if not has[i]: continue
    v = GWn[gwp[t]]; sims[i] = GWn @ v; sims[i, gwp[t]] = -np.inf   # exclude self
np.save(OUT / "sims.npy", sims.astype(np.float32))
print("top-1 neighbour corr by target (median over targets in GWPS):", round(float(np.median(np.sort(sims[has], 1)[:, -1])), 3), flush=True)
for K in (10, 25, 50, 100, 200, 500):
    for banks_tag, banks in (("gw", []), ("gw+ess+jk", ["k562_essential", "nadig_jurkat"])):
        SG = SG0.copy(); SGs = SG0.copy()
        for i in range(nT):
            if not has[i]: SG[i] = np.where(cb[i], sg4, SG0[i]); SGs[i] = SG[i]; continue
            nb = np.argsort(-sims[i])[:K]; P, nn = votes(nb, banks)
            SG[i] = np.where(cb[i], signs_from(P, sg4), SG0[i])
            # shrinkage toward the marginal tilt with prior weight m = 50 votes
            Pm = P_marg_3 if banks else P_marg_gw; m = 50.0
            Ps = np.where(np.isfinite(P) & np.isfinite(Pm), (nn * np.nan_to_num(P) + m * np.nan_to_num(Pm)) / (nn + m), np.where(np.isfinite(P), P, Pm))
            SGs[i] = np.where(cb[i], signs_from(Ps, sg4), SG0[i])
        eval_arm(f"kNN{K}_{banks_tag}", SG, out["arms"]); eval_arm(f"kNN{K}_{banks_tag}_shrink50", SGs, out["arms"])
    # random-neighbour control at this K (GWPS votes only), 3 draws
    for d in range(3):
        SG = SG0.copy()
        for i in range(nT):
            if not has[i]: SG[i] = np.where(cb[i], sg4, SG0[i]); continue
            nb = rng.choice(len(GWP), K, replace=False); P, nn = votes(nb, [])
            SG[i] = np.where(cb[i], signs_from(P, sg4), SG0[i])
        eval_arm(f"random{K}_gw_draw{d}", SG, out["arms"])
(OUT / "results397.json").write_text(json.dumps(out, indent=1))
print(f"done in {time.time()-t0:.0f}s", flush=True)
