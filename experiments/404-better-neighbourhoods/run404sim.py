#!/usr/bin/env python3
"""404 -- similarity variants for the kNN block tilt (LOLO, six-bank family, shrink 400).

(a) Pearson all genes [baseline, code path identical to lolo398d];
(b) Pearson on query row's top-500 |delta| genes;
(c) cosine on rows z-scored per gene across GWPS;
(d) Spearman (average-tie ranks per row);
(e) Pearson on rows projected onto GWPS top-50 PCs.
K in {10,25,50,100} for (a); all sims at K=25; winner gets K/shrink sweep.
Controls: marginal, all-down, random-K. Endpoints: unw + |delta|-weighted block
accuracy; paired per-target bootstrap (10k resamples) for d vs (a)K25 and vs marginal.
"""
import numpy as np, json, time
from pathlib import Path
from scipy.stats import rankdata
REPO = Path("~/Code/virtual-cell-challenge-2026")
OUT = REPO / "experiments/404-better-neighbourhoods/work/sim"
OUT.mkdir(parents=True, exist_ok=True)

def load_bank(name):
    d = (REPO / "data/donors" / name) if name != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
    P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    B = np.load(d / "b_ctrl.npy") if (d / "b_ctrl.npy").exists() else None
    return D, P, G, B

t0 = time.time()
GWm, GWP, GWG, _ = load_bank("k562_gwps")
GW = np.asarray(GWm, dtype=np.float32)  # dense 9787x8248 for row slicing
print(f"GW {GW.shape} zeros_frac {float((GW==0).mean()):.4f}", flush=True)
gwp = {p: i for i, p in enumerate(GWP)}; gwg = {g: i for i, g in enumerate(GWG)}
GWc = GW - GW.mean(1, keepdims=True)
GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)  # (a), bit-identical to 398d
banks = {}
for b in ("k562_essential", "nadig_jurkat", "rpe1", "nadig_hepg2", "orion_hct116", "orion_hek293t"):
    banks[b] = load_bank(b); print(b, banks[b][0].shape, flush=True)

rng_neigh = np.random.default_rng(20260917)
rng_boot = np.random.default_rng(20260916)
results = {}

def rownorm(X):
    Xc = X - X.mean(1, keepdims=True)
    return Xc / np.maximum(np.linalg.norm(Xc, axis=1, keepdims=True), 1e-12)

for L in ("nadig_jurkat", "nadig_hepg2", "rpe1"):
    DL, PL, GL, BL = banks[L]
    shared = [g for g in GL if g in gwg]; tset = set(PL)
    order = np.argsort(-BL, kind="stable"); gl = {g: i for i, g in enumerate(GL)}
    block = [GL[j] for j in order if GL[j] in gwg and GL[j] not in tset][:2200]
    bl_idx = np.array([gl[g] for g in block]); gw_idx = np.array([gwg[g] for g in block])
    vote_banks = [b for b in banks if b != L and b != "h1_2025"]
    vb = []
    for b in vote_banks:
        D, P, G, _ = banks[b]; pi = {p: i for i, p in enumerate(P)}; gi = {g: i for i, g in enumerate(G)}
        vb.append((D, pi, np.array([gi.get(g, -1) for g in block])))

    def votes(rows):
        sub = GW[rows][:, gw_idx]; nz = sub != 0; dn = ((sub < 0) & nz).sum(0).astype(float); nn = nz.sum(0).astype(float)
        for D, pi, gm in vb:
            ridx = [pi[GWP[j]] for j in rows if GWP[j] in pi]
            if not ridx: continue
            cv = gm >= 0; s2 = D[ridx][:, np.where(cv, gm, 0)]; nz2 = (s2 != 0) & cv[None, :]
            dn += ((s2 < 0) & nz2).sum(0); nn += nz2.sum(0)
        return np.where(nn > 0, dn / np.maximum(nn, 1), np.nan), nn

    P_marg, _ = votes(np.arange(len(GWP)))
    targets = [i for i, t in enumerate(PL) if t in gwp]
    tg = np.array([gwp[PL[i]] for i in targets])  # GWPS row per target
    print(f"{L}: block {len(block)}, {len(targets)}/{len(PL)} targets", flush=True)
    truth = DL[targets][:, bl_idx]; T_nz = truth != 0; T_sign = np.sign(truth); T_w = np.abs(truth)

    def acc(pred_sign):
        ok = T_nz; m = (pred_sign == T_sign) & ok
        per_t = m.sum(1) / np.maximum(ok.sum(1), 1); per_tw = (m * T_w).sum(1) / np.maximum((ok * T_w).sum(1), 1e-12)
        return per_t, per_tw

    def sg(P, fb): return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fb)
    fb_marg = sg(P_marg, -1.0)
    arms = {"all_down": np.full((len(targets), len(block)), -1.0),
            "marginal": np.tile(fb_marg, (len(targets), 1))}

    # ---- neighbour lists per similarity (nT x Kmax), self excluded ----
    KMAX = 100
    NB, TOP1 = {}, {}
    # (a) baseline Pearson, identical path to 398d (default quicksort argsort, NaN warning suppressed)
    import warnings; warnings.filterwarnings("ignore", category=RuntimeWarning)
    sims_a = GWn @ GWn[tg].T
    NB["a"] = np.zeros((len(targets), KMAX), dtype=np.int32)
    t1 = np.zeros(len(targets))
    for j, i in enumerate(targets):
        s_ = sims_a[:, j].copy(); s_[gwp[PL[i]]] = -np.inf
        NB["a"][j] = np.argsort(-s_)[:KMAX]; t1[j] = s_[NB["a"][j, 0]]
    TOP1["a"] = float(np.median(t1)); del sims_a
    print(f"  {L} sim-a median top1 r = {TOP1['a']:.4f}", flush=True)

    # (b) Pearson on query top-500 |delta|
    NB["b"] = np.zeros((len(targets), KMAX), dtype=np.int32); t1 = np.zeros(len(targets))
    tb = time.time()
    for j, i in enumerate(targets):
        q = GW[tg[j]]
        idx = np.argpartition(-np.abs(q), 500)[:500]
        Q = GW[:, idx]; Qc = Q - Q.mean(1, keepdims=True)
        Qn = Qc / np.maximum(np.linalg.norm(Qc, axis=1, keepdims=True), 1e-12)
        qc = q[idx] - q[idx].mean(); qn = qc / max(np.linalg.norm(qc), 1e-12)
        s_ = Qn @ qn; s_[tg[j]] = -np.inf
        NB["b"][j] = np.argsort(-s_, kind="stable")[:KMAX]; t1[j] = s_[NB["b"][j, 0]]
    TOP1["b"] = float(np.median(t1)); print(f"  {L} sim-b median top1 r = {TOP1['b']:.4f} ({time.time()-tb:.0f}s)", flush=True)

    # (c) cosine on per-gene z-scored rows
    mu = GW.mean(0); sd = GW.std(0); sd[sd == 0] = 1.0
    Z = rownorm(((GW - mu) / sd).astype(np.float32)); del mu, sd
    sims_c = Z @ Z[tg].T; del Z
    NB["c"] = np.zeros((len(targets), KMAX), dtype=np.int32); t1 = np.zeros(len(targets))
    for j, i in enumerate(targets):
        s_ = sims_c[:, j].copy(); s_[tg[j]] = -np.inf
        NB["c"][j] = np.argsort(-s_, kind="stable")[:KMAX]; t1[j] = s_[NB["c"][j, 0]]
    TOP1["c"] = float(np.median(t1)); del sims_c
    print(f"  {L} sim-c median top1 r = {TOP1['c']:.4f}", flush=True)

    # (d) Spearman
    tb = time.time()
    R = np.empty_like(GW)
    for r in range(GW.shape[0]):
        R[r] = rankdata(GW[r], method="average").astype(np.float32)
    Rn = rownorm(R); del R
    sims_d = Rn @ Rn[tg].T; del Rn
    NB["d"] = np.zeros((len(targets), KMAX), dtype=np.int32); t1 = np.zeros(len(targets))
    for j, i in enumerate(targets):
        s_ = sims_d[:, j].copy(); s_[tg[j]] = -np.inf
        NB["d"][j] = np.argsort(-s_, kind="stable")[:KMAX]; t1[j] = s_[NB["d"][j, 0]]
    TOP1["d"] = float(np.median(t1)); del sims_d
    print(f"  {L} sim-d median top1 rho = {TOP1['d']:.4f} ({time.time()-tb:.0f}s)", flush=True)

    # (e) top-50 PCs
    U = np.load(REPO / "data/donors/k562_gwps/svd_U.npy", mmap_mode="r")[:, :50].astype(np.float32)
    S = np.load(REPO / "data/donors/k562_gwps/svd_S.npy")[:50].astype(np.float32)
    Cn = rownorm(U * S[None, :]); del U, S
    sims_e = Cn @ Cn[tg].T; del Cn
    NB["e"] = np.zeros((len(targets), KMAX), dtype=np.int32); t1 = np.zeros(len(targets))
    for j, i in enumerate(targets):
        s_ = sims_e[:, j].copy(); s_[tg[j]] = -np.inf
        NB["e"][j] = np.argsort(-s_, kind="stable")[:KMAX]; t1[j] = s_[NB["e"][j, 0]]
    TOP1["e"] = float(np.median(t1)); del sims_e
    print(f"  {L} sim-e median top1 r = {TOP1['e']:.4f}", flush=True)

    def knn_arm(nb_mat, K, mm=400.0):
        Smat = np.zeros((len(targets), len(block)))
        for j in range(len(targets)):
            P, nn = votes(nb_mat[j, :K])
            Psh = np.where(np.isfinite(P) & np.isfinite(P_marg),
                           (nn * np.nan_to_num(P) + mm * np.nan_to_num(P_marg)) / (nn + mm),
                           np.where(np.isfinite(P), P, P_marg))
            Smat[j] = sg(Psh, -1.0)
        return Smat

    # phase 1: (a) x K sweep + all sims at K25 + random controls
    arms["a_kNN10"] = knn_arm(NB["a"], 10); arms["a_kNN25"] = knn_arm(NB["a"], 25)
    arms["a_kNN50"] = knn_arm(NB["a"], 50); arms["a_kNN100"] = knn_arm(NB["a"], 100)
    for s in "bcde":
        arms[f"{s}_kNN25"] = knn_arm(NB[s], 25)
    for K, nm in ((25, "random25"), (50, "random50")):
        Rm = np.zeros((len(targets), len(block)))
        for j in range(len(targets)):
            rb = rng_neigh.choice(len(GWP), K, replace=False); P2, _ = votes(rb)
            Rm[j] = sg(P2, fb_marg)
        arms[nm] = Rm

    # phase-1 scoring to pick winner for K/shrink sweep (decided in analysis; run all K/shrink for b..e at 25 only here)
    # phase 2: K/shrink sweep for every sim (cheap enough: 4 sims x 3 K + shrinks)
    for s in "bcde":
        for K in (10, 50, 100):
            arms[f"{s}_kNN{K}"] = knn_arm(NB[s], K)
    for s in "abcde":
        for mm, nm in ((200.0, "sh200"), (800.0, "sh800")):
            Smat = np.zeros((len(targets), len(block)))
            for j in range(len(targets)):
                P, nn = votes(NB[s][j, :25])
                Psh = np.where(np.isfinite(P) & np.isfinite(P_marg),
                               (nn * np.nan_to_num(P) + mm * np.nan_to_num(P_marg)) / (nn + mm),
                               np.where(np.isfinite(P), P, P_marg))
                Smat[j] = sg(Psh, -1.0)
            arms[f"{s}_kNN25_{nm}"] = Smat

    base_u, base_w = acc(arms["a_kNN25"])
    marg_u, _ = acc(arms["marginal"])
    nT = len(targets)
    D = rng_boot.integers(0, nT, (10000, nT)).astype(np.int32)
    res = {"n_targets": nT, "median_top1": TOP1}
    for name, Smat in arms.items():
        u, w = acc(Smat)
        du = u - base_u; dw = w - base_w; dm = u - marg_u
        bu = du[D].mean(1); bw = dw[D].mean(1); bm = dm[D].mean(1)
        res[name] = {"acc_unw": round(float(u.mean()), 5), "acc_w": round(float(w.mean()), 5),
                     "d_vs_a_unw": round(float(du.mean()), 5),
                     "ci_vs_a": [round(float(np.quantile(bu, .025)), 5), round(float(np.quantile(bu, .975)), 5)],
                     "d_vs_a_w": round(float(dw.mean()), 5),
                     "d_vs_marg_unw": round(float(dm.mean()), 5),
                     "ci_vs_marg": [round(float(np.quantile(bm, .025)), 5), round(float(np.quantile(bm, .975)), 5)]}
        print(f"  {L:12s} {name:14s} unw {u.mean():.4f} (vsA {du.mean():+.4f} {res[name]['ci_vs_a']}) w {w.mean():.4f} ({dw.mean():+.4f})", flush=True)
    results[L] = res
(OUT / "results404sim.json").write_text(json.dumps(results, indent=1))
print(f"done {time.time()-t0:.0f}s", flush=True)
