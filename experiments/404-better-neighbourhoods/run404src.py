#!/usr/bin/env python3
"""404 -- neighbour SOURCE variants: neighbours by GWPS rows vs by a recipient-like
bank's rows (shared genes, GWPS-overlap perts only) vs union of both top-25 lists.
Votes stay the six-bank family minus recipient; shrink 400. All arms scored on the
same target subset (L targets in GWPS∩S); Deltas vs subset baseline (a)K25 with
paired 10k bootstrap CI. Surrogates: Jurkat->nadig_hepg2, HepG2->orion_hct116,
RPE1->orion_hek293t (+ swaps Jurkat->k562_essential, HepG2->hek293t, RPE1->hct116).
"""
import numpy as np, json, time
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026")
OUT = REPO / "experiments/404-better-neighbourhoods/work/src"
OUT.mkdir(parents=True, exist_ok=True)

def load_bank(name):
    d = (REPO / "data/donors" / name) if name != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
    P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    B = np.load(d / "b_ctrl.npy") if (d / "b_ctrl.npy").exists() else None
    return D, P, G, B

t0 = time.time()
GWm, GWP, GWG, _ = load_bank("k562_gwps")
GW = np.asarray(GWm, dtype=np.float32)
gwp = {p: i for i, p in enumerate(GWP)}; gwg = {g: i for i, g in enumerate(GWG)}
GWc = GW - GW.mean(1, keepdims=True)
GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
banks = {}
for b in ("k562_essential", "nadig_jurkat", "rpe1", "nadig_hepg2", "orion_hct116", "orion_hek293t"):
    banks[b] = load_bank(b); print(b, banks[b][0].shape, flush=True)

rng_neigh = np.random.default_rng(20260918)
rng_boot = np.random.default_rng(20260916)
import warnings; warnings.filterwarnings("ignore", category=RuntimeWarning)

def rownorm(X):
    Xc = X - X.mean(1, keepdims=True)
    return Xc / np.maximum(np.linalg.norm(Xc, axis=1, keepdims=True), 1e-12)

PAIRS = [("nadig_jurkat", "nadig_hepg2"), ("nadig_hepg2", "orion_hct116"), ("rpe1", "orion_hek293t"),
         ("nadig_jurkat", "k562_essential"), ("nadig_hepg2", "orion_hek293t"), ("rpe1", "orion_hct116")]
results = {}
for L, S in PAIRS:
    DL, PL, GL, BL = banks[L]
    DS, PS, GS, _ = banks[S]
    order = np.argsort(-BL, kind="stable"); gl = {g: i for i, g in enumerate(GL)}
    tset = set(PL)
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
    # search space in S: shared genes, overlap perts
    sg_genes = [g for g in GS if g in gwg]
    s_gi = {g: i for i, g in enumerate(GS)}
    s_cols = np.array([s_gi[g] for g in sg_genes]); gw_cols = np.array([gwg[g] for g in sg_genes])
    ps_set = set(PS.tolist())
    ov = [p for p in GWP if p in ps_set]  # GWPS-ordered overlap perts
    s_pi = {p: i for i, p in enumerate(PS)}
    s_rows = np.array([s_pi[p] for p in ov]); gw_rows_ov = np.array([gwp[p] for p in ov])
    pos_in_ov = {p: k for k, p in enumerate(ov)}
    # load S submatrix dense (orion: ~9125 x 7976 = 291MB float32 -- ok, freed per pair)
    DSm = np.asarray(DS[s_rows][:, s_cols], dtype=np.float32)
    print(f"{L}<-{S}: block {len(block)}, shared genes {len(sg_genes)}, overlap perts {len(ov)}", flush=True)
    # targets: L perts in GWPS and in S
    targets = [i for i, t in enumerate(PL) if t in gwp and t in pos_in_ov]
    print(f"  {len(targets)}/{len(PL)} targets usable", flush=True)
    truth = DL[targets][:, bl_idx]; T_nz = truth != 0; T_sign = np.sign(truth); T_w = np.abs(truth)

    def acc(pred_sign):
        ok = T_nz; m = (pred_sign == T_sign) & ok
        per_t = m.sum(1) / np.maximum(ok.sum(1), 1); per_tw = (m * T_w).sum(1) / np.maximum((ok * T_w).sum(1), 1e-12)
        return per_t, per_tw

    def sg1(P, fb): return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fb)
    fb_marg = sg1(P_marg, -1.0)

    # GWPS-side Pearson neighbours (recomputed for subset; same metric as (a))
    tg = np.array([gwp[PL[i]] for i in targets])
    sims_g = GWn @ GWn[tg].T
    G25 = np.zeros((len(targets), 25), dtype=np.int32); G50 = np.zeros((len(targets), 50), dtype=np.int32)
    for j, i in enumerate(targets):
        s_ = sims_g[:, j].copy(); s_[tg[j]] = -np.inf
        nb = np.argsort(-s_)[:50]; G25[j] = nb[:25]; G50[j] = nb
    del sims_g

    # S-side neighbours: Pearson (p), z-scored cosine (z), query-top500 Pearson (t)
    Sn_p = rownorm(DSm)
    mu = DSm.mean(0); sd = DSm.std(0); sd[sd == 0] = 1.0
    Sn_z = rownorm(((DSm - mu) / sd).astype(np.float32)); del mu, sd
    q_all = np.array([pos_in_ov[PL[i]] for i in targets])  # query rows in DSm
    sims_sp = Sn_p @ Sn_p[q_all].T; sims_sz = Sn_z @ Sn_z[q_all].T
    del Sn_p, Sn_z
    S25p = np.zeros((len(targets), 25), dtype=np.int32); S25z = np.zeros_like(S25p); S25t = np.zeros_like(S25p)
    GS25p = np.zeros_like(S25p)  # mapped to GWPS row indices
    GS25z = np.zeros_like(S25p); GS25t = np.zeros_like(S25p)
    for j, i in enumerate(targets):
        q = q_all[j]
        for sims, M, GM in ((sims_sp, S25p, GS25p), (sims_sz, S25z, GS25z)):
            s_ = sims[:, j].copy(); s_[q] = -np.inf
            nb = np.argsort(-s_)[:25]; M[j] = nb; GM[j] = gw_rows_ov[nb]
        # top500 by query |delta| in S over shared genes
        qq = DSm[q]; idx = np.argpartition(-np.abs(qq), 500)[:500]
        Q = DSm[:, idx]; Qc = Q - Q.mean(1, keepdims=True)
        Qn = Qc / np.maximum(np.linalg.norm(Qc, axis=1, keepdims=True), 1e-12)
        qc = qq[idx] - qq[idx].mean(); qn = qc / max(np.linalg.norm(qc), 1e-12)
        s_ = Qn @ qn; s_[q] = -np.inf
        nb = np.argsort(-s_)[:25]; S25t[j] = nb; GS25t[j] = gw_rows_ov[nb]
    del sims_sp, sims_sz, DSm

    def knn_arm(nb_gwrows, mm=400.0):
        Smat = np.zeros((len(targets), len(block)))
        for j in range(len(targets)):
            P, nn = votes(nb_gwrows[j])
            Psh = np.where(np.isfinite(P) & np.isfinite(P_marg),
                           (nn * np.nan_to_num(P) + mm * np.nan_to_num(P_marg)) / (nn + mm),
                           np.where(np.isfinite(P), P, P_marg))
            Smat[j] = sg1(Psh, -1.0)
        return Smat

    U_p = [np.union1d(G25[j], GS25p[j]) for j in range(len(targets))]
    U_z = [np.union1d(G25[j], GS25z[j]) for j in range(len(targets))]
    print(f"  union sizes: p {np.mean([len(u) for u in U_p]):.1f}, z {np.mean([len(u) for u in U_z]):.1f}", flush=True)
    arms = {"marginal": np.tile(fb_marg, (len(targets), 1)),
            "G_a25": knn_arm(G25), "G_a50": knn_arm(G50),
            "S_p25": knn_arm(GS25p), "S_z25": knn_arm(GS25z), "S_t25": knn_arm(GS25t),
            "U_p": knn_arm(U_p), "U_z": knn_arm(U_z)}
    for K, nm in ((25, "random25"), (50, "random50")):
        Rm = np.zeros((len(targets), len(block)))
        for j in range(len(targets)):
            rb = rng_neigh.choice(len(GWP), K, replace=False); P2, _ = votes(rb)
            Rm[j] = sg1(P2, fb_marg)
        arms[nm] = Rm

    base_u, base_w = acc(arms["G_a25"])
    marg_u, _ = acc(arms["marginal"])
    nT = len(targets)
    D = rng_boot.integers(0, nT, (10000, nT)).astype(np.int32)
    res = {"n_targets": nT, "n_overlap_perts": len(ov), "n_shared_genes": len(sg_genes)}
    for name, Smat in arms.items():
        u, w = acc(Smat)
        du = u - base_u; dw = w - base_w; dm = u - marg_u
        bu = du[D].mean(1); bw = dw[D].mean(1); bm = dm[D].mean(1)
        res[name] = {"acc_unw": round(float(u.mean()), 5), "acc_w": round(float(w.mean()), 5),
                     "d_vs_Ga25_unw": round(float(du.mean()), 5),
                     "ci_vs_Ga25": [round(float(np.quantile(bu, .025)), 5), round(float(np.quantile(bu, .975)), 5)],
                     "d_vs_Ga25_w": round(float(dw.mean()), 5),
                     "d_vs_marg_unw": round(float(dm.mean()), 5)}
        print(f"  {L}<-{S} {name:10s} unw {u.mean():.4f} (vsGa25 {du.mean():+.4f} {res[name]['ci_vs_Ga25']}) w {w.mean():.4f} ({dw.mean():+.4f})", flush=True)
    results[f"{L}<-{S}"] = res
(OUT / "results404src.json").write_text(json.dumps(results, indent=1))
print(f"done {time.time()-t0:.0f}s", flush=True)
