#!/usr/bin/env python3
"""398 -- neighbour-conditioned block tilt, leave-one-bank-out on the public banks (371's framework, 370's block rule).

For held-out bank L (nadig_jurkat, k562_essential, h1_2025, rpe1 if present): block_L = top-2200 genes by L's own b_ctrl
among genes shared with GWPS; for each L target t that is also a GWPS perturbation, neighbours = top-K GWPS rows by Pearson
correlation with t's GWPS row (t itself excluded; when L is a K562 bank the GWPS row of t is still excluded);
votes = neighbours' GWPS rows (+ the same perturbations' rows in the other non-held banks); truth = sign of L's delta for t.
Compared: all-down, marginal tilt from the same vote banks, kNN tilt, kNN shrunk toward the marginal (m=50), random-K control.
Endpoints: unweighted block accuracy (fid reads this) and |delta|-weighted (nmae proxy); per-target paired bootstrap CI for
kNN minus marginal."""
import numpy as np, json, time, sys
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); OUT = REPO / "experiments/398-knn-tilt-lolo/work/orion_m"; OUT.mkdir(parents=True, exist_ok=True)
def load_bank(name):
    d = (REPO / "data/donors" / name) if name != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
    P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    B = np.load(d / "b_ctrl.npy") if (d / "b_ctrl.npy").exists() else None
    return D, P, G, B
t0 = time.time()
GW, GWP, GWG, _ = load_bank("k562_gwps"); gwp = {p: i for i, p in enumerate(GWP)}; gwg = {g: i for i, g in enumerate(GWG)}
GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
banks = {}
for b in ("k562_essential", "nadig_jurkat", "rpe1", "nadig_hepg2", "orion_hct116", "orion_hek293t"):
    try: banks[b] = load_bank(b); print(b, banks[b][0].shape, flush=True)
    except Exception as ex: print("skip", b, ex, flush=True)
rng = np.random.default_rng(20260916)
results = {}
for L in ("rpe1", "nadig_hepg2", "nadig_jurkat"):
    DL, PL, GL, BL = banks[L]
    if BL is None: print("no b_ctrl for", L); continue
    # block: top-2200 by L's control among genes shared with GWPS, minus L's own targets (as the emitter excludes panel targets)
    shared = [g for g in GL if g in gwg]; tset = set(PL)
    order = np.argsort(-BL, kind="stable"); gl = {g: i for i, g in enumerate(GL)}
    block = [GL[j] for j in order if GL[j] in gwg and GL[j] not in tset][:2200]
    bl_idx = np.array([gl[g] for g in block]); gw_idx = np.array([gwg[g] for g in block])
    vote_banks = [b for b in banks if b != L and b != "h1_2025"]   # H1 never votes (it is the H1 harness truth elsewhere)
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
    print(f"{L}: block {len(block)} genes, {len(targets)}/{len(PL)} targets in GWPS, vote banks {vote_banks}", flush=True)
    truth = DL[targets][:, bl_idx]; T_nz = truth != 0; T_sign = np.sign(truth); T_w = np.abs(truth)
    def acc(pred_sign):   # pred_sign: (nT, nG) or (nG,)
        ok = T_nz; m = (pred_sign == T_sign) & ok
        per_t = m.sum(1) / np.maximum(ok.sum(1), 1); per_tw = (m * T_w).sum(1) / np.maximum((ok * T_w).sum(1), 1e-12)
        return per_t, per_tw
    def sg(P, fb): return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fb)
    arms = {"all_down": np.full((len(targets), len(block)), -1.0), "marginal": np.tile(sg(P_marg, -1.0), (len(targets), 1))}
    sims = GWn @ GWn[[gwp[PL[i]] for i in targets]].T   # (9787, nT)
    for K in (25,):
        Pk = {mm: np.zeros((len(targets), len(block))) for mm in (50, 100, 200, 400, 800, 1600)}
        for j, i in enumerate(targets):
            s_ = sims[:, j].copy(); s_[gwp[PL[i]]] = -np.inf
            nb = np.argsort(-s_)[:K]; P, nn = votes(nb)
            for mm in Pk:
                Psh = np.where(np.isfinite(P) & np.isfinite(P_marg), (nn * np.nan_to_num(P) + float(mm) * np.nan_to_num(P_marg)) / (nn + float(mm)), np.where(np.isfinite(P), P, P_marg))
                Pk[mm][j] = sg(Psh, -1.0)
        for mm in Pk: arms[f"kNN{K}_shrink{mm}"] = Pk[mm]
    res = {}
    base_u, base_w = acc(arms["marginal"])
    for name, S in arms.items():
        u, w = acc(S); du = u - base_u; dw = w - base_w
        D = rng.integers(0, len(u), (4000, len(u))); bu = du[D].mean(1); bw = dw[D].mean(1)
        res[name] = {"acc_unw": round(float(u.mean()), 5), "acc_w": round(float(w.mean()), 5),
                     "d_vs_marginal_unw": round(float(du.mean()), 5), "ci_unw": [round(float(np.quantile(bu, .025)), 5), round(float(np.quantile(bu, .975)), 5)],
                     "d_vs_marginal_w": round(float(dw.mean()), 5), "ci_w": [round(float(np.quantile(bw, .025)), 5), round(float(np.quantile(bw, .975)), 5)]}
        print(f"  {L:16s} {name:18s} unw {res[name]['acc_unw']:.4f} ({res[name]['d_vs_marginal_unw']:+.4f} {res[name]['ci_unw']})  w {res[name]['acc_w']:.4f} ({res[name]['d_vs_marginal_w']:+.4f})", flush=True)
    results[L] = {"n_targets": len(targets), "n_block": len(block), "vote_banks": vote_banks, "arms": res}
(OUT / "results398.json").write_text(json.dumps(results, indent=1)); print(f"done {time.time()-t0:.0f}s", flush=True)
