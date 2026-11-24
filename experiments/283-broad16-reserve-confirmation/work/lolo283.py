"""283: 268's leave-one-line-out protocol (`lolo268.py`), run in two modes.

  KA       — folds = the 12 WORKING lines, train = the other 11.  Must reproduce 268's lolo268_K2.csv at 0.0e+00.
  RESERVE  — folds = the 4 RESERVED lines, train = all 12 working lines.  The preregistered single access.

The main-loop body is character-identical to lolo268.py; everything new (the permutation null, the OR-control null,
the guide-half residual) draws from a SEPARATE rng so the `rng` stream that feeds the SHUFFLED_labels arms is
untouched and the KA reproduces exactly.   usage: lolo283.py KA|RESERVE
"""
import json, sys, numpy as np, pandas as pd
from lib283 import *
rng = np.random.default_rng(268)          # 268's stream: SHUFFLED_labels + bootstrap.  Do not consume elsewhere.
rng2 = np.random.default_rng(283)         # 283's stream: the new secondaries only.
MODE = sys.argv[1] if len(sys.argv) > 1 else "KA"
MIN_CPM = 5.0; K = 2; J = 5; LAM = 1.0; NPERM = 5
MIN_CELLS = 6
Hw, genes = load_halves(WORKING)
T82 = usable_targets(Hw, WORKING, MIN_CELLS)          # 268's frozen 82 targets
print(f"MODE={MODE}  frozen usable targets (working lines) {len(T82)}", flush=True)
FOLDS = WORKING if MODE == "KA" else RESERVED
batchmap = {l: META.loc[l, "batch"] for l in WORKING + RESERVED}
def cos(a, b):
    return (a * b).sum(-1) / np.maximum(np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1), 1e-12)
def ridge_fit(Z, Y, lam):  # Z (n,k) standardized, Y (n, ...) -> coefficients (k, ...)
    A = Z.T @ Z + lam * np.eye(Z.shape[1]); return np.linalg.solve(A, Z.T @ Y.reshape(Y.shape[0], -1)).reshape((Z.shape[1],) + Y.shape[1:])
def retrieval(P, truth):
    """PDS-like: for each target, 1 - rank of own truth among all truths by cosine to the prediction /(n-1)."""
    Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12); Tn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    S = Pn @ Tn.T; own = np.diag(S); rank = (S > own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1); return 1 - rank / (len(P) - 1)
rows = []; fold_summary = []
for h in FOLDS:
    if MODE == "KA":
        H = Hw; lines = list(WORKING); train = [l for l in lines if l != h]; T = list(T82); dropped = []
    else:
        H, genes = load_halves(WORKING + [h]); lines = list(WORKING) + [h]; train = list(WORKING)
        nc = {t: int(H[h]["ncell"][H[h]["targets"].index(t), 0]) for t in T82 if t in H[h]["targets"]}
        T = [t for t in T82 if nc.get(t, 0) >= MIN_CELLS]; dropped = [(t, nc.get(t, 0)) for t in T82 if nc.get(t, 0) < MIN_CELLS]
    nT = len(T); hi = lines.index(h); ti = [lines.index(l) for l in train]
    R, ok, avail = build_responses(H, genes, lines, T, keys=("full", "A", "B"), min_cpm=MIN_CPM, filter_lines=train); G = int(ok.sum())
    RF, RA, RB = R["full"], R["A"], R["B"]; assert R["avail"]["full"].all() and R["avail"]["A"].all() and R["avail"]["B"].all()
    m_tr = RF[:, ti].mean(1); r_true = RF[:, hi] - m_tr; truth_full = RF[:, hi]
    rA_h = RA[:, hi] - RA[:, ti].mean(1); rB_h = RB[:, hi] - RB[:, ti].mean(1); ceiling = cos(rA_h, rB_h)
    r_tr = RF[:, ti] - m_tr[:, None]  # (nT, |train|, G)
    d_tr = r_tr.mean(0)
    e_tr = r_tr - d_tr[None]
    d_true = r_true.mean(0)  # ORACLE held line offset (uses held data; reference only)
    C = control_profiles(H, genes, lines, ok); mu = C[ti].mean(0); Ct = C[ti] - mu; ch = C[hi] - mu
    U, s, Vt = np.linalg.svd(Ct, full_matrices=False); comps = Vt[:K]
    Ztr = Ct @ comps.T; sd = Ztr.std(0, ddof=1); Ztr = Ztr / sd; zh = (ch @ comps.T) / sd
    basal_var_share = float((s[:K] ** 2).sum() / (s ** 2).sum())
    simb = cos(Ct, ch[None]); nn = int(np.argmax(simb))
    E2 = e_tr.transpose(1, 0, 2).reshape(len(ti) * nT, G); _, _, Vte = np.linalg.svd(E2, full_matrices=False); Vp = Vte[:J]
    S_tr = e_tr @ Vp.T
    def predict(Z, z, S=S_tr, e=e_tr, d=d_tr):
        Wd = ridge_fit(Z, d, LAM); dhat = z @ Wd.reshape(K, -1)
        Bc = ridge_fit(Z, e.transpose(1, 0, 2), LAM)
        c1 = dhat[None] + np.einsum("k,ktg->tg", z, Bc)
        Bs = ridge_fit(Z, S.transpose(1, 0, 2), LAM)
        c2 = dhat[None] + np.einsum("k,ktj->tj", z, Bs) @ Vp
        return dhat, c1, c2
    dhat, c1, c2 = predict(Ztr, zh)
    preds = {"a_target_only": np.zeros_like(r_true), "b_recipient_offset": np.broadcast_to(dhat, r_true.shape), "c1_target_x_basal_ridge": c1, "c2_target_x_basal_programmes": c2,
             "nn_basal_nearest_line": r_tr[:, nn], "oracle_line_offset": np.broadcast_to(d_true, r_true.shape)}
    wrong = {"b": [], "c1": [], "c2": [], "nn": []}
    for w in range(len(ti)):
        dw, c1w, c2w = predict(Ztr, Ztr[w]); wrong["b"].append(np.broadcast_to(dw, r_true.shape)); wrong["c1"].append(c1w); wrong["c2"].append(c2w); wrong["nn"].append(r_tr[:, w])
    shuf = {"b": [], "c1": [], "c2": []}
    for p in range(NPERM):
        Zp = Ztr[rng.permutation(len(ti))]; dp, c1p, c2p = predict(Zp, zh); shuf["b"].append(np.broadcast_to(dp, r_true.shape)); shuf["c1"].append(c1p); shuf["c2"].append(c2p)
    def metrics(P):
        full = m_tr + P; return cos(P, r_true), cos(full, truth_full), retrieval(full, truth_full)
    base_cf, base_ret = cos(m_tr, truth_full), retrieval(m_tr, truth_full)
    rec = {}
    for name, P in preds.items():
        cr, cf, ret = metrics(P); rec[name] = (cr, cf - base_cf, ret - base_ret)
    for key, name in (("b", "b_recipient_offset"), ("c1", "c1_target_x_basal_ridge"), ("c2", "c2_target_x_basal_programmes"), ("nn", "nn_basal_nearest_line")):
        ms = [metrics(P) for P in wrong[key]]; rec[name + "__WRONG_recipient"] = tuple(np.mean([m[i] for m in ms], 0) - (0, base_cf, base_ret)[i] for i in range(3))
        if key in shuf:
            ms = [metrics(P) for P in shuf[key]]; rec[name + "__SHUFFLED_labels"] = tuple(np.mean([m[i] for m in ms], 0) - (0, base_cf, base_ret)[i] for i in range(3))
    # ---- 283 secondaries (separate rng; no effect on anything above) ----
    pi = rng2.permutation(nT)
    null_cos = cos(rA_h, rB_h[pi])                       # S2: within-line permutation null of the residual cosine
    ORr, orav = or_control_responses(H, genes, lines, ok, keys=("A", "B"))   # S3: OR-control residual, same construction
    om = orav.all(1)
    orA = ORr["A"][om][:, hi] - ORr["A"][om][:, ti].mean(1); orB = ORr["B"][om][:, hi] - ORr["B"][om][:, ti].mean(1)
    or_cos = cos(orA, orB)
    Rg, okg, _ = build_responses(H, genes, lines, T, keys=("v1", "v2"), min_cpm=MIN_CPM, filter_lines=train)   # S4: guide halves
    gm = (Rg["avail"]["v1"] & Rg["avail"]["v2"]).all(1)
    g1 = Rg["v1"][:, hi] - Rg["v1"][:, ti].mean(1); g2 = Rg["v2"][:, hi] - Rg["v2"][:, ti].mean(1)
    guide_cos = np.where(gm, cos(g1, g2), np.nan)
    for i, t in enumerate(T):
        row = {"held": h, "target": t, "ceiling_cos_AB": float(ceiling[i]), "base_cos_full": float(base_cf[i]), "base_retrieval": float(base_ret[i]),
               "res_norm_share": float(np.linalg.norm(r_true[i]) ** 2 / max(np.linalg.norm(truth_full[i]) ** 2, 1e-12)),
               "null_cos_AB": float(null_cos[i]), "guide_cos_v1v2": float(guide_cos[i]), "n_cells": int(H[h]["ncell"][H[h]["targets"].index(t), 0])}
        for name, (cr, dcf, dret) in rec.items():
            row[f"{name}|cos_res"] = float(cr[i]); row[f"{name}|dcos_full"] = float(dcf[i]); row[f"{name}|dretrieval"] = float(dret[i])
        rows.append(row)
    fs = {"held": h, "batch": batchmap[h], "genes": G, "n_targets": nT, "n_train": len(ti), "dropped_targets": dropped, "basal_var_share_K": basal_var_share,
          "nn_line": train[nn], "nn_same_batch": batchmap[train[nn]] == batchmap[h], "zh": zh.round(2).tolist(),
          "ceiling_median": float(np.median(ceiling)), "base_cos_full_mean": float(base_cf.mean()), "base_cos_full_median": float(np.median(base_cf)), "base_retrieval_mean": float(base_ret.mean()),
          "null_cos_median": float(np.median(null_cos)), "null_cos_p95": float(np.quantile(null_cos, 0.95)),
          "frac_targets_above_null_p95": float((ceiling > np.quantile(null_cos, 0.95)).mean()),
          "or_control_cos_median": float(np.median(or_cos)), "n_or_controls": int(om.sum()),
          "guide_cos_median": float(np.nanmedian(guide_cos)), "n_guide_targets": int(gm.sum()),
          "cells_per_ko_median": float(np.median([H[h]["ncell"][H[h]["targets"].index(t), 0] for t in T]))}
    for name, (cr, dcf, dret) in rec.items(): fs[name] = (round(float(cr.mean()), 4), round(float(dcf.mean()), 4), round(float(dret.mean()), 4))
    fold_summary.append(fs)
    print(f"{h:9s} G={G} nT={nT} basal share {basal_var_share:.2f} nn={train[nn]}({'same' if fs['nn_same_batch'] else 'diff'} batch) ceiling med {fs['ceiling_median']:.4f} "
          f"null med {fs['null_cos_median']:+.4f} p95 {fs['null_cos_p95']:.4f} frac>p95 {fs['frac_targets_above_null_p95']:.3f} OR med {fs['or_control_cos_median']:+.4f} "
          f"guide med {fs['guide_cos_median']:.4f} (n={fs['n_guide_targets']}) base med {fs['base_cos_full_median']:.4f} | "
          + " ".join(f"{k}:{fs[k][0]:+.3f}/{fs[k][1]:+.4f}/{fs[k][2]:+.4f}" for k in ("c1_target_x_basal_ridge", "c2_target_x_basal_programmes")), flush=True)
df = pd.DataFrame(rows); df.to_csv(W / f"lolo283_{MODE}.csv", index=False)
json.dump(fold_summary, open(W / f"lolo283_{MODE}_folds.json", "w"), indent=1, default=float)
print("wrote", W / f"lolo283_{MODE}.csv", len(df), "rows", flush=True)
