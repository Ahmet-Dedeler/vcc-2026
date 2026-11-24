"""268 stage 1, step 4: leave-one-line-out prediction of the held line's residual from (a) target-only, (b) recipient-only,
(c) target + recipient basal programmes (linear, few parameters), (nn) nearest basal line, with correct / WRONG-recipient /
shuffled-label controls; per-target cosine to the held residual, split-half reproducibility ceiling, paired target bootstrap."""
import json, sys, numpy as np, pandas as pd
from load268 import *
rng = np.random.default_rng(268)
MIN_CPM = 5.0; K = int(sys.argv[1]) if len(sys.argv) > 1 else 2; J = 5; LAM = 1.0; NPERM = 5
MIN_CELLS = 6
H, genes = load_halves(WORKING); lines = WORKING; T = usable_targets(H, lines, MIN_CELLS); nT = len(T); nL = len(lines)
print(f"usable targets {nT}", flush=True)
batch = {l: META.loc[l, "batch"] for l in lines}
def cos(a, b):
    return (a * b).sum(-1) / np.maximum(np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1), 1e-12)
def ridge_fit(Z, Y, lam):  # Z (n,k) standardized, Y (n, ...) -> coefficients (k, ...)
    A = Z.T @ Z + lam * np.eye(Z.shape[1]); return np.linalg.solve(A, Z.T @ Y.reshape(Y.shape[0], -1)).reshape((Z.shape[1],) + Y.shape[1:])
def retrieval(P, truth):
    """PDS-like: for each target, 1 - rank of own truth among all truths by cosine to the prediction /(n-1)."""
    Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12); Tn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    S = Pn @ Tn.T; own = np.diag(S); rank = (S > own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1); return 1 - rank / (len(P) - 1)
rows = []; fold_summary = []
for h in lines:
    train = [l for l in lines if l != h]; hi = lines.index(h); ti = [lines.index(l) for l in train]
    R, ok, avail = build_responses(H, genes, lines, T, keys=("full", "A", "B"), min_cpm=MIN_CPM, filter_lines=train); G = int(ok.sum())
    RF, RA, RB = R["full"], R["A"], R["B"]; assert R["avail"]["full"].all() and R["avail"]["A"].all() and R["avail"]["B"].all()
    m_tr = RF[:, ti].mean(1); r_true = RF[:, hi] - m_tr; truth_full = RF[:, hi]
    rA_h = RA[:, hi] - RA[:, ti].mean(1); rB_h = RB[:, hi] - RB[:, ti].mean(1); ceiling = cos(rA_h, rB_h)
    r_tr = RF[:, ti] - m_tr[:, None]  # (nT, 11, G)
    d_tr = r_tr.mean(0)  # (11, G) line offsets
    e_tr = r_tr - d_tr[None]  # interaction
    d_true = r_true.mean(0)  # ORACLE held line offset (uses held data; reference only)
    # basal programmes fitted on training lines
    C = control_profiles(H, genes, lines, ok); mu = C[ti].mean(0); Ct = C[ti] - mu; ch = C[hi] - mu
    U, s, Vt = np.linalg.svd(Ct, full_matrices=False); comps = Vt[:K]
    Ztr = Ct @ comps.T; sd = Ztr.std(0, ddof=1); Ztr = Ztr / sd; zh = (ch @ comps.T) / sd
    basal_var_share = float((s[:K] ** 2).sum() / (s ** 2).sum())
    # nearest basal line (cosine of centred control profiles)
    simb = cos(Ct, ch[None]); nn = int(np.argmax(simb))
    # programme basis of the interaction (fit on training only)
    E2 = e_tr.transpose(1, 0, 2).reshape(len(ti) * nT, G); _, _, Vte = np.linalg.svd(E2, full_matrices=False); Vp = Vte[:J]  # (J, G)
    S_tr = e_tr @ Vp.T  # (nT, 11, J)
    def predict(Z, z, S=S_tr, e=e_tr, d=d_tr):
        Wd = ridge_fit(Z, d, LAM); dhat = z @ Wd.reshape(K, -1)
        Bc = ridge_fit(Z, e.transpose(1, 0, 2), LAM)  # (K, nT, G)
        c1 = dhat[None] + np.einsum("k,ktg->tg", z, Bc)
        Bs = ridge_fit(Z, S.transpose(1, 0, 2), LAM)  # (K, nT, J)
        c2 = dhat[None] + np.einsum("k,ktj->tj", z, Bs) @ Vp
        return dhat, c1, c2
    dhat, c1, c2 = predict(Ztr, zh)
    preds = {"a_target_only": np.zeros_like(r_true), "b_recipient_offset": np.broadcast_to(dhat, r_true.shape), "c1_target_x_basal_ridge": c1, "c2_target_x_basal_programmes": c2,
             "nn_basal_nearest_line": r_tr[:, nn], "oracle_line_offset": np.broadcast_to(d_true, r_true.shape)}
    # WRONG recipient: every other training line's basal scores / nearest-neighbour replaced by that line's own residual
    wrong = {"b": [], "c1": [], "c2": [], "nn": []}
    for w in range(len(ti)):
        dw, c1w, c2w = predict(Ztr, Ztr[w]); wrong["b"].append(np.broadcast_to(dw, r_true.shape)); wrong["c1"].append(c1w); wrong["c2"].append(c2w); wrong["nn"].append(r_tr[:, w])
    # SHUFFLED line labels: permute which training line carries which basal scores, refit, predict with the correct z_h
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
    for i, t in enumerate(T):
        row = {"held": h, "target": t, "ceiling_cos_AB": float(ceiling[i]), "base_cos_full": float(base_cf[i]), "base_retrieval": float(base_ret[i]), "res_norm_share": float(np.linalg.norm(r_true[i]) ** 2 / max(np.linalg.norm(truth_full[i]) ** 2, 1e-12))}
        for name, (cr, dcf, dret) in rec.items():
            row[f"{name}|cos_res"] = float(cr[i]); row[f"{name}|dcos_full"] = float(dcf[i]); row[f"{name}|dretrieval"] = float(dret[i])
        rows.append(row)
    fs = {"held": h, "batch": batch[h], "genes": G, "basal_var_share_K": basal_var_share, "nn_line": train[nn], "nn_same_batch": batch[train[nn]] == batch[h], "zh": zh.round(2).tolist(),
          "ceiling_median": float(np.median(ceiling)), "base_cos_full_mean": float(base_cf.mean()), "base_retrieval_mean": float(base_ret.mean())}
    for name, (cr, dcf, dret) in rec.items(): fs[name] = (round(float(cr.mean()), 4), round(float(dcf.mean()), 4), round(float(dret.mean()), 4))
    fold_summary.append(fs)
    print(f"{h:9s} G={G} basal share {basal_var_share:.2f} nn={train[nn]}({'same' if fs['nn_same_batch'] else 'diff'} batch) ceiling med {fs['ceiling_median']:.3f} | " + " ".join(f"{k}:{fs[k][0]:+.3f}/{fs[k][2]:+.4f}" for k in ("a_target_only", "b_recipient_offset", "c1_target_x_basal_ridge", "c2_target_x_basal_programmes", "nn_basal_nearest_line", "oracle_line_offset")), flush=True)
    for k, v in rec.items():
        if "WRONG" in k or "SHUFFLED" in k: print(f"          {k}: cos_res {v[0].mean():+.4f} dcos_full {v[1].mean():+.5f} dret {v[2].mean():+.5f}", flush=True)
df = pd.DataFrame(rows); df.to_csv(W / f"lolo268_K{K}.csv", index=False)
# pooled + paired target bootstrap (resample targets; each target carries all 12 folds)
arms = [c[:-8] for c in df.columns if c.endswith("|cos_res")]
tgt_means = {a: df.groupby("target")[[f"{a}|cos_res", f"{a}|dcos_full", f"{a}|dretrieval"]].mean() for a in arms}
def boot(fn, n=2000):
    tg = df.target.unique(); vals = []
    for _ in range(n):
        smp = rng.choice(tg, len(tg), replace=True); vals.append(fn(smp))
    return np.array(vals)
summary = {"K": K, "J": J, "lam": LAM, "n_targets": nT, "n_folds": nL, "ceiling_cos_AB_median": float(df.ceiling_cos_AB.median()), "ceiling_cos_AB_mean": float(df.ceiling_cos_AB.mean()),
           "base_cos_full_mean": float(df.base_cos_full.mean()), "base_retrieval_mean": float(df.base_retrieval.mean()), "res_norm_share_median": float(df.res_norm_share.median()), "arms": {}}
for a in arms:
    m = tgt_means[a]; cr = m[f"{a}|cos_res"]; dcf = m[f"{a}|dcos_full"]; dr = m[f"{a}|dretrieval"]
    b = boot(lambda s: [cr.loc[s].mean(), dcf.loc[s].mean(), dr.loc[s].mean()])
    summary["arms"][a] = {"cos_res": float(cr.mean()), "cos_res_ci": [float(np.quantile(b[:, 0], q)) for q in (0.025, 0.975)], "dcos_full": float(dcf.mean()), "dcos_full_ci": [float(np.quantile(b[:, 1], q)) for q in (0.025, 0.975)],
                          "dretrieval": float(dr.mean()), "dretrieval_ci": [float(np.quantile(b[:, 2], q)) for q in (0.025, 0.975)], "frac_targets_cos_res_pos": float((cr > 0).mean())}
# paired differences: correct vs wrong recipient, correct vs shuffled, c vs a
for a in ("b_recipient_offset", "c1_target_x_basal_ridge", "c2_target_x_basal_programmes", "nn_basal_nearest_line"):
    for ctrl in ("WRONG_recipient", "SHUFFLED_labels"):
        w = f"{a}__{ctrl}"
        if w not in tgt_means: continue
        for met in ("cos_res", "dcos_full", "dretrieval"):
            d = tgt_means[a][f"{a}|{met}"] - tgt_means[w][f"{w}|{met}"]; b = boot(lambda s: d.loc[s].mean())
            summary["arms"][a][f"minus_{ctrl}_{met}"] = [float(d.mean()), float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975)), float((b <= 0).mean())]
summary["folds"] = fold_summary
json.dump(summary, open(W / f"lolo268_K{K}.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in summary.items() if k != "folds"}, indent=1, default=float))
