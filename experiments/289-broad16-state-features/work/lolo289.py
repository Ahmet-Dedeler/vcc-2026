"""289 step 3: leave-one-line-out prediction of the held Broad-16 line's residual from the SOURCE PAPER's state
features, in 268's fold loop (`lolo268.py`), with 268's controls, readouts and paired target bootstrap.

Only the feature matrix Z changes: everything else (gene filter, centring, ridge lam=1, J=5 residual programmes,
split-half ceiling, retrieval, WRONG-recipient and SHUFFLED-label controls) is 268's, refit inside the 11 training
lines of every fold. KA-3: the arm `pc2` (the 2 basal PCs) must reproduce 268's (c1)/(c2) rows.

Also computes the per-target GATING result (PREREGISTRATION.md §7) inside the same fold loop.

Outputs: work/lolo289.csv (fold x target x arm readouts), work/lolo289.json (summary + bootstrap),
work/gating289.csv (target x direction x fold correctness), work/folds289.csv (per-fold means + depth).
"""
import json, sys, numpy as np, pandas as pd
from pathlib import Path

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E268 = ROOT / "experiments/268-broad16-context-residuals/work"
W = ROOT / "experiments/289-broad16-state-features/work"
sys.path.insert(0, str(E268))
from load268 import *  # noqa: WORKING, META, load_halves, build_responses, control_profiles, usable_targets
W = ROOT / "experiments/289-broad16-state-features/work"   # re-assert after the star import (load268 exports its own W)

rng = np.random.default_rng(289)
MIN_CPM = 5.0; KPC = 2; J = 5; LAM = 1.0; NPERM = 5; MIN_CELLS = 6; NBOOT = 4000
GMT = ROOT / "data/external/broad16line/metadata/hallmark_symbols.gmt"

ST = pd.read_csv(W / "state289.csv").set_index("line")
HM = {}
for ln in open(GMT):
    p = ln.rstrip("\n").split("\t"); HM[p[0]] = [g for g in p[2:] if g]
P53SET = set(HM["HALLMARK_P53_PATHWAY"]); E2FSET = set(HM["HALLMARK_E2F_TARGETS"]) | set(HM["HALLMARK_G2M_CHECKPOINT"])

H, genes = load_halves(WORKING); lines = list(WORKING)
T = usable_targets(H, lines, MIN_CELLS); nT = len(T); nL = len(lines)
print(f"usable targets {nT}", flush=True)
batch = {l: META.loc[l, "batch"] for l in lines}
lineage = {l: META.loc[l, "OncotreeLineage"] for l in lines}

def cos(a, b):
    return (a * b).sum(-1) / np.maximum(np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1), 1e-12)

def ridge_fit(Z, Y2, lam):
    """Z (n,k) standardised, Y2 (n, m) already flattened -> (k, m)."""
    A = Z.T @ Z + lam * np.eye(Z.shape[1])
    return np.linalg.solve(A, Z.T @ Y2)

def retrieval(P, truth):
    Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12)
    Tn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    S = Pn @ Tn.T; own = np.diag(S)
    rank = (S > own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)
    return 1 - rank / (len(P) - 1)

# ---- the frozen feature sets (PREREGISTRATION.md §3). `pc` = the fold's basal PCs, injected per fold.
STATE_COLS = ["p53_wt", "hpv", "e2f_g2m", "p53_path", "cn_18q21", "cn_11q23_24"]
ARMS = {
    "pc2": ["PC1", "PC2"],                                   # KA-3: 268's (c1)/(c2)
    "S1": STATE_COLS,
    "S1n": ["p53_wt", "hpv", "e2f_g2m", "p53_path"],
    "S3": STATE_COLS + ["PC1", "PC2"],
    "S4": ["p53_wt_mut", "hpv", "e2f_g2m", "p53_path", "cn_18q21", "cn_11q23_24"],
}
for c in STATE_COLS + ["p53_score", "codel"]:
    ARMS[f"S2_{c}"] = [c]

rows = []; fold_summary = []; gat_rows = []
for h in lines:
    train = [l for l in lines if l != h]; hi = lines.index(h); ti = [lines.index(l) for l in train]
    R, ok, avail = build_responses(H, genes, lines, T, keys=("full", "A", "B"), min_cpm=MIN_CPM, filter_lines=train)
    G = int(ok.sum())
    RF, RA, RB = R["full"], R["A"], R["B"]
    assert R["avail"]["full"].all() and R["avail"]["A"].all() and R["avail"]["B"].all()
    m_tr = RF[:, ti].mean(1); r_true = RF[:, hi] - m_tr; truth_full = RF[:, hi]
    rA_h = RA[:, hi] - RA[:, ti].mean(1); rB_h = RB[:, hi] - RB[:, ti].mean(1); ceiling = cos(rA_h, rB_h)
    r_tr = RF[:, ti] - m_tr[:, None]
    d_tr = r_tr.mean(0); e_tr = r_tr - d_tr[None]; d_true = r_true.mean(0)

    # basal control profiles on this fold's gene axis
    C = control_profiles(H, genes, lines, ok); mu = C[ti].mean(0); Ct = C[ti] - mu; ch = C[hi] - mu
    U, s, Vt = np.linalg.svd(Ct, full_matrices=False); comps = Vt[:KPC]
    basal_var_share = float((s[:KPC] ** 2).sum() / (s ** 2).sum())
    Zpc_tr = Ct @ comps.T; sdpc = Zpc_tr.std(0, ddof=1); Zpc_tr = Zpc_tr / sdpc; zpc_h = (ch @ comps.T) / sdpc
    simb = cos(Ct, ch[None]); nn = int(np.argmax(simb))
    same_lin = [k for k, l in enumerate(train) if lineage[l] == lineage[h]]
    nnL = same_lin[int(np.argmax(simb[same_lin]))] if same_lin else None

    # hallmark scores of the basal state, this fold's gene axis; standardised on the training rows
    gk = genes[ok]; mP = np.isin(gk, list(P53SET)); mE = np.isin(gk, list(E2FSET))
    Craw = C  # (12, G) uncentred log-profiles
    hall = np.stack([Craw[:, mE].mean(1) - Craw.mean(1), Craw[:, mP].mean(1) - Craw.mean(1)], 1)  # (12, 2)

    # full 12-line raw feature table for this fold
    FEAT = {}
    FEAT["e2f_g2m"] = hall[:, 0]; FEAT["p53_path"] = hall[:, 1]
    for c in ("p53_wt", "hpv", "cn_18q21", "cn_11q23_24", "p53_score", "codel", "p53_wt_mut"):
        FEAT[c] = ST.loc[lines, c].to_numpy(float)
    PCcols = {"PC1": (Zpc_tr[:, 0], zpc_h[0]), "PC2": (Zpc_tr[:, 1], zpc_h[1])}

    def build_Z(cols):
        """returns (Ztr (11,k), zh (k,)) standardised on the 11 training rows (PCs already standardised by 268's rule)."""
        tr, hh = [], []
        for c in cols:
            if c in PCcols:
                a, b = PCcols[c]; tr.append(a); hh.append(b)
            else:
                v = FEAT[c]; vt = v[ti]; m = vt.mean(); sd = vt.std(ddof=1)
                sd = sd if sd > 0 else 1.0
                tr.append((vt - m) / sd); hh.append((v[hi] - m) / sd)
        return np.array(tr).T.copy(), np.array(hh)

    # precomputed flattened targets for the ridge
    Yd = d_tr                                        # (11, G)
    Ye = np.ascontiguousarray(e_tr.transpose(1, 0, 2)).reshape(len(ti), -1)   # (11, nT*G)
    _, _, Vte = np.linalg.svd(np.ascontiguousarray(e_tr.transpose(1, 0, 2)).reshape(len(ti) * nT, G), full_matrices=False)
    Vp = Vte[:J]
    S_tr = e_tr @ Vp.T                               # (nT, 11, J)
    Ys = np.ascontiguousarray(S_tr.transpose(1, 0, 2)).reshape(len(ti), -1)   # (11, nT*J)

    def fit_Z(Z, Yd=Yd, Ye=Ye, Ys=Ys):
        """268's three ridge fits; independent of the recipient's feature vector, so fitted once per Z."""
        k = Z.shape[1]
        return (ridge_fit(Z, Yd, LAM), ridge_fit(Z, Ye, LAM).reshape(k, nT, G), ridge_fit(Z, Ys, LAM).reshape(k, nT, J))

    def apply_Z(F, z):
        Wd, Bc, Bs = F
        dhat = z @ Wd
        c1 = dhat[None] + np.einsum("k,ktg->tg", z, Bc)
        c2 = dhat[None] + np.einsum("k,ktj->tj", z, Bs) @ Vp
        return dhat, c1, c2

    def predict(Z, z):
        return apply_Z(fit_Z(Z), z)

    base_cf, base_ret = cos(m_tr, truth_full), retrieval(m_tr, truth_full)
    def metrics(P):
        full = m_tr + P
        return cos(P, r_true), cos(full, truth_full) - base_cf, retrieval(full, truth_full) - base_ret

    rec = {}
    rec["a_target_only"] = metrics(np.zeros_like(r_true))
    rec["oracle_line_offset"] = metrics(np.broadcast_to(d_true, r_true.shape))
    rec["nn_basal_nearest_line"] = metrics(r_tr[:, nn])
    ms = [metrics(r_tr[:, w]) for w in range(len(ti))]
    rec["nn_basal_nearest_line__WRONG_recipient"] = tuple(np.mean([m[i] for m in ms], 0) for i in range(3))
    if nnL is not None:
        rec["nnL_same_lineage_line"] = metrics(r_tr[:, nnL])

    for arm, cols in ARMS.items():
        Ztr, zh = build_Z(cols)
        F = fit_Z(Ztr)
        dh, c1, c2 = apply_Z(F, zh)
        rec[f"{arm}.b"] = metrics(np.broadcast_to(dh, r_true.shape))
        rec[f"{arm}.g"] = metrics(c1); rec[f"{arm}.p"] = metrics(c2)
        wr = {"b": [], "g": [], "p": []}
        for w in range(len(ti)):   # WRONG recipient: same fit, another training line's features
            dw, c1w, c2w = apply_Z(F, Ztr[w])
            wr["b"].append(np.broadcast_to(dw, r_true.shape)); wr["g"].append(c1w); wr["p"].append(c2w)
        sh = {"b": [], "g": [], "p": []}
        for _ in range(NPERM):     # SHUFFLED labels: refit on permuted rows, correct held features
            Zp = Ztr[rng.permutation(len(ti))]
            dp, c1p, c2p = predict(Zp, zh)
            sh["b"].append(np.broadcast_to(dp, r_true.shape)); sh["g"].append(c1p); sh["p"].append(c2p)
        for form in ("b", "g", "p"):
            for cname, lst in (("WRONG_recipient", wr[form]), ("SHUFFLED_labels", sh[form])):
                ms = [metrics(P) for P in lst]
                rec[f"{arm}.{form}__{cname}"] = tuple(np.mean([m[i] for m in ms], 0) for i in range(3))

    # ---- per-target gating (PREREGISTRATION.md §7)
    for dname, mask, cols in (("P53", mP, ["hpv", "p53_wt"]), ("E2F", mE, ["e2f_g2m"])):
        u = mask.astype(np.float64); u = u - u.mean(); u = u / np.linalg.norm(u)
        y_tr = r_tr @ u                      # (nT, 11)
        y_h = r_true @ u                     # (nT,)
        Ztr, zh = build_Z(cols)
        B = ridge_fit(Ztr, y_tr.T, LAM)      # (k, nT)
        pred = zh @ B                        # (nT,)
        shp = []
        for _ in range(NPERM):
            Zp = Ztr[rng.permutation(len(ti))]
            shp.append(zh @ ridge_fit(Zp, y_tr.T, LAM))
        for i, t in enumerate(T):
            gat_rows.append({"held": h, "target": t, "direction": dname,
                             "correct": int(np.sign(pred[i]) == np.sign(y_h[i])),
                             "correct_shuf": float(np.mean([np.sign(p[i]) == np.sign(y_h[i]) for p in shp])),
                             "y_true": float(y_h[i]), "y_pred": float(pred[i])})

    for i, t in enumerate(T):
        row = {"held": h, "target": t, "ceiling_cos_AB": float(ceiling[i]), "base_cos_full": float(base_cf[i]),
               "base_retrieval": float(base_ret[i]),
               "res_norm_share": float(np.linalg.norm(r_true[i]) ** 2 / max(np.linalg.norm(truth_full[i]) ** 2, 1e-12))}
        for name, (cr, dcf, dret) in rec.items():
            row[f"{name}|cos_res"] = float(cr[i]); row[f"{name}|dcos_full"] = float(dcf[i]); row[f"{name}|dretrieval"] = float(dret[i])
        rows.append(row)
    fs = {"held": h, "batch": batch[h], "lineage": lineage[h], "genes": G, "basal_var_share_K": basal_var_share,
          "nn_line": train[nn], "nnL_line": (train[nnL] if nnL is not None else ""),
          "ceiling_median": float(np.median(ceiling)), "base_cos_full_mean": float(base_cf.mean()),
          "base_retrieval_mean": float(base_ret.mean())}
    for c in ("avg_umis_per_ko", "total_umis_per_ko", "n_cells_per_ko"): fs[c] = float(ST.loc[h, c])
    for name, (cr, dcf, dret) in rec.items():
        fs[name] = (round(float(cr.mean()), 4), round(float(dcf.mean()), 4), round(float(dret.mean()), 4))
    fold_summary.append(fs)
    print(f"{h:9s} G={G} ceil {fs['ceiling_median']:.3f} | " +
          " ".join(f"{k}:{fs[k][0]:+.3f}/{fs[k][1]:+.4f}" for k in ("pc2.g", "pc2.p", "S1.g", "S1.p", "S3.g", "S3.p")), flush=True)

df = pd.DataFrame(rows); df.to_csv(W / "lolo289.csv", index=False)
pd.DataFrame(gat_rows).to_csv(W / "gating289.csv", index=False)
pd.DataFrame(fold_summary).to_csv(W / "folds289.csv", index=False)

arms = [c[:-8] for c in df.columns if c.endswith("|cos_res")]
tgt = {a: df.groupby("target")[[f"{a}|cos_res", f"{a}|dcos_full", f"{a}|dretrieval"]].mean() for a in arms}
tg = df.target.unique()
BS = np.array([rng.choice(tg, len(tg), replace=True) for _ in range(NBOOT)])
def boot_mean(series):
    v = series.reindex(tg).to_numpy()
    pos = {t: i for i, t in enumerate(tg)}
    idx = np.vectorize(pos.get)(BS)
    return v[idx].mean(1)
summary = {"K_pc": KPC, "J": J, "lam": LAM, "n_targets": nT, "n_folds": nL, "nboot": NBOOT,
           "ceiling_cos_AB_median": float(df.ceiling_cos_AB.median()), "ceiling_cos_AB_mean": float(df.ceiling_cos_AB.mean()),
           "base_cos_full_mean": float(df.base_cos_full.mean()), "base_retrieval_mean": float(df.base_retrieval.mean()), "arms": {}}
for a in arms:
    m = tgt[a]; e = {}
    for met in ("cos_res", "dcos_full", "dretrieval"):
        s = m[f"{a}|{met}"]; b = boot_mean(s)
        e[met] = float(s.mean()); e[met + "_ci"] = [float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))]
        e[met + "_P_le0"] = float((b <= 0).mean())
    e["folds_pos_cos_res"] = int((df.groupby("held")[f"{a}|cos_res"].mean() > 0).sum())
    e["folds_pos_dcos_full"] = int((df.groupby("held")[f"{a}|dcos_full"].mean() > 0).sum())
    e["frac_targets_cos_res_pos"] = float((m[f"{a}|cos_res"] > 0).mean())
    summary["arms"][a] = e
for a in arms:
    for ctrl in ("WRONG_recipient", "SHUFFLED_labels"):
        w = f"{a}__{ctrl}"
        if w not in tgt: continue
        for met in ("cos_res", "dcos_full", "dretrieval"):
            d = tgt[a][f"{a}|{met}"] - tgt[w][f"{w}|{met}"]; b = boot_mean(d)
            summary["arms"][a][f"minus_{ctrl}_{met}"] = [float(d.mean()), float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975)), float((b <= 0).mean())]
    d = tgt[a][f"{a}|dcos_full"] - tgt["a_target_only"]["a_target_only|dcos_full"]
    summary["arms"][a]["vs_target_only_dcos_full"] = float(d.mean())
summary["folds"] = fold_summary
json.dump(summary, open(W / "lolo289.json", "w"), indent=1, default=float)
print("\nARMS")
for a in sorted(arms):
    e = summary["arms"][a]
    print(f"{a:42s} cos_res {e['cos_res']:+.4f} [{e['cos_res_ci'][0]:+.4f},{e['cos_res_ci'][1]:+.4f}]  "
          f"dcos {e['dcos_full']:+.4f} [{e['dcos_full_ci'][0]:+.4f},{e['dcos_full_ci'][1]:+.4f}]  "
          f"dret {e['dretrieval']:+.4f} [{e['dretrieval_ci'][0]:+.4f},{e['dretrieval_ci'][1]:+.4f}]  folds+ {e['folds_pos_cos_res']}/12")
