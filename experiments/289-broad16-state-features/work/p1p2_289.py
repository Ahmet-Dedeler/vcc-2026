"""289 step 4: the source paper's two mechanistic predictions, tested on OUR residual object
(PREREGISTRATION.md §5 P1, §6 P2), plus the per-target gating with per-permutation shuffled controls (§7).

Outputs work/p1p2_289.json, work/gating289b.csv, work/p1_scores289.csv.
"""
import json, sys, itertools, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E268 = ROOT / "experiments/268-broad16-context-residuals/work"
sys.path.insert(0, str(E268))
from load268 import *  # noqa
W = ROOT / "experiments/289-broad16-state-features/work"

rng = np.random.default_rng(289)
MIN_CPM = 5.0; LAM = 1.0; NPERM = 5; MIN_CELLS = 6
ST = pd.read_csv(W / "state289.csv").set_index("line")
HM = {}
for ln in open(ROOT / "data/external/broad16line/metadata/hallmark_symbols.gmt"):
    q = ln.rstrip("\n").split("\t"); HM[q[0]] = [g for g in q[2:] if g]
P53SET = set(HM["HALLMARK_P53_PATHWAY"]); E2FSET = set(HM["HALLMARK_E2F_TARGETS"]) | set(HM["HALLMARK_G2M_CHECKPOINT"])

H, genes = load_halves(WORKING); lines = list(WORKING); nL = len(lines)
T82 = usable_targets(H, lines, MIN_CELLS)
out = {"lines": lines, "n_targets_82": len(T82)}

# =========================================================== P1
P1T = ["UBE3A", "MDM2"]
ncell = {t: {l: int(H[l]["ncell"][H[l]["targets"].index(t), 0]) if t in H[l]["targets"] else 0 for l in lines} for t in P1T}
out["P1_cells_per_line"] = ncell
R, ok, avail = build_responses(H, genes, lines, P1T, keys=("full",), min_cpm=MIN_CPM)
gk = genes[ok]; mP = np.isin(gk, list(P53SET)); mE = np.isin(gk, list(E2FSET))
RF = R["full"]; av = R["avail"]["full"]
out["P1_gene_axis"] = {"G": int(ok.sum()), "n_p53": int(mP.sum()), "n_e2f": int(mE.sum())}

scores = {}
for i, t in enumerate(P1T):
    incl = [j for j, l in enumerate(lines) if av[i, j] and ncell[t][l] >= MIN_CELLS]
    r = RF[i, incl] - RF[i, incl].mean(0)                    # residual against the included-line mean
    scores[t] = {"lines": [lines[j] for j in incl],
                 "p53": (r[:, mP].mean(1) - r.mean(1)).astype(float),
                 "e2f": (r[:, mE].mean(1) - r.mean(1)).astype(float)}
    out.setdefault("P1_included", {})[t] = [lines[j] for j in incl]

common = [l for l in scores[P1T[0]]["lines"] if l in scores[P1T[1]]["lines"]]
def stack(key):
    M = []
    for t in P1T:
        idx = [scores[t]["lines"].index(l) for l in common]
        M.append(scores[t][key][idx])
    return np.array(M)                                        # (2 targets, n_common)
P53M, E2FM = stack("p53"), stack("e2f")
prim = P53M.mean(0)                                           # mean of the two targets, per line
out["P1_n_common_lines"] = len(common); out["P1_common_lines"] = common

hpv = [l for l in common if ST.loc[l, "hpv"] == 1]
wt = [l for l in common if ST.loc[l, "hpv"] == 0 and ST.loc[l, "p53_wt"] == 1]
dfc = [l for l in common if ST.loc[l, "hpv"] == 0 and ST.loc[l, "p53_wt"] == 0]
gi = lambda names: [common.index(l) for l in names]
obs = float(prim[gi(hpv)].mean())
pairs = list(itertools.combinations(range(len(common)), 2))
null = np.array([prim[list(p)].mean() for p in pairs])
p_perm = float((null >= obs).mean())
out["P1"] = {
    "targets": P1T, "statistic": "mean over UBE3A,MDM2 of the HALLMARK_P53_PATHWAY score of the line-specific residual",
    "per_line": {l: float(prim[k]) for k, l in enumerate(common)},
    "per_line_per_target": {t: {l: float(P53M[i, k]) for k, l in enumerate(common)} for i, t in enumerate(P1T)},
    "hpv_lines": hpv, "p53wt_nonhpv_lines": wt, "p53def_lines": dfc,
    "mean_hpv": obs, "mean_p53wt_nonhpv": float(prim[gi(wt)].mean()) if wt else None,
    "mean_p53def": float(prim[gi(dfc)].mean()) if dfc else None,
    "n_pairs": len(pairs), "p_perm_one_sided": p_perm,
    "rank_of_hpv_pair": int((null >= obs).sum()),
    "order_hpv_gt_wt": bool(wt and obs > prim[gi(wt)].mean()),
    "order_hpv_gt_def": bool(dfc and obs > prim[gi(dfc)].mean()),
}
out["P1"]["PASS"] = bool(p_perm <= 0.05 and out["P1"]["order_hpv_gt_wt"] and out["P1"]["order_hpv_gt_def"])
# per-target secondary
for i, t in enumerate(P1T):
    o = float(P53M[i][gi(hpv)].mean()); n = np.array([P53M[i][list(p)].mean() for p in pairs])
    out["P1"][f"{t}_only"] = {"mean_hpv": o, "p_perm": float((n >= o).mean()),
                              "mean_p53wt_nonhpv": float(P53M[i][gi(wt)].mean()) if wt else None,
                              "mean_p53def": float(P53M[i][gi(dfc)].mean()) if dfc else None}
# specificity: the same statistic on the E2F/G2M direction
eprim = E2FM.mean(0); eo = float(eprim[gi(hpv)].mean()); en = np.array([eprim[list(p)].mean() for p in pairs])
out["P1_specificity_E2F"] = {"per_line": {l: float(eprim[k]) for k, l in enumerate(common)},
                             "mean_hpv": eo, "p_perm": float((en >= eo).mean())}
pd.DataFrame({"line": common, "p53_score_mean": prim, "e2f_score_mean": eprim,
              **{f"p53_{t}": P53M[i] for i, t in enumerate(P1T)},
              **{f"e2f_{t}": E2FM[i] for i, t in enumerate(P1T)}}).to_csv(W / "p1_scores289.csv", index=False)
print("P1", json.dumps({k: v for k, v in out["P1"].items() if k not in ("per_line", "per_line_per_target")}, indent=1, default=float), flush=True)

# =========================================================== P2
R2, ok2, av2 = build_responses(H, genes, lines, T82, keys=("full",), min_cpm=MIN_CPM)
RF2 = R2["full"]                                    # (82, 12, G)
res2 = RF2 - RF2.mean(1, keepdims=True)
i_ier = T82.index("IER3IP1")
nrm = np.linalg.norm(RF2, axis=2)                   # (82, 12)
rnrm = np.linalg.norm(res2, axis=2)
rel_resp = nrm[i_ier] / np.median(nrm, axis=0)
rel_res = rnrm[i_ier] / np.median(rnrm, axis=0)
codel = ST.loc[lines, "codel"].to_numpy(int)
cscore = (ST.loc[lines, "cn_18q21"].to_numpy(float) + ST.loc[lines, "cn_11q23_24"].to_numpy(float)) / 2
fitness = ST.loc[lines, "IER3IP1_effect"].to_numpy(float)
def mw(x, alt):  # one-sided exact Mann-Whitney, codel group vs the rest
    a, b = x[codel == 1], x[codel == 0]
    if min(len(a), len(b)) < 2: return None
    return float(stats.mannwhitneyu(a, b, alternative=alt).pvalue)
p2 = {"codel_lines": [l for l, c in zip(lines, codel) if c], "n_codel": int(codel.sum()),
      "per_line": {l: {"fitness_IER3IP1": float(fitness[j]), "rel_response_norm": float(rel_resp[j]),
                       "rel_residual_norm": float(rel_res[j]), "codel": int(codel[j]), "codel_score": float(cscore[j])}
                   for j, l in enumerate(lines)},
      "mean_fitness_codel": float(fitness[codel == 1].mean()), "mean_fitness_other": float(fitness[codel == 0].mean()),
      "mean_relresp_codel": float(rel_resp[codel == 1].mean()), "mean_relresp_other": float(rel_resp[codel == 0].mean()),
      "mean_relres_codel": float(rel_res[codel == 1].mean()), "mean_relres_other": float(rel_res[codel == 0].mean()),
      "mw_fitness_p_less": mw(fitness, "less"), "mw_relresp_p_greater": mw(rel_resp, "greater"),
      "mw_relres_p_greater": mw(rel_res, "greater"),
      "spearman_cscore_fitness": [float(v) for v in stats.spearmanr(cscore, fitness)],
      "spearman_cscore_relresp": [float(v) for v in stats.spearmanr(cscore, rel_resp)],
      "spearman_cscore_relres": [float(v) for v in stats.spearmanr(cscore, rel_res)]}
p2["PASS"] = bool(p2["mw_fitness_p_less"] is not None and p2["mw_fitness_p_less"] <= 0.05
                  and p2["mean_relresp_codel"] > p2["mean_relresp_other"])
out["P2"] = p2
print("P2", json.dumps({k: v for k, v in p2.items() if k != "per_line"}, indent=1, default=float), flush=True)

# =========================================================== gating with per-permutation shuffles
def ridge(Z, Y2, lam=LAM):
    return np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ Y2)
gat = []
for h in lines:
    train = [l for l in lines if l != h]; hi = lines.index(h); ti = [lines.index(l) for l in train]
    R3, ok3, _ = build_responses(H, genes, lines, T82, keys=("full",), min_cpm=MIN_CPM, filter_lines=train)
    RF3 = R3["full"]; m_tr = RF3[:, ti].mean(1)
    r_true = RF3[:, hi] - m_tr; r_tr = RF3[:, ti] - m_tr[:, None]
    C = control_profiles(H, genes, lines, ok3)
    gk3 = genes[ok3]; mP3 = np.isin(gk3, list(P53SET)); mE3 = np.isin(gk3, list(E2FSET))
    FEAT = {"e2f_g2m": C[:, mE3].mean(1) - C.mean(1),
            "hpv": ST.loc[lines, "hpv"].to_numpy(float), "p53_wt": ST.loc[lines, "p53_wt"].to_numpy(float)}
    def Z_of(cols):
        tr, hh = [], []
        for c in cols:
            v = FEAT[c]; vt = v[ti]; m = vt.mean(); sd = vt.std(ddof=1) or 1.0
            tr.append((vt - m) / sd); hh.append((v[hi] - m) / sd)
        return np.array(tr).T.copy(), np.array(hh)
    for dname, mask, cols in (("P53", mP3, ["hpv", "p53_wt"]), ("E2F", mE3, ["e2f_g2m"])):
        u = mask.astype(np.float64); u -= u.mean(); u /= np.linalg.norm(u)
        y_tr = r_tr @ u; y_h = r_true @ u
        Ztr, zh = Z_of(cols)
        pred = zh @ ridge(Ztr, y_tr.T)
        sh = [zh @ ridge(Ztr[rng.permutation(len(ti))], y_tr.T) for _ in range(NPERM)]
        for i, t in enumerate(T82):
            row = {"held": h, "target": t, "direction": dname,
                   "correct": int(np.sign(pred[i]) == np.sign(y_h[i])), "y_true": float(y_h[i]), "y_pred": float(pred[i])}
            for p in range(NPERM): row[f"shuf{p}"] = int(np.sign(sh[p][i]) == np.sign(y_h[i]))
            gat.append(row)
    print("gating fold", h, flush=True)
gd = pd.DataFrame(gat); gd.to_csv(W / "gating289b.csv", index=False)

G = {}
for dname, sub in gd.groupby("direction"):
    acc = sub.groupby("target").correct.mean()
    nk = sub.groupby("target").correct.sum()
    shf = {p: sub.groupby("target")[f"shuf{p}"].sum() for p in range(NPERM)}
    frac9 = float((nk >= 9).mean())
    frac9_sh = float(np.mean([(shf[p] >= 9).mean() for p in range(NPERM)]))
    d = acc - sub.groupby("target")[[f"shuf{p}" for p in range(NPERM)]].mean().mean(1)
    tg = acc.index.to_numpy()
    b = np.array([d.loc[rng.choice(tg, len(tg), replace=True)].mean() for _ in range(4000)])
    G[dname] = {"pooled_accuracy": float(sub.correct.mean()),
                "pooled_accuracy_shuffled": float(sub[[f"shuf{p}" for p in range(NPERM)]].to_numpy().mean()),
                "frac_targets_ge9of12": frac9, "frac_targets_ge9of12_shuffled": frac9_sh,
                "binomial_expected_frac": float(stats.binom.sf(8, 12, 0.5)),
                "n_targets": int(len(acc)),
                "paired_vs_shuffled": [float(d.mean()), float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975)), float((b <= 0).mean())],
                "binom_p_pooled": float(stats.binomtest(int(sub.correct.sum()), int(len(sub)), 0.5, alternative="greater").pvalue)}
out["gating"] = G
print("GATING", json.dumps(G, indent=1), flush=True)
json.dump(out, open(W / "p1p2_289.json", "w"), indent=1, default=float)
