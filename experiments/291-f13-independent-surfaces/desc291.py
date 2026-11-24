"""291 — POST-HOC, labelled: written and run AFTER the two frozen verdicts were read.
Nothing here can change a verdict.

  A1  is `C_broad16` above its own wrong-target chance level, and where do the OR controls sit?
  A2  the most over- and under-predicted Broad-16 targets
  B1  is 290's DNA-binder exception present on the GWCD4i surface?
  B2  are the response and annotation models two views of the same thing (partials, combined)?

    PYTHONDONTWRITEBYTECODE=1 uv run python desc291.py   -> work/desc291.json
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

import lib291 as L

sys.path.insert(0, str(L.W268))
import load268 as B                                             # noqa: E402

T0 = time.time()
LOG: list[str] = []


def P(m=""):
    line = f"[{time.time() - T0:6.1f}s] {m}" if m else ""
    print(line, flush=True)
    LOG.append(line)


def main():
    h = L.check_prereg()
    P(f"prereg {h} MATCH (POST-HOC script)")
    out = {}
    A = pd.read_csv(L.WORK / "pred291_broad16.csv")
    Ad = A[~A["is_OR_control"]].reset_index(drop=True)
    Ao = A[A["is_OR_control"]].reset_index(drop=True)

    # ---- A1: the wrong-target chance level of the conservation statistic ------------------------
    H, genes = B.load_halves(B.WORKING)
    lines = B.WORKING
    tg = B.ALL_TARGETS
    R, ok, avail = B.build_responses(H, genes, lines, tg, keys=("full",), min_cpm=5.0)
    tpos = {t: i for i, t in enumerate(tg)}
    unit = list(Ad["target"])
    uidx = np.array([tpos[t] for t in unit])
    Ru = R["full"][uidx].astype(np.float64)
    U = Ru - Ru.mean(0)[None, :, :]
    U /= np.maximum(np.linalg.norm(U, axis=2), 1e-12)[:, :, None]
    pairs = [(i, j) for i in range(len(lines)) for j in range(i + 1, len(lines))]
    rng = np.random.default_rng(L.SEED)
    nulls = []
    for _ in range(20):                       # 20 independent wrong-target permutations
        perm = rng.permutation(len(unit))
        while (perm == np.arange(len(unit))).any():
            perm = rng.permutation(len(unit))
        nulls.append(np.stack([(U[:, i] * U[perm][:, j]).sum(1) for i, j in pairs], 1).mean(1))
    nulls = np.stack(nulls)
    obs = Ad["C_broad16"].to_numpy(float)
    a1 = dict(observed_median=float(np.median(obs)),
              null_median=float(np.median(nulls)), null_p95=float(np.percentile(nulls, 95)),
              frac_targets_above_null_p95=float((obs > np.percentile(nulls, 95)).mean()),
              or_median=float(np.median(Ao["C_broad16"])),
              mwu_p_targets_vs_or=float(stats.mannwhitneyu(obs, Ao["C_broad16"].to_numpy(float),
                                                           alternative="greater").pvalue),
              n_or=len(Ao))
    P(f"A1 C_broad16 median {a1['observed_median']:+.5f}; wrong-target null median "
      f"{a1['null_median']:+.5f} (p95 {a1['null_p95']:+.5f}); "
      f"{a1['frac_targets_above_null_p95']:.1%} of targets above the null p95; "
      f"OR controls median {a1['or_median']:+.5f}, MWU(targets > OR) p {a1['mwu_p_targets_vs_or']:.3f}")
    out["A1"] = a1
    del R, U, Ru, H

    # ---- A2: the extremes --------------------------------------------------------------------
    r_obs = stats.rankdata(obs) / len(obs)
    r_pred = stats.rankdata(Ad["pred_M_resp"].to_numpy(float)) / len(obs)
    resid = r_obs - r_pred
    o = np.argsort(resid)
    a2 = dict(most_over_predicted=[dict(target=unit[i], C=float(obs[i]),
                                        rank_resid=float(resid[i])) for i in o[:12]],
              most_under_predicted=[dict(target=unit[i], C=float(obs[i]),
                                         rank_resid=float(resid[i])) for i in o[::-1][:12]])
    P("A2 most UNDER-predicted (more conserved than the geometry says): "
      + ", ".join(x["target"] for x in a2["most_under_predicted"]))
    P("A2 most OVER-predicted: " + ", ".join(x["target"] for x in a2["most_over_predicted"]))
    out["A2"] = a2

    # ---- A3: partials between the two models on surface A ---------------------------------------
    a3 = dict(
        resp_given_ann=L.partial_spearman(Ad["pred_M_resp"], obs, Ad["pred_M_ann"]),
        ann_given_resp=L.partial_spearman(Ad["pred_M_ann"], obs, Ad["pred_M_resp"]),
        rho_resp_ann=L.spearman(Ad["pred_M_resp"], Ad["pred_M_ann"]))
    P(f"A3 partials: rho(M_resp, C | M_ann) {a3['resp_given_ann']:+.4f}; "
      f"rho(M_ann, C | M_resp) {a3['ann_given_resp']:+.4f}; "
      f"rho(M_resp, M_ann) {a3['rho_resp_ann']:+.4f}")
    out["A3"] = a3

    # ---- B1: the DNA-binder exception on the GWCD4i surface -------------------------------------
    Bt = pd.read_csv(L.WORK / "pred291_gwcd4i.csv")
    ann = L.load_annotation()
    tf = dict(zip(ann["gene"], (ann["is_tf_collectri"].fillna(0) > 0)
                  | (ann["is_tf_dorothea"].fillna(0) > 0)))
    b1 = {}
    for cond in L.CONDS:
        s = Bt[Bt["condition"] == cond].reset_index(drop=True)
        y = stats.rankdata(s["logfc_pearson_r"]) / len(s)
        p = stats.rankdata(s["pred_M_resp"]) / len(s)
        res = y - p
        znf = s["znf_class"].astype(str).str.startswith("ZNF").to_numpy()
        istf = np.array([bool(tf.get(g, False)) for g in s["gene"]])
        row = {}
        for lab, m in (("ZNF", znf), ("TF(CollecTRI|DoRothEA)", istf), ("other", ~(znf | istf))):
            row[lab] = dict(n=int(m.sum()), mean_rank_resid=float(res[m].mean()),
                            median_rank_resid=float(np.median(res[m])),
                            mwu_p_vs_rest=float(stats.mannwhitneyu(
                                res[m], res[~m], alternative="greater").pvalue) if 0 < m.sum() < len(m) else float("nan"))
        b1[cond] = row
        P(f"B1 {cond:9s} " + "; ".join(
            f"{k} n {v['n']} mean rank-resid {v['mean_rank_resid']:+.4f} (p {v['mwu_p_vs_rest']:.3f})"
            for k, v in row.items()))
    out["B1"] = b1

    # ---- B2: are the two models two views? partials + the combined fit --------------------------
    store = L.Store()
    rem = np.load(L.REM284, allow_pickle=True)
    fold = np.asarray(rem["fold"], int)
    comp = np.asarray(rem["comp"], int)
    z = np.load(L.HELD279, allow_pickle=True)
    c_held = np.asarray(z["c_held"], float)
    agg = store.aggregate(store.subset(exclude=L.K562_SOURCES))
    Xc = np.column_stack([agg[L.FEATS290].to_numpy(float), ann[L.FEATS].to_numpy(float)])
    pred_c, info_c = L.fit_apply_grouped(Xc, c_held, fold, Xc)
    gpos = {g: i for i, g in enumerate(store.genes)}
    b2 = {}
    for cond in L.CONDS:
        s = Bt[Bt["condition"] == cond].reset_index(drop=True)
        gi = np.array([gpos[g] for g in s["gene"]])
        y = s["logfc_pearson_r"].to_numpy(float)
        blk = comp[gi]
        preds = {"M_resp": s["pred_M_resp"].to_numpy(float),
                 "M_ann": s["pred_M_ann"].to_numpy(float),
                 "M_both": pred_c[gi]}
        bb = L.block_boot(y, preds, blk, np.random.default_rng(L.SEED))
        row = dict(
            rho_both=bb["M_both"], rho_resp=bb["M_resp"], rho_ann=bb["M_ann"],
            paired_both_minus_ann=bb["_paired_rho"]["M_both-M_ann"],
            paired_both_minus_resp=bb["_paired_rho"]["M_both-M_resp"],
            resp_given_ann=L.partial_spearman(preds["M_resp"], y, preds["M_ann"]),
            ann_given_resp=L.partial_spearman(preds["M_ann"], y, preds["M_resp"]),
            rho_resp_ann=L.spearman(preds["M_resp"], preds["M_ann"]))
        b2[cond] = row
        P(f"B2 {cond:9s} M_both {row['rho_both']['point']:+.4f} "
          f"[{row['rho_both']['lo']:+.4f}, {row['rho_both']['hi']:+.4f}]; "
          f"vs M_ann {row['paired_both_minus_ann']['mean']:+.4f} "
          f"[{row['paired_both_minus_ann']['lo']:+.4f}, {row['paired_both_minus_ann']['hi']:+.4f}]; "
          f"rho(M_resp,y|M_ann) {row['resp_given_ann']:+.4f}; "
          f"rho(M_ann,y|M_resp) {row['ann_given_resp']:+.4f}; "
          f"rho(M_resp,M_ann) {row['rho_resp_ann']:+.4f}")
    out["B2"] = dict(per_condition=b2, combined_model_info=info_c)

    (L.WORK / "desc291.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "desc291.log").write_text("\n".join(LOG))
    P("wrote work/desc291.json")


if __name__ == "__main__":
    main()
