"""281 — POST-HOC, run after the preregistered held numbers were read (fit281.json).

Nothing here can change the verdict. It does five things:

 (a) KA-4 CORRECTED — the shuffled-target null over MANY shuffle draws (the preregistered version
     took a gene bootstrap conditional on ONE draw, which is not the null distribution).
 (b) a CALIBRATED single-best-feature baseline, so its R^2 is on the same scale as the models'.
 (c) the reliability partial — how much of the held prediction survives removing `r_bar`.
 (d) feature-BLOCK ablation on the family-held split (each block alone / each block dropped).
 (e) out-of-fold predictions for all 17,106 held genes under grouped (family-block) 5-fold CV,
     the residual tail, and the family concentration of the residual.

Outputs: work/extra281.json, work/extra281.log, and extra columns in work/pred281.csv.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold

import fit281 as F                       # reuse the frozen feature list / helpers / seed

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SEED = F.SEED

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


BLOCKS = {
    "essentiality": ["dm_mean_effect", "dm_sd_effect", "dm_frac_dependent", "dm_frac_strong",
                     "dm_min_effect", "dm_common_essential", "replogle_essential"],
    "expression": ["ccle_breadth", "ccle_median_tpm", "ccle_sd"],
    "regulon": ["log_collectri_regulon_size", "collectri_indegree", "is_tf_collectri",
                "log_dorothea_regulon_size", "is_tf_dorothea"],
    "complexes": ["humap_n_complexes", "humap_max_complex_size", "reactome_n_complexes",
                  "reactome_max_complex_size"],
    "evolution": ["log_n_paralogues", "ortholog_depth"],
}


def main() -> None:
    out: dict = {"seed": SEED}
    df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    panel = df[df["set"] == "panel"].reset_index(drop=True)
    fit = json.loads((WORK / "fit281.json").read_text())

    block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                     held["gene_group1"].fillna("").to_numpy(),
                     "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
    held["block"] = block
    rng = np.random.default_rng(SEED)
    blocks = np.array(sorted(set(block)))
    order = rng.permutation(len(blocks))
    counts = held["block"].value_counts()
    cum, take = 0, []
    for k in blocks[order]:
        if cum >= F.HELD_FRAC * len(held):
            break
        take.append(k)
        cum += int(counts[k])
    fam_test = held["block"].isin(set(take)).to_numpy()
    y = held["y_raw"].to_numpy(float)
    X = held[F.FEATS].to_numpy(float)
    tr, te = ~fam_test, fam_test
    P(f"family-held split reproduced: {int(te.sum())} test genes "
      f"(fit281 recorded {fit['splits']['family']['n_genes_held']}) -> "
      f"{'MATCH' if int(te.sum()) == fit['splits']['family']['n_genes_held'] else 'MISMATCH'}")

    # -------------------------------------------------------------- (a) KA-4 corrected
    P("\n(a) KA-4 CORRECTED: held Spearman under 200 independent shuffles of the TRAINING target")
    obs_r = fit["results"]["family_raw"]["ridge"]["spearman"]
    obs_g = fit["results"]["family_raw"]["gbm"]["spearman"]
    nulls_r, nulls_g = [], []
    gbm_params = {**fit["results"]["family_raw"]["gbm"]["params"], "random_state": SEED}
    alpha = fit["results"]["family_raw"]["ridge"]["alpha"]
    st = F.Standardiser().fit(X[tr])
    Xtr_s, Xte_s = st.transform(X[tr]), st.transform(X[te])
    for d in range(200):
        ys = y[tr].copy()
        np.random.default_rng(SEED + 1000 + d).shuffle(ys)
        nulls_r.append(F.spearman(y[te], Ridge(alpha=alpha).fit(Xtr_s, ys).predict(Xte_s)))
        if d < 40:                                     # the GBM is the expensive one
            nulls_g.append(F.spearman(
                y[te], HistGradientBoostingRegressor(**gbm_params).fit(X[tr], ys).predict(X[te])))
    nr, ng = np.array(nulls_r), np.array(nulls_g)
    ka4 = dict(
        ridge=dict(observed=obs_r, null_mean=float(nr.mean()), null_sd=float(nr.std()),
                   null_p2_5=float(np.percentile(nr, 2.5)), null_p97_5=float(np.percentile(nr, 97.5)),
                   n_draws=len(nr), draws_ge_observed=int((nr >= obs_r).sum()),
                   z=float((obs_r - nr.mean()) / max(nr.std(), 1e-12))),
        gbm=dict(observed=obs_g, null_mean=float(ng.mean()), null_sd=float(ng.std()),
                 null_p2_5=float(np.percentile(ng, 2.5)), null_p97_5=float(np.percentile(ng, 97.5)),
                 n_draws=len(ng), draws_ge_observed=int((ng >= obs_g).sum()),
                 z=float((obs_g - ng.mean()) / max(ng.std(), 1e-12))))
    ka4["pass_"] = bool(ka4["ridge"]["draws_ge_observed"] == 0 and ka4["gbm"]["draws_ge_observed"] == 0
                        and abs(ka4["ridge"]["null_mean"]) < 0.02 and abs(ka4["gbm"]["null_mean"]) < 0.02)
    for k in ("ridge", "gbm"):
        v = ka4[k]
        P(f"  {k:5s} observed {v['observed']:+.4f}; null over {v['n_draws']} draws "
          f"{v['null_mean']:+.4f} ± {v['null_sd']:.4f} "
          f"[{v['null_p2_5']:+.4f}, {v['null_p97_5']:+.4f}]; draws >= observed "
          f"{v['draws_ge_observed']}; z {v['z']:+.1f}")
    P(f"  KA-4 (corrected) -> {'PASS' if ka4['pass_'] else 'FAIL'}")
    out["ka4_corrected"] = ka4

    # -------------------------------------------------------------- (b) calibrated best single
    P("\n(b) single best feature, CALIBRATED (univariate ridge fitted on the training genes)")
    rhos = np.array([F.spearman(X[tr][:, j], y[tr]) for j in range(X.shape[1])])
    j = int(np.nanargmax(np.abs(rhos)))
    name = F.FEATS[j]
    st1 = F.Standardiser().fit(X[tr][:, [j]])
    m1 = Ridge(alpha=1.0).fit(st1.transform(X[tr][:, [j]]), y[tr])
    p1 = m1.predict(st1.transform(X[te][:, [j]]))
    r = np.random.default_rng(SEED + 21)
    out["best_single_calibrated"] = dict(feature=name, train_spearman=float(rhos[j]),
                                         **F.report(y[te], p1, r))
    v = out["best_single_calibrated"]
    P(f"  {name}: rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]  "
      f"R2 {v['r2']:+.4f} [{v['r2_lo']:+.4f}, {v['r2_hi']:+.4f}]  AUC10 {v['auc_bottom_decile']:.4f}")

    # -------------------------------------------------------------- (c) the reliability partial
    P("\n(c) partial association after removing measured reliability `r_bar` (rank scale)")
    # fit281 stripped the raw predictions from the json; recompute the ridge prediction exactly
    # (same alpha, same standardiser, same training rows -> identical numbers).
    pred_te = Ridge(alpha=alpha).fit(Xtr_s, y[tr]).predict(Xte_s)
    chk = F.spearman(y[te], pred_te)
    P(f"  regression check: recomputed held ridge rho {chk:+.6f} vs fit281's "
      f"{fit['results']['family_raw']['ridge']['spearman']:+.6f} "
      f"(|Δ| {abs(chk - fit['results']['family_raw']['ridge']['spearman']):.2e})")
    rb = held["r_bar"].to_numpy(float)[te]
    ok = np.isfinite(rb)

    def rank(a):
        return stats.rankdata(a)

    def resid_on(a, b):
        A, B = rank(a), rank(b)
        B = (B - B.mean()) / B.std()
        return A - A.mean() - float((A * B).mean() / (B * B).mean()) * B

    yt, pt = y[te][ok], pred_te[ok]
    part = dict(
        n=int(ok.sum()),
        raw_rho=F.spearman(yt, pt),
        rho_r_bar_vs_y=F.spearman(yt, rb[ok]),
        rho_pred_vs_r_bar=F.spearman(pt, rb[ok]),
        partial_rho=float(np.corrcoef(resid_on(yt, rb[ok]), resid_on(pt, rb[ok]))[0, 1]))
    P(f"  n {part['n']}: rho(pred, y) {part['raw_rho']:+.4f}; rho(r_bar, y) "
      f"{part['rho_r_bar_vs_y']:+.4f}; rho(pred, r_bar) {part['rho_pred_vs_r_bar']:+.4f}; "
      f"PARTIAL rho(pred, y | r_bar) {part['partial_rho']:+.4f}")
    out["reliability_partial"] = part

    # -------------------------------------------------------------- (d) feature-block ablation
    P("\n(d) feature-block ablation, family-held split, ridge + gbm")
    abl = {}
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    for label, cols in ([(f"only::{b}", c) for b, c in BLOCKS.items()]
                        + [(f"drop::{b}", [f for f in F.FEATS if f not in c])
                           for b, c in BLOCKS.items()]
                        + [("all", list(F.FEATS))]):
        Xi = held[cols].to_numpy(float)
        pr_, mr_ = F.fit_ridge(Xi[tr], y[tr], Xi[te])
        pg_, mg_, _ = F.fit_gbm(Xi[tr], y[tr], Xi[te])
        abl[label] = dict(n_features=len(cols),
                          ridge_rho=F.spearman(y[te], pr_), ridge_r2=F.r2(y[te], pr_),
                          gbm_rho=F.spearman(y[te], pg_), gbm_r2=F.r2(y[te], pg_))
        a = abl[label]
        P(f"  {label:22s} ({a['n_features']:2d} feats)  ridge rho {a['ridge_rho']:+.4f} "
          f"R2 {a['ridge_r2']:+.4f}   gbm rho {a['gbm_rho']:+.4f} R2 {a['gbm_r2']:+.4f}")
    out["block_ablation"] = abl

    # -------------------------------------------------------------- (e) out-of-fold + residuals
    P("\n(e) grouped 5-fold out-of-fold predictions over all 17,106 held genes (blocks = families)")
    gkf = GroupKFold(n_splits=5)
    oof_r = np.full(len(held), np.nan)
    oof_g = np.full(len(held), np.nan)
    for k, (a, b) in enumerate(gkf.split(X, y, groups=block)):
        pr_, mr_ = F.fit_ridge(X[a], y[a], X[b])
        oof_r[b] = pr_
        mg = HistGradientBoostingRegressor(**gbm_params).fit(X[a], y[a])
        oof_g[b] = mg.predict(X[b])
        P(f"  fold {k}: train {len(a)}, test {len(b)}, alpha {mr_['alpha']:.3g}, "
          f"rho {F.spearman(y[b], pr_):+.4f} / {F.spearman(y[b], oof_g[b]):+.4f}")
    rr = np.random.default_rng(SEED + 31)
    out["oof"] = dict(ridge=F.report(y, oof_r, rr), gbm=F.report(y, oof_g, rr))
    P(f"  OOF overall: ridge rho {out['oof']['ridge']['spearman']:+.4f} "
      f"[{out['oof']['ridge']['spearman_lo']:+.4f}, {out['oof']['ridge']['spearman_hi']:+.4f}] "
      f"R2 {out['oof']['ridge']['r2']:+.4f} AUC10 {out['oof']['ridge']['auc_bottom_decile']:.4f}; "
      f"gbm rho {out['oof']['gbm']['spearman']:+.4f} R2 {out['oof']['gbm']['r2']:+.4f}")

    resid = stats.rankdata(y) - stats.rankdata(oof_r)          # rank residual, ridge
    resid = resid / len(y)                                      # in units of the gene fraction
    held["oof_ridge"] = oof_r
    held["oof_gbm"] = oof_g
    held["rank_resid"] = resid
    o = np.argsort(resid)
    cols = ["gene", "gene_group1", "y_raw", "oof_ridge", "rank_resid", "r_bar", "n_pairs"]
    worse = held.iloc[o[:20]][cols]                             # far BELOW prediction
    better = held.iloc[o[-20:]][cols]                           # far ABOVE prediction
    out["residual_tail"] = dict(
        below=worse.to_dict("records"), above=better.to_dict("records"))
    P("\n  20 genes FAR BELOW prediction (state-dependent given their biology):")
    for _, rrow in worse.iterrows():
        P(f"    {rrow['gene']:12s} C {rrow['y_raw']:+.4f}  pred {rrow['oof_ridge']:+.4f}  "
          f"rankres {rrow['rank_resid']:+.3f}  [{str(rrow['gene_group1'])[:40]}]")
    P("\n  20 genes FAR ABOVE prediction (conserved beyond their biology):")
    for _, rrow in better.iloc[::-1].iterrows():
        P(f"    {rrow['gene']:12s} C {rrow['y_raw']:+.4f}  pred {rrow['oof_ridge']:+.4f}  "
          f"rankres {rrow['rank_resid']:+.3f}  [{str(rrow['gene_group1'])[:40]}]")

    # family concentration of the residual
    fam = held.loc[held["gene_group1"].fillna("") != ""].groupby("gene_group1")
    tabf = fam["rank_resid"].agg(["mean", "median", "count"])
    tabf = tabf[tabf["count"] >= 30].sort_values("median")
    out["residual_by_family"] = dict(
        bottom=tabf.head(15).reset_index().to_dict("records"),
        top=tabf.tail(15).reset_index().to_dict("records"),
        n_families=int(len(tabf)))
    P(f"\n  residual by HGNC family ({len(tabf)} families with >= 30 held genes), "
      f"median rank residual:")
    P("    LOWEST (less conserved than their biology predicts):")
    for _, rrow in tabf.head(12).reset_index().iterrows():
        P(f"      {str(rrow['gene_group1'])[:46]:46s} n {int(rrow['count']):4d}  "
          f"median {rrow['median']:+.3f}  mean {rrow['mean']:+.3f}")
    P("    HIGHEST (more conserved than their biology predicts):")
    for _, rrow in tabf.tail(12).reset_index().iloc[::-1].iterrows():
        P(f"      {str(rrow['gene_group1'])[:46]:46s} n {int(rrow['count']):4d}  "
          f"median {rrow['median']:+.3f}  mean {rrow['mean']:+.3f}")

    # the ZNF class, post-hoc: does the model already explain 279's deficit?
    isz = held["znf_class"].to_numpy() == "ZNF"
    isn = held["znf_class"].to_numpy() == "NONZNF"
    znf = dict(
        d_median_observed=float(np.median(y[isz]) - np.median(y[isn])),
        d_median_predicted=float(np.median(oof_r[isz]) - np.median(oof_r[isn])),
        d_median_residual=float(np.median(resid[isz]) - np.median(resid[isn])))
    P(f"\n  279's ZNF deficit through this model: observed Δmedian "
      f"{znf['d_median_observed']:+.5f}; PREDICTED from gene features "
      f"{znf['d_median_predicted']:+.5f}; left in the rank residual "
      f"{znf['d_median_residual']:+.4f}")
    out["znf_through_the_model"] = znf

    # -------------------------------------------------------------- write
    pred = pd.read_csv(WORK / "pred281.csv")
    pred["pred_oof_family_ridge"] = oof_r
    pred["pred_oof_family_gbm"] = oof_g
    pred["rank_resid_oof_ridge"] = resid
    prow = panel[["gene", "y_raw", "y_rel", "r_bar", "gene_group1"]].copy()
    Xp = panel[F.FEATS].to_numpy(float)
    stf = F.Standardiser().fit(X[tr])
    prow["pred_family_ridge"] = Ridge(alpha=alpha).fit(stf.transform(X[tr]),
                                                       y[tr]).predict(stf.transform(Xp))
    prow["pred_family_gbm"] = HistGradientBoostingRegressor(**gbm_params).fit(
        X[tr], y[tr]).predict(Xp)
    prow = prow.rename(columns={"y_raw": "observed", "y_rel": "observed_rel"})
    prow["set"] = "panel"
    prow["split_family"] = "held(second set)"
    prow["split_random"] = "held(second set)"
    pred = pd.concat([pred, prow], ignore_index=True)
    pred.to_csv(WORK / "pred281.csv", index=False)
    P(f"\nwrote work/pred281.csv ({len(pred)} rows: 17,106 held genes + 300 panel targets)")

    out["log"] = LOG
    (WORK / "extra281.json").write_text(json.dumps(out, indent=1, default=float))
    P("wrote work/extra281.json")


if __name__ == "__main__":
    main()
