"""281 — the measurement-design control, POST-HOC and decisive.

`C_held` is built from however many cross-line pairs happen to measure a gene: 5,420 of the 17,106
held genes have 1-2 pairs (sd of C 0.126) and 2,438 have >= 17 (sd 0.072); genes with 3-16 pairs
have sd ~0.04. Features that track "is this gene expressed / essential" also track how many sources
measured it, so part of the preregistered result can be prediction of the MAP'S NOISE rather than of
conservation.

279's `c3_balanced_core` removes that confound by construction: the conservation of every gene
computed on EXACTLY the same three cross-line pairs {K562-GWPS x Orion HCT116, K562-GWPS x Orion
HEK293T, Orion HCT116 x Orion HEK293T}, on the 8,541 held genes present in all three sources
(279 §4 c3). This script repeats the preregistered analysis on that target, and stratifies the
primary target by pair count.

Outputs: work/core281.json, work/core281.log
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
from sklearn.model_selection import GroupKFold

import boot281 as B
import fit281 as F

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


def main() -> None:
    out: dict = {"seed": SEED}
    df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    z = pd.read_csv(ROOT / "experiments/279-znf-state-dependence/work/znf279.csv",
                    usecols=["gene", "c3_balanced_core", "c1_orion_pair"])
    held = held.merge(z, on="gene", how="left")

    # KA-5: 279 §4(c3) reports the balanced core's ZNF contrast; reproduce it off the column.
    isz = held["znf_class"].to_numpy() == "ZNF"
    isn = held["znf_class"].to_numpy() == "NONZNF"
    c3 = held["c3_balanced_core"].to_numpy(float)
    ok3 = np.isfinite(c3)
    ka5 = dict(n=int(ok3.sum()),
               znf_median=float(np.median(c3[isz & ok3])),
               nonznf_median=float(np.median(c3[isn & ok3])))
    ka5["d_median"] = ka5["znf_median"] - ka5["nonznf_median"]
    ka5["pass_"] = bool(abs(ka5["znf_median"] - 0.01181) < 5e-5
                        and abs(ka5["nonznf_median"] - 0.02966) < 5e-5
                        and abs(ka5["d_median"] + 0.01785) < 5e-5)
    P(f"KA-5 (279 §4c3 balanced core): n {ka5['n']}, ZNF median {ka5['znf_median']:+.5f} "
      f"(req +0.01181), non-ZNF {ka5['nonznf_median']:+.5f} (req +0.02966), Δ "
      f"{ka5['d_median']:+.5f} (req -0.01785) -> {'PASS' if ka5['pass_'] else 'FAIL'}")
    out["ka5"] = ka5

    block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                     held["gene_group1"].fillna("").to_numpy(),
                     "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
    rng = np.random.default_rng(SEED)
    blocks = np.array(sorted(set(block)))
    order = rng.permutation(len(blocks))
    counts = pd.Series(block).value_counts()
    cum, take = 0, []
    for k in blocks[order]:
        if cum >= F.HELD_FRAC * len(held):
            break
        take.append(k)
        cum += int(counts[k])
    te_all = pd.Series(block).isin(set(take)).to_numpy()

    X = held[F.FEATS].to_numpy(float)
    XR = held[F.FEATS + F.RELIABILITY].to_numpy(float)
    XO = held[F.RELIABILITY].to_numpy(float)

    def arm(y, mask_ok, label):
        te = te_all & mask_ok
        tr = (~te_all) & mask_ok
        P(f"\n=== {label}: train {int(tr.sum())}, held-family test {int(te.sum())} "
          f"({len(set(block[te]))} blocks) ===")

        def R(M):
            p, m = F.fit_ridge(M[tr], y[tr], M[te])
            return p, m

        def G(M):
            p, m, _ = F.fit_gbm(M[tr], y[tr], M[te])
            return p, m

        pr, mr = R(X)
        pg, mg = G(X)
        po, mo = R(XO)
        pog, mog = G(XO)
        pp, mp = R(XR)
        ppg, mpg = G(XR)
        j = F.FEATS.index("replogle_essential")
        pb, _ = R(X[:, [j]])
        preds = {"features_ridge": pr, "features_gbm": pg, "reliability_only": po,
                 "reliability_only_gbm": pog, "features_plus_reliability": pp,
                 "features_plus_reliability_gbm": ppg, "best_single_calibrated": pb,
                 "mean": np.full(int(te.sum()), float(y[tr].mean()))}
        res = B.block_boot(y[te], preds, block[te], np.random.default_rng(SEED + 201), n=2000)
        for k in [k for k in res if not k.startswith("_")]:
            v = res[k]
            P(f"  {k:30s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, "
              f"{v['rho']['hi']:+.4f}]   R2 {v['r2']['point']:+.4f} [{v['r2']['lo']:+.4f}, "
              f"{v['r2']['hi']:+.4f}]   AUC10 {v['auc']['point']:.4f} "
              f"[{v['auc']['lo']:.4f}, {v['auc']['hi']:.4f}]")
        P("  paired margins:")
        for k in ("features_ridge-reliability_only",
                  "features_plus_reliability-reliability_only",
                  "features_plus_reliability_gbm-reliability_only_gbm",
                  "features_ridge-best_single_calibrated"):
            v = res["_paired_rho"][k]
            P(f"    Δrho {k:52s} {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  "
              f"P(Δ<=0) {v['p_le_0']:.3f}")
        res["_hyper"] = dict(ridge_alpha=mr["alpha"], gbm=mg["params"])
        res["_n"] = dict(train=int(tr.sum()), test=int(te.sum()),
                         blocks=int(len(set(block[te]))))
        return res

    # -------------------------------------------------------- 1. the balanced core
    out["balanced_core"] = arm(c3, ok3, "BALANCED CORE c3 (identical 3 pairs for every gene)")

    # -------------------------------------------------------- 2. y_raw stratified by pair count
    P("\n=== the primary target stratified by pair count (why the core matters) ===")
    y = held["y_raw"].to_numpy(float)
    npairs = held["n_pairs"].to_numpy(float)
    te = te_all
    tr = ~te_all
    pr, mr = F.fit_ridge(X[tr], y[tr], X[te])
    strat = {}
    for lab, lo, hi in (("1-2", 1, 2), ("3-8", 3, 8), ("9-16", 9, 16), ("17+", 17, 1e9)):
        m = (npairs[te] >= lo) & (npairs[te] <= hi)
        if m.sum() < 50:
            continue
        strat[lab] = dict(n=int(m.sum()), sd_y=float(np.std(y[te][m])),
                          rho=F.spearman(y[te][m], pr[m]),
                          auc=F.auc_bottom_decile(y[te][m], pr[m]))
        P(f"  n_pairs {lab:5s}: n {strat[lab]['n']:5d}  sd(C) {strat[lab]['sd_y']:.4f}  "
          f"rho {strat[lab]['rho']:+.4f}  AUC10 {strat[lab]['auc']:.4f}")
    P(f"  spearman(C_held, n_pairs) over all held genes "
      f"{F.spearman(y, npairs):+.4f}; the ridge prediction's own "
      f"spearman with n_pairs on the test genes {F.spearman(pr, npairs[te]):+.4f}")
    out["stratified_by_pairs"] = strat
    out["rho_C_vs_npairs"] = F.spearman(y, npairs)
    out["rho_pred_vs_npairs"] = F.spearman(pr, npairs[te])

    # -------------------------------------------------------- 3. OOF on the core + residual tail
    P("\n=== grouped 5-fold OOF on the balanced core, and the residual tail ===")
    idx = np.flatnonzero(ok3)
    Xc, yc, bc = X[idx], c3[idx], block[idx]
    oof = np.full(len(idx), np.nan)
    gkf = GroupKFold(n_splits=5)
    for k, (a, b) in enumerate(gkf.split(Xc, yc, groups=bc)):
        p, m = F.fit_ridge(Xc[a], yc[a], Xc[b])
        oof[b] = p
        P(f"  fold {k}: test {len(b)}, alpha {m['alpha']:.3g}, rho {F.spearman(yc[b], p):+.4f}")
    rr = np.random.default_rng(SEED + 211)
    out["core_oof"] = F.report(yc, oof, rr)
    v = out["core_oof"]
    P(f"  OOF overall: rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
      f"  R2 {v['r2']:+.4f}  AUC10 {v['auc_bottom_decile']:.4f}  n {v['n']}")

    sub = held.iloc[idx].copy()
    sub["oof"] = oof
    sub["rank_resid"] = (stats.rankdata(yc) - stats.rankdata(oof)) / len(yc)
    o = np.argsort(sub["rank_resid"].to_numpy())
    cols = ["gene", "gene_group1", "c3_balanced_core", "oof", "rank_resid", "r_bar", "n_pairs"]
    below = sub.iloc[o[:20]][cols]
    above = sub.iloc[o[-20:]][cols]
    out["core_residual_tail"] = dict(below=below.to_dict("records"),
                                     above=above.to_dict("records"))
    P("\n  20 genes FAR BELOW prediction on the balanced core "
      "(state-dependent given their biology):")
    for _, r_ in below.iterrows():
        P(f"    {r_['gene']:12s} C3 {r_['c3_balanced_core']:+.4f}  pred {r_['oof']:+.4f}  "
          f"rankres {r_['rank_resid']:+.3f}  r {r_['r_bar']:+.3f}  [{str(r_['gene_group1'])[:38]}]")
    P("\n  20 genes FAR ABOVE prediction on the balanced core:")
    for _, r_ in above.iloc[::-1].iterrows():
        P(f"    {r_['gene']:12s} C3 {r_['c3_balanced_core']:+.4f}  pred {r_['oof']:+.4f}  "
          f"rankres {r_['rank_resid']:+.3f}  r {r_['r_bar']:+.3f}  [{str(r_['gene_group1'])[:38]}]")

    fam = sub.loc[sub["gene_group1"].fillna("") != ""].groupby("gene_group1")
    tabf = fam["rank_resid"].agg(["mean", "median", "count"])
    tabf = tabf[tabf["count"] >= 20].sort_values("median")
    out["core_residual_by_family"] = dict(
        n_families=int(len(tabf)),
        bottom=tabf.head(15).reset_index().to_dict("records"),
        top=tabf.tail(15).reset_index().to_dict("records"))
    P(f"\n  residual by HGNC family on the core ({len(tabf)} families with >= 20 genes):")
    P("    LOWEST (less conserved than their biology predicts):")
    for _, r_ in tabf.head(12).reset_index().iterrows():
        P(f"      {str(r_['gene_group1'])[:48]:48s} n {int(r_['count']):4d}  "
          f"median {r_['median']:+.3f}")
    P("    HIGHEST (more conserved than their biology predicts):")
    for _, r_ in tabf.tail(12).reset_index().iloc[::-1].iterrows():
        P(f"      {str(r_['gene_group1'])[:48]:48s} n {int(r_['count']):4d}  "
          f"median {r_['median']:+.3f}")

    znf = dict(observed=float(np.median(yc[held["znf_class"].to_numpy()[idx] == "ZNF"])
                              - np.median(yc[held["znf_class"].to_numpy()[idx] == "NONZNF"])),
               predicted=float(np.median(oof[held["znf_class"].to_numpy()[idx] == "ZNF"])
                               - np.median(oof[held["znf_class"].to_numpy()[idx] == "NONZNF"])),
               residual=float(np.median(sub["rank_resid"].to_numpy()[
                   held["znf_class"].to_numpy()[idx] == "ZNF"])
                   - np.median(sub["rank_resid"].to_numpy()[
                       held["znf_class"].to_numpy()[idx] == "NONZNF"])))
    P(f"\n  279's ZNF deficit on the balanced core: observed Δmedian {znf['observed']:+.5f}; "
      f"predicted from gene features {znf['predicted']:+.5f}; left in the rank residual "
      f"{znf['residual']:+.4f}")
    out["znf_on_core"] = znf

    out["log"] = LOG
    (WORK / "core281.json").write_text(json.dumps(out, indent=1, default=float))
    P("\nwrote work/core281.json")


if __name__ == "__main__":
    main()
