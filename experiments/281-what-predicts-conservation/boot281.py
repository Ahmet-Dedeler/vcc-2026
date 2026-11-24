"""281 — honest intervals: a BLOCK (family) bootstrap of the held numbers, and the paired margin
over the reliability baseline.

Why: under a family-held split the test genes are not independent — they arrive in whole HGNC
families. The gene bootstrap in fit281 therefore understates the uncertainty (extra281's corrected
KA-4 measures the same thing from the other side: the shuffled-target null has SD 0.098 on ridge,
far larger than the gene bootstrap's ±0.031). Everything here resamples BLOCKS with replacement.

Outputs: work/boot281.json, work/boot281.log
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge

import fit281 as F

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SEED = F.SEED
NB = 2000

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def block_boot(y, preds: dict, groups, rng, n=NB):
    """Resample whole blocks with replacement; return per-predictor CIs and paired deltas."""
    uniq = np.array(sorted(set(groups)))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    keys = list(preds)
    acc = {k: {"rho": [], "r2": [], "auc": []} for k in keys}
    pair = {f"{a}-{b}": [] for a in keys for b in keys if a != b}
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        yy = y[sel]
        cur = {}
        for k in keys:
            pp = preds[k][sel]
            cur[k] = F.spearman(yy, pp)
            acc[k]["rho"].append(cur[k])
            acc[k]["r2"].append(F.r2(yy, pp))
            acc[k]["auc"].append(F.auc_bottom_decile(yy, pp))
        for a in keys:
            for b in keys:
                if a != b:
                    pair[f"{a}-{b}"].append(cur[a] - cur[b])
    out = {}
    for k in keys:
        out[k] = {m: dict(point=float(F.spearman(y, preds[k])) if m == "rho" else
                          (float(F.r2(y, preds[k])) if m == "r2"
                           else float(F.auc_bottom_decile(y, preds[k]))),
                          lo=float(np.nanpercentile(acc[k][m], 2.5)),
                          hi=float(np.nanpercentile(acc[k][m], 97.5)),
                          sd=float(np.nanstd(acc[k][m])))
                  for m in ("rho", "r2", "auc")}
    out["_paired_rho"] = {k: dict(mean=float(np.nanmean(v)),
                                  lo=float(np.nanpercentile(v, 2.5)),
                                  hi=float(np.nanpercentile(v, 97.5)),
                                  p_le_0=float(np.mean(np.asarray(v) <= 0)))
                          for k, v in pair.items()}
    return out


def main() -> None:
    df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    panel = df[df["set"] == "panel"].reset_index(drop=True)
    fit = json.loads((WORK / "fit281.json").read_text())

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
    te = pd.Series(block).isin(set(take)).to_numpy()
    tr = ~te
    y = held["y_raw"].to_numpy(float)
    X = held[F.FEATS].to_numpy(float)
    XR = held[F.FEATS + F.RELIABILITY].to_numpy(float)
    XO = held[F.RELIABILITY].to_numpy(float)
    alpha = fit["results"]["family_raw"]["ridge"]["alpha"]
    gbm_params = {**fit["results"]["family_raw"]["gbm"]["params"], "random_state": SEED}
    P(f"family-held: {int(te.sum())} test genes in {len(set(block[te]))} blocks; "
      f"alpha {alpha:.4g}, gbm {gbm_params}")

    def ridge_pred(M, a):
        st = F.Standardiser().fit(M[tr])
        return Ridge(alpha=a).fit(st.transform(M[tr]), y[tr]).predict(st.transform(M[te]))

    # the genes on which every arm is comparable: r_bar must exist
    okr = np.isfinite(held["r_bar"].to_numpy(float))
    P(f"r_bar present on {int(okr[te].sum())}/{int(te.sum())} held-family test genes")

    preds = {
        "features_ridge": ridge_pred(X, alpha),
        "features_gbm": HistGradientBoostingRegressor(**gbm_params).fit(X[tr], y[tr]).predict(X[te]),
        "reliability_only": ridge_pred(XO, fit["results"]["family_raw_rel_only"]["ridge"]["alpha"]),
        "features_plus_reliability":
            ridge_pred(XR, fit["results"]["family_raw_feats_plus_rel"]["ridge"]["alpha"]),
        "best_single_calibrated": ridge_pred(X[:, [F.FEATS.index("replogle_essential")]], 1.0),
        "mean": np.full(int(te.sum()), float(y[tr].mean())),
    }
    # a gbm arm on features + reliability, so the margin is also read on the stronger model
    preds["features_plus_reliability_gbm"] = HistGradientBoostingRegressor(
        **{**fit["results"]["family_raw_feats_plus_rel"]["gbm"]["params"],
           "random_state": SEED}).fit(XR[tr], y[tr]).predict(XR[te])
    preds["reliability_only_gbm"] = HistGradientBoostingRegressor(
        **{**fit["results"]["family_raw_rel_only"]["gbm"]["params"],
           "random_state": SEED}).fit(XO[tr], y[tr]).predict(XO[te])

    # regression check against fit281
    for k, ref in (("features_ridge", fit["results"]["family_raw"]["ridge"]["spearman"]),
                   ("features_gbm", fit["results"]["family_raw"]["gbm"]["spearman"]),
                   ("reliability_only", fit["results"]["family_raw_rel_only"]["ridge"]["spearman"]),
                   ("features_plus_reliability",
                    fit["results"]["family_raw_feats_plus_rel"]["ridge"]["spearman"])):
        got = F.spearman(y[te], preds[k])
        P(f"  regression check {k:28s} {got:+.6f} vs fit281 {ref:+.6f} |Δ| {abs(got - ref):.2e}")

    rb = np.random.default_rng(SEED + 101)
    res = block_boot(y[te], preds, block[te], rb)
    P("\nBLOCK (family) bootstrap, 2,000 draws over the held families:")
    for k in [k for k in res if not k.startswith("_")]:
        v = res[k]
        P(f"  {k:30s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, {v['rho']['hi']:+.4f}] "
          f"(sd {v['rho']['sd']:.4f})   R2 {v['r2']['point']:+.4f} "
          f"[{v['r2']['lo']:+.4f}, {v['r2']['hi']:+.4f}]   AUC10 {v['auc']['point']:.4f} "
          f"[{v['auc']['lo']:.4f}, {v['auc']['hi']:.4f}]")
    P("\nPaired margins (block bootstrap, same resample):")
    for k in ("features_ridge-reliability_only", "features_plus_reliability-reliability_only",
              "features_plus_reliability_gbm-reliability_only_gbm",
              "features_ridge-best_single_calibrated"):
        v = res["_paired_rho"][k]
        P(f"  Δrho {k:52s} {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  "
          f"P(Δ<=0) {v['p_le_0']:.3f}")

    # ---------------- the panel, same treatment (families as blocks) --------------------------
    Xp = panel[F.FEATS].to_numpy(float)
    stf = F.Standardiser().fit(X[tr])
    pp = {
        "features_ridge": Ridge(alpha=alpha).fit(stf.transform(X[tr]), y[tr]).predict(
            stf.transform(Xp)),
        "features_gbm": HistGradientBoostingRegressor(**gbm_params).fit(X[tr], y[tr]).predict(Xp),
        "mean": np.full(len(panel), float(y[tr].mean())),
    }
    pblock = np.where(panel["gene_group1"].fillna("").to_numpy() != "",
                      panel["gene_group1"].fillna("").to_numpy(),
                      "NOFAM:" + panel["gene"].astype(str).to_numpy())
    yp = panel["y_raw"].to_numpy(float)
    rbp = panel["r_bar"].to_numpy(float)
    okp = np.isfinite(yp)
    pres = block_boot(yp[okp], {k: v[okp] for k, v in pp.items()}, pblock[okp],
                      np.random.default_rng(SEED + 103), n=2000)
    P(f"\nPANEL (second held set, n {int(okp.sum())}), block bootstrap over "
      f"{len(set(pblock[okp]))} family blocks:")
    for k in [k for k in pres if not k.startswith("_")]:
        v = pres[k]
        P(f"  {k:20s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, {v['rho']['hi']:+.4f}] "
          f"  R2 {v['r2']['point']:+.4f} [{v['r2']['lo']:+.4f}, {v['r2']['hi']:+.4f}]   "
          f"AUC10 {v['auc']['point']:.4f} [{v['auc']['lo']:.4f}, {v['auc']['hi']:.4f}]")
    # the panel's own reliability baseline (271 ref_rel, hek293t reserved)
    okp2 = okp & np.isfinite(rbp)
    pres2 = block_boot(yp[okp2], {"features_ridge": pp["features_ridge"][okp2],
                                  "reliability_only": rbp[okp2]},
                       pblock[okp2], np.random.default_rng(SEED + 105), n=2000)
    P(f"\nPANEL vs its own reliability (271 `ref_rel`, k562_gwps+orion_hct116; n {int(okp2.sum())}):")
    for k in [k for k in pres2 if not k.startswith("_")]:
        v = pres2[k]
        P(f"  {k:20s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, {v['rho']['hi']:+.4f}]")
    v = pres2["_paired_rho"]["features_ridge-reliability_only"]
    P(f"  Δrho features - reliability  {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  "
      f"P(Δ<=0) {v['p_le_0']:.3f}")

    out = dict(seed=SEED, n_boot=NB, held=res, panel=pres, panel_vs_reliability=pres2, log=LOG)
    (WORK / "boot281.json").write_text(json.dumps(out, indent=1, default=float))
    P("wrote work/boot281.json")


if __name__ == "__main__":
    main()
