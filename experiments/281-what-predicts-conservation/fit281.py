"""281 — fit the preregistered models and report ONLY held-split numbers.

Run AFTER `PREREGISTRATION.md` is hashed. Everything in §"frozen" below is quoted from it.

  uv run python fit281.py

Outputs: work/fit281.json, work/pred281.csv, work/fit281.log
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

SEED = 20260906
N_BOOT = 2000
HELD_FRAC = 0.20
RNG = np.random.default_rng(SEED)

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


# ---- frozen feature set (PREREGISTRATION §3) --------------------------------------------------
FEATS = [
    # essentiality / fitness (DepMap 24Q2 Public, Chronos gene effect)
    "dm_mean_effect", "dm_sd_effect", "dm_frac_dependent", "dm_frac_strong", "dm_min_effect",
    "dm_common_essential", "replogle_essential",
    # expression (CCLE 2025, 1,600 lines)
    "ccle_breadth", "ccle_median_tpm", "ccle_sd",
    # regulatory role (CollecTRI / DoRothEA A-C)
    "log_collectri_regulon_size", "collectri_indegree", "is_tf_collectri",
    "log_dorothea_regulon_size", "is_tf_dorothea",
    # protein complexes (hu.MAP 2.0, Reactome)
    "humap_n_complexes", "humap_max_complex_size", "reactome_n_complexes",
    "reactome_max_complex_size",
    # evolution (Ensembl BioMart)
    "log_n_paralogues", "ortholog_depth",
]
NUISANCE = ["log_n_pairs", "n_lines"]        # declared variant only
RELIABILITY = ["r_bar"]                      # the confounder baseline


# ---- metrics ----------------------------------------------------------------------------------
def spearman(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[m], b[m]).statistic)


def r2(y: np.ndarray, p: np.ndarray) -> float:
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    sst = float(((y - y.mean()) ** 2).sum())
    return float(1.0 - ((y - p) ** 2).sum() / sst) if sst > 0 else float("nan")


def auc_bottom_decile(y: np.ndarray, p: np.ndarray) -> float:
    """P(a bottom-decile gene is predicted lower than a non-bottom-decile gene)."""
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    if len(y) < 20:
        return float("nan")
    thr = np.quantile(y, 0.10)
    lab = y <= thr
    if lab.sum() == 0 or (~lab).sum() == 0:
        return float("nan")
    r = stats.rankdata(-p)                     # low prediction -> high rank score
    n1, n0 = int(lab.sum()), int((~lab).sum())
    return float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def boot_ci(fn, y: np.ndarray, p: np.ndarray, rng: np.random.Generator, n: int = N_BOOT):
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    if len(y) < 20:
        return float("nan"), float("nan"), float("nan")
    point = fn(y, p)
    idx = rng.integers(0, len(y), size=(n, len(y)))
    vals = np.array([fn(y[i], p[i]) for i in idx])
    return point, float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def report(y: np.ndarray, p: np.ndarray, rng: np.random.Generator) -> dict:
    rho, lo, hi = boot_ci(lambda a, b: spearman(a, b), y, p, rng)
    rr, rlo, rhi = boot_ci(lambda a, b: r2(a, b), y, p, rng)
    au, alo, ahi = boot_ci(lambda a, b: auc_bottom_decile(a, b), y, p, rng)
    return dict(n=int((np.isfinite(y) & np.isfinite(p)).sum()),
                spearman=rho, spearman_lo=lo, spearman_hi=hi,
                r2=rr, r2_lo=rlo, r2_hi=rhi,
                auc_bottom_decile=au, auc_lo=alo, auc_hi=ahi)


# ---- models -----------------------------------------------------------------------------------
ALPHAS = np.logspace(-3, 4, 15)
GBM_GRID = [dict(learning_rate=lr, max_iter=it, max_depth=3, min_samples_leaf=40,
                 l2_regularization=1.0, early_stopping=False, random_state=SEED)
            for lr in (0.03, 0.10) for it in (200, 500)]


class Standardiser:
    def fit(self, X):
        self.med = np.nanmedian(X, axis=0)
        Xi = np.where(np.isfinite(X), X, self.med)
        self.mu = Xi.mean(axis=0)
        self.sd = Xi.std(axis=0)
        self.sd[self.sd < 1e-12] = 1.0
        return self

    def transform(self, X):
        Xi = np.where(np.isfinite(X), X, self.med)
        return (Xi - self.mu) / self.sd


def fit_ridge(Xtr, ytr, Xte):
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    best, best_a = -np.inf, ALPHAS[0]
    for a in ALPHAS:
        sc = []
        for tr, va in kf.split(Xtr):
            st = Standardiser().fit(Xtr[tr])
            m = Ridge(alpha=a).fit(st.transform(Xtr[tr]), ytr[tr])
            sc.append(spearman(ytr[va], m.predict(st.transform(Xtr[va]))))
        s = float(np.mean(sc))
        if s > best:
            best, best_a = s, a
    st = Standardiser().fit(Xtr)
    m = Ridge(alpha=best_a).fit(st.transform(Xtr), ytr)
    return m.predict(st.transform(Xte)), dict(alpha=float(best_a), inner_spearman=best,
                                              coef=m.coef_.tolist())


def fit_gbm(Xtr, ytr, Xte):
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    best, best_p = -np.inf, GBM_GRID[0]
    for prm in GBM_GRID:
        sc = []
        for tr, va in kf.split(Xtr):
            m = HistGradientBoostingRegressor(**prm).fit(Xtr[tr], ytr[tr])
            sc.append(spearman(ytr[va], m.predict(Xtr[va])))
        s = float(np.mean(sc))
        if s > best:
            best, best_p = s, prm
    m = HistGradientBoostingRegressor(**best_p).fit(Xtr, ytr)
    return m.predict(Xte), dict(params={k: v for k, v in best_p.items() if k != "random_state"},
                                inner_spearman=best), m


def fit_best_single(Xtr, ytr, Xte, names):
    rhos = np.array([spearman(Xtr[:, j], ytr) for j in range(Xtr.shape[1])])
    j = int(np.nanargmax(np.abs(rhos)))
    sgn = float(np.sign(rhos[j])) or 1.0
    col = Xte[:, j].astype(float)
    col = np.where(np.isfinite(col), col, np.nanmedian(Xtr[:, j]))
    return sgn * col, dict(feature=names[j], train_spearman=float(rhos[j]))


# ---- data -------------------------------------------------------------------------------------
def main() -> None:
    prereg = HERE / "PREREGISTRATION.md"
    h = hashlib.sha256(prereg.read_bytes()).hexdigest()
    frozen = (WORK / "PREREGISTRATION.sha256").read_text().split()[0]
    P(f"PREREGISTRATION.md sha256 {h}")
    P(f"frozen digest            {frozen}  -> {'MATCH' if h == frozen else 'MISMATCH'}")
    if h != frozen:
        sys.exit("preregistration digest mismatch — aborting")

    df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    panel = df[df["set"] == "panel"].reset_index(drop=True)
    P(f"held {len(held)} genes, panel {len(panel)} targets")

    out: dict = {"seed": SEED, "prereg_sha256": h, "feats": FEATS, "ka": {}, "splits": {},
                 "results": {}}

    # ---- splits (PREREGISTRATION §2) ----------------------------------------------------------
    block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                     held["gene_group1"].fillna("").to_numpy(),
                     "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
    held["block"] = block
    rng = np.random.default_rng(SEED)
    blocks = np.array(sorted(set(block)))
    order = rng.permutation(len(blocks))
    cum, take = 0, []
    counts = held["block"].value_counts()
    for k in blocks[order]:
        if cum >= HELD_FRAC * len(held):
            break
        take.append(k)
        cum += int(counts[k])
    fam_test = held["block"].isin(set(take)).to_numpy()
    P(f"family-held split: {len(take)} of {len(blocks)} blocks held, "
      f"{int(fam_test.sum())} genes ({fam_test.mean():.3f})")
    out["splits"]["family"] = dict(n_blocks=int(len(blocks)), n_blocks_held=len(take),
                                   n_genes_held=int(fam_test.sum()),
                                   frac=float(fam_test.mean()))

    rng2 = np.random.default_rng(SEED + 1)
    perm = rng2.permutation(len(held))
    rnd_test = np.zeros(len(held), dtype=bool)
    rnd_test[perm[: int(round(HELD_FRAC * len(held)))]] = True
    P(f"random-gene split: {int(rnd_test.sum())} genes held ({rnd_test.mean():.3f})")
    out["splits"]["random"] = dict(n_genes_held=int(rnd_test.sum()), frac=float(rnd_test.mean()))

    # family one-hot: the ~N HGNC families with >= 30 held genes (random split only)
    fam_counts = held.loc[held["gene_group1"] != "", "gene_group1"].value_counts()
    big_fams = sorted(fam_counts[fam_counts >= 30].index)
    P(f"family one-hot: {len(big_fams)} HGNC families with >= 30 held genes")
    out["splits"]["n_big_families"] = len(big_fams)
    FAM_COLS = [f"fam::{f}" for f in big_fams]
    for f in big_fams:
        held[f"fam::{f}"] = (held["gene_group1"] == f).astype(float)
        panel[f"fam::{f}"] = (panel["gene_group1"] == f).astype(float)

    # ---- KA-4: the shuffled-target control must read at chance -------------------------------
    def run(name: str, feats: list[str], test_mask: np.ndarray, ytarget: str,
            shuffle_y: bool = False, extra_eval: list[tuple[str, pd.DataFrame]] | None = None,
            rng_seed: int = SEED) -> dict:
        y = held[ytarget].to_numpy(float)
        ok = np.isfinite(y)
        tr = ok & ~test_mask
        te = ok & test_mask
        X = held[feats].to_numpy(float)
        ytr = y[tr].copy()
        if shuffle_y:
            np.random.default_rng(rng_seed).shuffle(ytr)
        res: dict = {"n_train": int(tr.sum()), "n_test": int(te.sum()), "features": feats}
        r = np.random.default_rng(rng_seed + 7)

        p_ridge, meta_r = fit_ridge(X[tr], ytr, X[te])
        p_gbm, meta_g, gbm_model = fit_gbm(X[tr], ytr, X[te])
        p_best, meta_b = fit_best_single(X[tr], ytr, X[te], feats)
        p_mean = np.full(int(te.sum()), float(ytr.mean()))

        res["ridge"] = {**report(y[te], p_ridge, r), **meta_r}
        res["gbm"] = {**report(y[te], p_gbm, r), **meta_g}
        res["best_single"] = {**report(y[te], p_best, r), **meta_b}
        res["mean"] = report(y[te], p_mean, r)

        # permutation importance on the HELD split only (GBM and ridge)
        if not shuffle_y:
            imp = {}
            base_rho = spearman(y[te], p_gbm)
            base_r2 = r2(y[te], p_gbm)
            Xte = X[te].copy()
            pr = np.random.default_rng(rng_seed + 13)
            for j, f in enumerate(feats):
                drho, dr2 = [], []
                for _ in range(10):
                    Xp = Xte.copy()
                    pr.shuffle(Xp[:, j])
                    pp = gbm_model.predict(Xp)
                    drho.append(base_rho - spearman(y[te], pp))
                    dr2.append(base_r2 - r2(y[te], pp))
                imp[f] = dict(d_spearman=float(np.mean(drho)), d_r2=float(np.mean(dr2)))
            res["permutation_importance"] = imp

        # extra held sets (the panel), scored with the model refitted on ALL held-set training rows
        if extra_eval and not shuffle_y:
            st = Standardiser().fit(X[tr])
            mr = Ridge(alpha=meta_r["alpha"]).fit(st.transform(X[tr]), ytr)
            mg = HistGradientBoostingRegressor(**{**meta_g["params"], "random_state": SEED}).fit(
                X[tr], ytr)
            for ename, edf in extra_eval:
                Xe = edf[feats].to_numpy(float)
                ye = edf[ytarget].to_numpy(float)
                res[f"extra::{ename}::ridge"] = report(ye, mr.predict(st.transform(Xe)), r)
                res[f"extra::{ename}::gbm"] = report(ye, mg.predict(Xe), r)
                res[f"extra::{ename}::mean"] = report(
                    ye, np.full(len(ye), float(ytr.mean())), r)
        res["_pred"] = dict(idx=np.flatnonzero(te).tolist(), ridge=p_ridge.tolist(),
                            gbm=p_gbm.tolist())
        return res

    # ============================ the preregistered table =====================================
    P("\n=== PRIMARY: family-held split, y = C_held (raw) ===")
    out["results"]["family_raw"] = run("family_raw", FEATS, fam_test, "y_raw",
                                       extra_eval=[("panel", panel)])
    for k in ("ridge", "gbm", "best_single", "mean"):
        v = out["results"]["family_raw"][k]
        P(f"  {k:12s} rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
          f"  R2 {v['r2']:+.4f}  AUC10 {v['auc_bottom_decile']:.4f}  n {v['n']}")

    P("\n=== KA-4: shuffled training target, same pipeline (must read at chance) ===")
    ka4 = run("shuffled", FEATS, fam_test, "y_raw", shuffle_y=True)
    P(f"  ridge rho {ka4['ridge']['spearman']:+.4f} "
      f"[{ka4['ridge']['spearman_lo']:+.4f}, {ka4['ridge']['spearman_hi']:+.4f}]  "
      f"gbm rho {ka4['gbm']['spearman']:+.4f} "
      f"[{ka4['gbm']['spearman_lo']:+.4f}, {ka4['gbm']['spearman_hi']:+.4f}]")
    ka4_pass = bool(ka4["ridge"]["spearman_lo"] < 0 < ka4["ridge"]["spearman_hi"]
                    and ka4["gbm"]["spearman_lo"] < 0 < ka4["gbm"]["spearman_hi"])
    out["ka"]["KA-4_shuffled_target"] = dict(
        ridge=ka4["ridge"]["spearman"], ridge_lo=ka4["ridge"]["spearman_lo"],
        ridge_hi=ka4["ridge"]["spearman_hi"], gbm=ka4["gbm"]["spearman"],
        gbm_lo=ka4["gbm"]["spearman_lo"], gbm_hi=ka4["gbm"]["spearman_hi"], pass_=ka4_pass)
    P(f"  KA-4 -> {'PASS' if ka4_pass else 'FAIL'}")

    P("\n=== random-gene split, y = C_held (no family one-hot) ===")
    out["results"]["random_raw"] = run("random_raw", FEATS, rnd_test, "y_raw")
    for k in ("ridge", "gbm", "best_single", "mean"):
        v = out["results"]["random_raw"][k]
        P(f"  {k:12s} rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
          f"  R2 {v['r2']:+.4f}  AUC10 {v['auc_bottom_decile']:.4f}")

    P("\n=== random-gene split, y = C_held, WITH the family one-hot ===")
    out["results"]["random_raw_fam"] = run("random_raw_fam", FEATS + FAM_COLS, rnd_test, "y_raw")
    for k in ("ridge", "gbm", "mean"):
        v = out["results"]["random_raw_fam"][k]
        P(f"  {k:12s} rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
          f"  R2 {v['r2']:+.4f}  AUC10 {v['auc_bottom_decile']:.4f}")

    P("\n=== the confounder: reliability alone, and features on top of it (family-held) ===")
    out["results"]["family_raw_rel_only"] = run("rel_only", RELIABILITY, fam_test, "y_raw")
    out["results"]["family_raw_feats_plus_rel"] = run("feats+rel", FEATS + RELIABILITY,
                                                      fam_test, "y_raw")
    for nm in ("family_raw_rel_only", "family_raw_feats_plus_rel"):
        for k in ("ridge", "gbm"):
            v = out["results"][nm][k]
            P(f"  {nm:28s} {k:5s} rho {v['spearman']:+.4f} "
              f"[{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]  R2 {v['r2']:+.4f}")

    P("\n=== declared variant: + the nuisance covariate (n pairs / n lines) ===")
    out["results"]["family_raw_nuisance"] = run("nuisance", FEATS + NUISANCE, fam_test, "y_raw")
    for k in ("ridge", "gbm"):
        v = out["results"]["family_raw_nuisance"][k]
        P(f"  {k:5s} rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
          f"  R2 {v['r2']:+.4f}")

    P("\n=== secondary target: y = C_rel (reliability-corrected), family-held ===")
    out["results"]["family_rel"] = run("family_rel", FEATS, fam_test, "y_rel",
                                       extra_eval=[("panel", panel)])
    for k in ("ridge", "gbm", "mean"):
        v = out["results"]["family_rel"][k]
        P(f"  {k:12s} rho {v['spearman']:+.4f} [{v['spearman_lo']:+.4f}, {v['spearman_hi']:+.4f}]"
          f"  R2 {v['r2']:+.4f}  AUC10 {v['auc_bottom_decile']:.4f}  n {v['n']}")

    # ---- per-feature Spearman on the held split (descriptive) --------------------------------
    P("\n=== per-feature Spearman (family-held test genes, y = C_held) ===")
    yte = held.loc[fam_test, "y_raw"].to_numpy(float)
    per_feat = {}
    rr = np.random.default_rng(SEED + 3)
    for f in FEATS + NUISANCE + RELIABILITY:
        v = held.loc[fam_test, f].to_numpy(float)
        rho, lo, hi = boot_ci(lambda a, b: spearman(a, b), yte, v, rr, n=1000)
        per_feat[f] = dict(rho=rho, lo=lo, hi=hi)
    for f, v in sorted(per_feat.items(), key=lambda kv: -abs(kv[1]["rho"])):
        P(f"  {f:30s} rho {v['rho']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]")
    out["per_feature_spearman_held"] = per_feat

    # ---- the deliverable csv -----------------------------------------------------------------
    pred = pd.DataFrame(dict(gene=held["gene"], set="held", observed=held["y_raw"],
                             observed_rel=held["y_rel"], block=held["block"],
                             gene_group1=held["gene_group1"],
                             split_family=np.where(fam_test, "held", "train"),
                             split_random=np.where(rnd_test, "held", "train")))
    for nm, key in (("pred_family_ridge", "family_raw"), ("pred_random_ridge", "random_raw")):
        col = np.full(len(held), np.nan)
        col[out["results"][key]["_pred"]["idx"]] = out["results"][key]["_pred"]["ridge"]
        pred[nm] = col
    for nm, key in (("pred_family_gbm", "family_raw"), ("pred_random_gbm", "random_raw")):
        col = np.full(len(held), np.nan)
        col[out["results"][key]["_pred"]["idx"]] = out["results"][key]["_pred"]["gbm"]
        pred[nm] = col
    pred.to_csv(WORK / "pred281.csv", index=False)
    P(f"\nwrote work/pred281.csv ({len(pred)} rows)")

    for k in list(out["results"]):
        out["results"][k].pop("_pred", None)
    out["log"] = LOG
    (WORK / "fit281.json").write_text(json.dumps(out, indent=1, default=float))
    P("wrote work/fit281.json")


if __name__ == "__main__":
    main()
