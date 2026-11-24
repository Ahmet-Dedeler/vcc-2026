"""291 — shared code. COPIED, never imported across experiments (290's modules check their own
preregistration digest; importing would couple two frozen protocols).

Sources of every copy, verbatim unless marked:
  lib290.py:39-50    FEATS (284's 21 annotation columns)
  lib290.py:66-75    FEATS290 / PE_COLS / EXTRA_COLS
  lib290.py:90-115   spearman, r2, auc_bottom_decile
  lib290.py:119-153  ALPHAS, Standardiser, fit_ridge
  lib290.py:173-183  oof_ridge
  lib290.py:200-251  block_boot, boot_rho
  lib290.py:259-279  load_held, component_folds
  build290.py:143-161 agg()  -> reimplemented here as `aggregate()` over any source subset (KA-5)

New here: the surface loaders (Broad-16 conservation per DEFINITION-276 §3 on 268's responses; the
GWCD4i published table), the single-fit-then-apply protocol (PREREG §2), and the partial Spearman.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import Ridge

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

PERSRC = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"
FEAT290 = ROOT / "experiments/290-remainder-response-features/work/feat290.csv"
FEAT281 = ROOT / "experiments/281-what-predicts-conservation/work/feat281.csv"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
ZNF279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
REM284 = ROOT / "experiments/284-conservation-family-split/work/rem284.npz"
PRED284 = ROOT / "experiments/284-conservation-family-split/work/pred284.csv"
W268 = ROOT / "experiments/268-broad16-context-residuals/work"
K562TAB = ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24/K562_comparison.suppl_table.csv"
CQCOV = ROOT / "data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv"

SEED = 20260906
N_BOOT = 2000
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
PRIMARY_COND = "Stim48hr"          # PREREG §4

# ---- feature sets (lib290.py:59-86) -----------------------------------------------------------
PROGRAMMES = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_P53_PATHWAY", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_APOPTOSIS",
              "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
PE_COLS = [f"pe_{p.replace('HALLMARK_', '')}" for p in PROGRAMMES]
FEATS290 = (["logmag", "breadth", "gen_cos", "gen_abscos", "gen_pc5", "gen_pc20"]
            + PE_COLS + ["pe_prolif_abs", "ontarget", "on_axis"])
EXTRA_COLS = ["pe_stress_abs", "n_src", "log_n_src"]

FEATS = [
    "dm_mean_effect", "dm_sd_effect", "dm_frac_dependent", "dm_frac_strong", "dm_min_effect",
    "dm_common_essential", "replogle_essential",
    "ccle_breadth", "ccle_median_tpm", "ccle_sd",
    "log_collectri_regulon_size", "collectri_indegree", "is_tf_collectri",
    "log_dorothea_regulon_size", "is_tf_dorothea",
    "humap_n_complexes", "humap_max_complex_size", "reactome_n_complexes",
    "reactome_max_complex_size",
    "log_n_paralogues", "ortholog_depth",
]

K562_SOURCES = ["gasperini_k562", "jiang_K562", "k562_essential", "k562_gwps", "perturbmatch"]


# ---- metrics (lib290.py:90-115) ---------------------------------------------------------------
def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[m], b[m]).statistic)


def r2(y, p) -> float:
    y, p = np.asarray(y, float), np.asarray(p, float)
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    sst = float(((y - y.mean()) ** 2).sum())
    return float(1.0 - ((y - p) ** 2).sum() / sst) if sst > 0 else float("nan")


# ---- models (lib290.py:119-153) ---------------------------------------------------------------
ALPHAS = np.logspace(-3, 4, 15)


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
    """lib290.fit_ridge — random inner KFold for alpha (used only by the KA-1 rebuild)."""
    from sklearn.model_selection import KFold
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
    return m.predict(st.transform(Xte)), dict(alpha=float(best_a), inner_spearman=best)


def oof_ridge(X, y, fold, min_train: int = 50):
    """lib290.oof_ridge (rem284.py:33-44)."""
    p = np.full(len(y), np.nan)
    info = {}
    ok = np.isfinite(y)
    for k in range(int(np.nanmax(fold)) + 1):
        te = ok & (fold == k)
        tr = ok & (fold != k)
        if te.sum() == 0 or tr.sum() < min_train:
            continue
        p[te], info[k] = fit_ridge(X[tr], y[tr], X[te])
    return p, info


def fit_apply_grouped(X, y, fold, Xapply):
    """PREREG §2: alpha by grouped CV over the family-component folds, then ONE fit on all rows."""
    ok = np.isfinite(y) & np.isfinite(X).any(axis=1)
    Xf, yf, gf = X[ok], y[ok], fold[ok]
    best, best_a, curve = -np.inf, ALPHAS[0], {}
    for a in ALPHAS:
        sc = []
        for k in sorted(set(gf.tolist())):
            tr, va = gf != k, gf == k
            if tr.sum() < 50 or va.sum() < 3:
                continue
            st = Standardiser().fit(Xf[tr])
            m = Ridge(alpha=a).fit(st.transform(Xf[tr]), yf[tr])
            sc.append(spearman(yf[va], m.predict(st.transform(Xf[va]))))
        s = float(np.mean(sc))
        curve[float(a)] = s
        if s > best:
            best, best_a = s, a
    st = Standardiser().fit(Xf)
    m = Ridge(alpha=best_a).fit(st.transform(Xf), yf)
    return (m.predict(st.transform(Xapply)),
            dict(alpha=float(best_a), grouped_cv_spearman=best, n_train=int(ok.sum()),
                 alpha_curve=curve))


# ---- bootstraps (lib290.py:200-251) -----------------------------------------------------------
def block_boot(y, preds: dict, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    keys = list(preds)
    acc = {k: [] for k in keys}
    pair = {f"{a}-{b}": [] for a in keys for b in keys if a != b}
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        yy = y[sel]
        cur = {k: spearman(yy, preds[k][sel]) for k in keys}
        for k in keys:
            acc[k].append(cur[k])
        for a in keys:
            for b in keys:
                if a != b:
                    pair[f"{a}-{b}"].append(cur[a] - cur[b])
    out = {k: dict(point=spearman(y, preds[k]),
                   lo=float(np.nanpercentile(acc[k], 2.5)),
                   hi=float(np.nanpercentile(acc[k], 97.5)),
                   sd=float(np.nanstd(acc[k])),
                   n=int((np.isfinite(y) & np.isfinite(preds[k])).sum()))
           for k in keys}
    out["_paired_rho"] = {k: dict(mean=float(np.nanmean(v)),
                                  lo=float(np.nanpercentile(v, 2.5)),
                                  hi=float(np.nanpercentile(v, 97.5)),
                                  p_le_0=float(np.mean(np.asarray(v) <= 0)))
                          for k, v in pair.items()}
    return out


def boot_rho(a, b, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    vals = [spearman(a[sel], b[sel]) for sel in
            (np.concatenate([idx_of[uniq[p]] for p in rng.integers(0, len(uniq), size=len(uniq))])
             for _ in range(n))]
    return dict(point=spearman(a, b), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)),
                n=int((np.isfinite(a) & np.isfinite(b)).sum()))


def _rank(x):
    x = np.asarray(x, float)
    r = np.full(len(x), np.nan)
    m = np.isfinite(x)
    r[m] = stats.rankdata(x[m])
    return r


def partial_spearman(a, b, Z) -> float:
    """Rank-transform everything, residualise a and b on Z by OLS, Pearson of the residuals."""
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    m = np.isfinite(a) & np.isfinite(b) & np.isfinite(Z).all(axis=1)
    if m.sum() < 5:
        return float("nan")
    ra, rb = _rank(np.asarray(a)[m]), _rank(np.asarray(b)[m])
    RZ = np.column_stack([_rank(Z[m, j]) for j in range(Z.shape[1])] + [np.ones(int(m.sum()))])
    ea = ra - RZ @ np.linalg.lstsq(RZ, ra, rcond=None)[0]
    eb = rb - RZ @ np.linalg.lstsq(RZ, rb, rcond=None)[0]
    if ea.std() < 1e-12 or eb.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(ea, eb)[0, 1])


def boot_partial(a, b, Z, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    vals = []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        vals.append(partial_spearman(a[sel], b[sel], Z[sel]))
    return dict(point=partial_spearman(a, b, Z), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)))


# ---- the per-source feature store, aggregated over any subset (build290.py:143-161) -----------
class Store:
    def __init__(self):
        z = np.load(PERSRC, allow_pickle=True)
        self.genes = [str(g) for g in z["genes"]]
        self.keys = [str(k) for k in z["keys"]]
        self.lines = [str(l) for l in z["lines"]]
        self.half1 = np.asarray(z["half1"])
        self.half2 = np.asarray(z["half2"])
        self.on_axis = np.asarray(z["on_axis"])
        self.Q = {k: np.asarray(z[k]) for k in z.files
                  if k not in ("genes", "keys", "lines", "half1", "half2", "on_axis")}
        self.gidx = {g: i for i, g in enumerate(self.genes)}

    def subset(self, exclude: list[str] | None = None) -> np.ndarray:
        ex = set(exclude or [])
        return np.array([i for i, k in enumerate(self.keys) if k not in ex])

    def aggregate(self, rows: np.ndarray) -> pd.DataFrame:
        """build290.agg, verbatim arithmetic, over an arbitrary source subset."""
        Q = self.Q
        out = pd.DataFrame({"gene": self.genes})
        with np.errstate(invalid="ignore"):
            mag = np.nanmean(Q["mag"][rows], axis=0)
            out["logmag"] = np.log(np.where(mag > 0, mag, np.nan))
            out["breadth"] = np.nanmean(Q["breadth"][rows], axis=0)
            gc = Q["gen_cos"][rows]
            out["gen_cos"] = np.nanmean(gc, axis=0)
            out["gen_abscos"] = np.nanmean(np.abs(gc), axis=0)
            out["gen_pc5"] = np.nanmean(Q["pc5"][rows], axis=0)
            out["gen_pc20"] = np.nanmean(Q["pc20"][rows], axis=0)
            for name in PE_COLS:
                out[name] = np.nanmean(Q[name][rows], axis=0)
            out["pe_prolif_abs"] = np.nanmean(np.abs(Q["proj_prolif"][rows]), axis=0)
            out["pe_stress_abs"] = np.nanmean(np.abs(Q["proj_stress"][rows]), axis=0)
            out["ontarget"] = np.nanmean(Q["ontarget"][rows], axis=0)
            ns = np.isfinite(Q["mag"][rows]).sum(axis=0)
            out["n_src"] = ns.astype(float)
            out["log_n_src"] = np.log1p(ns)
        out["on_axis"] = self.on_axis.astype(float)
        return out


def load_annotation() -> pd.DataFrame:
    """284's 21 annotation features on the held genes (lib290.load_held)."""
    df = pd.read_csv(FEAT281, low_memory=False)
    return df[df["set"] == "held"].reset_index(drop=True)


def check_prereg() -> str:
    h = hashlib.sha256((HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    want = (WORK / "PREREGISTRATION.sha256").read_text().split()[0]
    if h != want:
        sys.exit(f"preregistration digest mismatch ({h} != {want}) — aborting")
    return h
