"""284 — functions COPIED from 281 (never imported: `fit281.main` checks its own preregistration
digest, and importing across experiments couples two frozen protocols).

Sources of every copy, verbatim unless marked:
  fit281.py:44-58   FEATS / NUISANCE / RELIABILITY
  fit281.py:64-111  spearman, r2, auc_bottom_decile, boot_ci, report
  fit281.py:115-175 ALPHAS, GBM_GRID, Standardiser, fit_ridge, fit_gbm, fit_best_single
  boot281.py:40-76  block_boot  (generalised: the block label is passed in, so the cluster can be
                    the CONNECTED COMPONENT rather than the first family — the only change)
  fit281.py:197-211 the 281 block construction and its check-then-append draw (as `draw_blocks`)

New here: the all-membership split (S1), the gene-family graph and its connected components (S2),
and the component 5-fold.
"""
from __future__ import annotations

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
FEAT281 = ROOT / "experiments/281-what-predicts-conservation/work/feat281.csv"
PRED281 = ROOT / "experiments/281-what-predicts-conservation/work/pred281.csv"
ZNF279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"

SEED = 20260906
HELD_FRAC = 0.20
N_BOOT = 2000

# ---- frozen feature set (281 PREREGISTRATION §3, copied from fit281.py:44-58) -----------------
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
NUISANCE = ["log_n_pairs", "n_lines"]
RELIABILITY = ["r_bar"]


# ---- metrics (fit281.py:64-111) ---------------------------------------------------------------
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
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    if len(y) < 20:
        return float("nan")
    thr = np.quantile(y, 0.10)
    lab = y <= thr
    if lab.sum() == 0 or (~lab).sum() == 0:
        return float("nan")
    r = stats.rankdata(-p)
    n1, n0 = int(lab.sum()), int((~lab).sum())
    return float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def boot_ci(fn, y, p, rng, n: int = N_BOOT):
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    if len(y) < 20:
        return float("nan"), float("nan"), float("nan")
    point = fn(y, p)
    idx = rng.integers(0, len(y), size=(n, len(y)))
    vals = np.array([fn(y[i], p[i]) for i in idx])
    return point, float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def report(y, p, rng) -> dict:
    rho, lo, hi = boot_ci(lambda a, b: spearman(a, b), y, p, rng)
    rr, rlo, rhi = boot_ci(lambda a, b: r2(a, b), y, p, rng)
    au, alo, ahi = boot_ci(lambda a, b: auc_bottom_decile(a, b), y, p, rng)
    return dict(n=int((np.isfinite(y) & np.isfinite(p)).sum()),
                spearman=rho, spearman_lo=lo, spearman_hi=hi,
                r2=rr, r2_lo=rlo, r2_hi=rhi,
                auc_bottom_decile=au, auc_lo=alo, auc_hi=ahi)


# ---- models (fit281.py:115-175) ---------------------------------------------------------------
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
    return m.predict(st.transform(Xte)), dict(alpha=float(best_a), inner_spearman=best)


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


# ---- block bootstrap (boot281.py:40-76; the block label is now an argument) --------------------
def block_boot(y, preds: dict, groups, rng, n: int = N_BOOT) -> dict:
    """Resample whole blocks with replacement; per-predictor CIs and paired deltas."""
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
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
            cur[k] = spearman(yy, pp)
            acc[k]["rho"].append(cur[k])
            acc[k]["r2"].append(r2(yy, pp))
            acc[k]["auc"].append(auc_bottom_decile(yy, pp))
        for a in keys:
            for b in keys:
                if a != b:
                    pair[f"{a}-{b}"].append(cur[a] - cur[b])
    out = {}
    for k in keys:
        out[k] = {m: dict(point=float(spearman(y, preds[k])) if m == "rho" else
                          (float(r2(y, preds[k])) if m == "r2"
                           else float(auc_bottom_decile(y, preds[k]))),
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


def paired_boot_rho(a1, a2, b1, b2, groups, rng, n: int = N_BOOT) -> dict:
    """Paired block bootstrap of Spearman(a1,a2) - Spearman(b1,b2) on the same resample."""
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    vals = []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        vals.append(spearman(a1[sel], a2[sel]) - spearman(b1[sel], b2[sel]))
    v = np.asarray(vals)
    return dict(point=float(spearman(a1, a2) - spearman(b1, b2)),
                lo=float(np.nanpercentile(v, 2.5)), hi=float(np.nanpercentile(v, 97.5)),
                sd=float(np.nanstd(v)))


def boot_rho(a, b, groups, rng, n: int = N_BOOT) -> dict:
    """Block bootstrap of a single Spearman."""
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    vals = []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        vals.append(spearman(a[sel], b[sel]))
    v = np.asarray(vals)
    return dict(point=float(spearman(a, b)), lo=float(np.nanpercentile(v, 2.5)),
                hi=float(np.nanpercentile(v, 97.5)), sd=float(np.nanstd(v)),
                n=int((np.isfinite(a) & np.isfinite(b)).sum()))


# ---- data + splits ----------------------------------------------------------------------------
def cmap(Ccross: np.ndarray, pair_idx: np.ndarray) -> np.ndarray:
    """279's conservation statistic restricted to a set of cross-line pairs (build279.py:198,
    `c_held = np.nanmean(Ccross, axis=0)`) — same dtype, same call."""
    return np.nanmean(Ccross[np.asarray(pair_idx)], axis=0)


def load_held() -> pd.DataFrame:
    df = pd.read_csv(FEAT281, low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    z = pd.read_csv(ZNF279, usecols=["gene", "c3_balanced_core"])
    held = held.merge(z, on="gene", how="left")
    return held


def memberships(held: pd.DataFrame) -> list[list[str]]:
    ga = held["gene_group_all"].fillna("").astype(str).to_numpy()
    return [[] if s == "" else s.split("|") for s in ga]


def block281(held: pd.DataFrame) -> np.ndarray:
    """fit281.py:197-199 verbatim."""
    return np.where(held["gene_group1"].fillna("").to_numpy() != "",
                    held["gene_group1"].fillna("").to_numpy(),
                    "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())


def draw_281(held: pd.DataFrame) -> tuple[np.ndarray, list[str], np.ndarray]:
    """281's family-held draw, verbatim (fit281.py:200-211)."""
    block = block281(held)
    rng = np.random.default_rng(SEED)
    blocks = np.array(sorted(set(block)))
    order = rng.permutation(len(blocks))
    counts = pd.Series(block).value_counts()
    cum, take = 0, []
    for k in blocks[order]:
        if cum >= HELD_FRAC * len(held):
            break
        take.append(k)
        cum += int(counts[k])
    te = pd.Series(block).isin(set(take)).to_numpy()
    return te, take, block


def draw_s1(held: pd.DataFrame) -> tuple[np.ndarray, list[str], np.ndarray]:
    """S1: the SAME permutation of the SAME blocks, but a gene is TEST if ANY of its HGNC groups is
    a taken family (family-less genes keep 281's cytoband blocks). The check-then-append loop stops
    when the ALL-MEMBERSHIP test set reaches 20 %."""
    block = block281(held)
    mem = memberships(held)
    rng = np.random.default_rng(SEED)
    blocks = np.array(sorted(set(block)))
    order = rng.permutation(len(blocks))
    fam_of: dict[str, list[int]] = {}
    for i, ms in enumerate(mem):
        for m in ms:
            fam_of.setdefault(m, []).append(i)
    nofam_of: dict[str, list[int]] = {}
    for i, b in enumerate(block):
        if b.startswith("NOFAM:") and not mem[i]:
            nofam_of.setdefault(b, []).append(i)
    te = np.zeros(len(held), dtype=bool)
    cum, take = 0, []
    for k in blocks[order]:
        if cum >= HELD_FRAC * len(held):
            break
        take.append(str(k))
        idx = fam_of.get(str(k), []) if not str(k).startswith("NOFAM:") else nofam_of.get(str(k), [])
        if idx:
            new = np.asarray(idx)[~te[np.asarray(idx)]]
            te[new] = True
            cum += int(len(new))
    return te, take, block


def components(held: pd.DataFrame) -> np.ndarray:
    """Connected components of: gene-gene edge when they share any HGNC group; plus, for
    family-less genes only, an edge when they share a cytoband (preserves 281's NOFAM blocking)."""
    n = len(held)
    parent = np.arange(n)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    mem = memberships(held)
    groups: dict[str, int] = {}
    for i, ms in enumerate(mem):
        for m in ms:
            if m in groups:
                union(groups[m], i)
            else:
                groups[m] = i
    block = block281(held)
    bands: dict[str, int] = {}
    for i, b in enumerate(block):
        if b.startswith("NOFAM:") and not mem[i]:
            if b in bands:
                union(bands[b], i)
            else:
                bands[b] = i
    roots = np.array([find(i) for i in range(n)])
    _, comp = np.unique(roots, return_inverse=True)
    return comp.astype(np.int64)


def draw_components(comp: np.ndarray, n_genes: int, seed: int = SEED) -> np.ndarray:
    """S2: permute components, take with 281's check-then-append loop until >= 20 % of genes."""
    rng = np.random.default_rng(seed)
    uniq = np.arange(comp.max() + 1)
    order = rng.permutation(len(uniq))
    counts = np.bincount(comp, minlength=len(uniq))
    cum, take = 0, []
    for c in uniq[order]:
        if cum >= HELD_FRAC * n_genes:
            break
        take.append(int(c))
        cum += int(counts[c])
    return np.isin(comp, take)


def draw_random(n: int) -> np.ndarray:
    """281's comparison split, verbatim (fit281.py:218-221): rng(SEED+1) = 20260907."""
    rng2 = np.random.default_rng(SEED + 1)
    perm = rng2.permutation(n)
    te = np.zeros(n, dtype=bool)
    te[perm[: int(round(HELD_FRAC * n))]] = True
    return te


def component_folds(comp: np.ndarray, k: int = 5) -> np.ndarray:
    """Assign whole components to k folds, greedily by descending component size."""
    sizes = np.bincount(comp)
    order = np.argsort(-sizes, kind="stable")
    load = np.zeros(k, dtype=np.int64)
    fold_of = np.zeros(len(sizes), dtype=np.int64)
    for c in order:
        j = int(np.argmin(load))
        fold_of[c] = j
        load[j] += int(sizes[c])
    return fold_of[comp]
