"""290 — shared code, COPIED (never imported across experiments: `rem284.main` / `build279.main`
check their own preregistration digests, and importing couples two frozen protocols).

Sources of every copy, verbatim unless marked:
  lib284.py:39-50    FEATS / NUISANCE / RELIABILITY
  lib284.py:54-100   spearman, r2, auc_bottom_decile, boot_ci, report
  lib284.py:104-163  ALPHAS, GBM_GRID, Standardiser, fit_ridge, fit_gbm
  lib284.py:167-236  block_boot, paired_boot_rho, boot_rho
  lib284.py:240-263  cmap, load_held, memberships, block281
  lib284.py:376-386  component_folds
  rem284.py:33-44    oof_ridge (the grouped out-of-fold ridge)
  build279.py:80-146 residual_rows, load_unit, pair_cos   (the row reader; KA-1 asserts it)
  lib276.py:68-141   read_bank, common_axis_paths, load_source_path, cosine_rows (imported by path
                     from 276, as 279 does — same object, no copy, so KA-1 tests the real reader)

New here: the half-source split, the response-derived feature block (PREREG §5) and the hallmark
indicator directions.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

FEAT281 = ROOT / "experiments/281-what-predicts-conservation/work/feat281.csv"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
SPLITS284 = ROOT / "experiments/284-conservation-family-split/work/splits284.npz"
REM284 = ROOT / "experiments/284-conservation-family-split/work/rem284.npz"
REM284J = ROOT / "experiments/284-conservation-family-split/work/rem284.json"
GMT = ROOT / "data/external/broad16line/metadata/hallmark_symbols.gmt"

sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(ROOT / "experiments/279-znf-state-dependence"))
import lib276 as L276  # noqa: E402
import inv279 as I279  # noqa: E402

SEED = 20260906
N_BOOT = 2000
BLK = 8192

# ---- PREREG §2: the primary partition's two halves, declared before any 290 number ------------
HALF1_LINES = ["DLD-1", "HCT116", "HEK293T", "HT-29", "K562", "MCF7", "RPE1", "THP-1",
               "iPSC-Tian19", "melanoma"]
HALF2_LINES = ["A549", "BxPC-3", "Calu-3", "H1", "HAP1", "HEK293", "HepG2", "Jurkat",
               "TeloHAEC", "cardiomyocyte", "iPSC-neuron"]

# ---- PREREG §5(d): the frozen programme list --------------------------------------------------
PROGRAMMES = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_P53_PATHWAY", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_APOPTOSIS",
              "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
PROLIF = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1"]
STRESS = ["HALLMARK_P53_PATHWAY", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_APOPTOSIS"]

PE_COLS = [f"pe_{p.replace('HALLMARK_', '')}" for p in PROGRAMMES]
# PREREG §5: the frozen 17-column model feature set
FEATS290 = (["logmag", "breadth", "gen_cos", "gen_abscos", "gen_pc5", "gen_pc20"]
            + PE_COLS + ["pe_prolif_abs", "ontarget", "on_axis"])
EXTRA_COLS = ["pe_stress_abs", "n_src", "log_n_src"]   # reported, NOT in the model
ALL_COLS = FEATS290 + EXTRA_COLS

PC_MIN_ROWS = 200      # PREREG §5(c3)
BREADTH_Q = 99.0       # PREREG §5(b)

# ---- 281's frozen annotation feature set (lib284.py:39-50) ------------------------------------
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


# ---- metrics (lib284.py:54-100) ---------------------------------------------------------------
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


# ---- models (lib284.py:104-163) ---------------------------------------------------------------
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


# ---- out-of-fold predictors (rem284.py:33-44) -------------------------------------------------
def oof_ridge(X, y, fold, min_train: int = 50):
    p = np.full(len(y), np.nan)
    info = {}
    ok = np.isfinite(y)
    for k in range(int(fold.max()) + 1):
        te = ok & (fold == k)
        tr = ok & (fold != k)
        if te.sum() == 0 or tr.sum() < min_train:
            continue
        p[te], info[k] = fit_ridge(X[tr], y[tr], X[te])
    return p, info


def oof_gbm(X, y, fold, min_train: int = 50):
    p = np.full(len(y), np.nan)
    info = {}
    ok = np.isfinite(y)
    for k in range(int(fold.max()) + 1):
        te = ok & (fold == k)
        tr = ok & (fold != k)
        if te.sum() == 0 or tr.sum() < min_train:
            continue
        p[te], info[k], _ = fit_gbm(X[tr], y[tr], X[te])
    return p, info


# ---- bootstraps (lib284.py:167-236) -----------------------------------------------------------
def block_boot(y, preds: dict, groups, rng, n: int = N_BOOT) -> dict:
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


def boot_rho(a, b, groups, rng, n: int = N_BOOT) -> dict:
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


# ---- data + splits (lib284.py:240-263, 376-386) ------------------------------------------------
def cmap(Ccross: np.ndarray, pair_idx: np.ndarray) -> np.ndarray:
    return np.nanmean(Ccross[np.asarray(pair_idx)], axis=0)


def load_held() -> pd.DataFrame:
    df = pd.read_csv(FEAT281, low_memory=False)
    return df[df["set"] == "held"].reset_index(drop=True)


def block281(held: pd.DataFrame) -> np.ndarray:
    return np.where(held["gene_group1"].fillna("").to_numpy() != "",
                    held["gene_group1"].fillna("").to_numpy(),
                    "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())


def component_folds(comp: np.ndarray, k: int = 5) -> np.ndarray:
    sizes = np.bincount(comp)
    order = np.argsort(-sizes, kind="stable")
    load = np.zeros(k, dtype=np.int64)
    fold_of = np.zeros(len(sizes), dtype=np.int64)
    for c in order:
        j = int(np.argmin(load))
        fold_of[c] = j
        load[j] += int(sizes[c])
    return fold_of[comp]


# ---- the row reader (build279.py:80-146, verbatim) --------------------------------------------
def residual_rows(rel: str, axis: list[str], genes: list[str], row_blk: int = 1024):
    perts, bank_genes, p = L276.read_bank(rel)
    gidx = {g: i for i, g in enumerate(bank_genes)}
    cols = np.asarray([gidx[g] for g in axis], dtype=np.int64)
    pidx = {t: i for i, t in enumerate(perts)}
    have = [(t, pidx[t]) for t in genes if t in pidx]
    if not have:
        return np.zeros(0, np.int64), np.zeros((0, len(axis)), np.float64), np.zeros(len(axis))
    tidx = {t: i for i, t in enumerate(genes)}
    order = np.argsort([r for _, r in have])
    rows_sorted = np.asarray([have[i][1] for i in order], dtype=np.int64)
    idx = np.asarray([tidx[have[i][0]] for i in order], dtype=np.int64)
    dest = np.empty(idx.size, np.int64)
    dest[np.argsort(idx)] = np.arange(idx.size)
    D = np.load(p / "delta.npy", mmap_mode="r")
    blk = np.empty((rows_sorted.size, len(axis)), dtype=np.float64)
    for lo in range(0, rows_sorted.size, row_blk):
        r = rows_sorted[lo:lo + row_blk]
        blk[dest[lo:lo + r.size]] = np.asarray(D[r.min():r.max() + 1][:, cols],
                                               dtype=np.float64)[r - r.min()]
    D._mmap.close(); del D
    mu = blk.mean(axis=0)                            # NEW: the generic direction, kept (PREREG §5c1)
    blk -= mu[None, :]                               # centre on this source's across-gene mean
    return np.sort(idx), blk, mu


def load_unit(rel: str, axis: list[str], genes: list[str]):
    """residual rows (DEFINITION-276 §2) restricted to `genes`, unit-normalised, plus their index.

    Returns additionally the CENTRED unnormalised block (for magnitude / breadth) and the bank's
    mean raw row `mu`; `build279.load_unit` returns only (idx, rows, nrm) and the arithmetic of
    those three is unchanged (KA-1)."""
    idx, blk, mu = residual_rows(rel, axis, genes)
    if idx.size == 0:
        return idx, np.zeros((0, len(axis)), np.float32), np.zeros(0), blk, mu
    rows = blk.astype(np.float32)
    nrm = np.linalg.norm(rows, axis=1)
    ok = nrm > 0
    idx, rows, nrm, blk = idx[ok], rows[ok], nrm[ok], blk[ok]
    rows /= nrm[:, None]
    return idx, rows, nrm, blk, mu


def pair_cos(idx_i, U_i, idx_j, U_j):
    common = np.intersect1d(idx_i, idx_j, assume_unique=True)
    if common.size == 0:
        return common, np.zeros(0, np.float32)
    pi = np.searchsorted(idx_i, common)
    pj = np.searchsorted(idx_j, common)
    out = np.empty(common.size, np.float32)
    for lo in range(0, common.size, BLK):
        hi = min(lo + BLK, common.size)
        out[lo:hi] = (U_i[pi[lo:hi]] * U_j[pj[lo:hi]]).sum(axis=1)
    return common, out


# ---- 279's source scan, reproduced (build279.main:157-169) -------------------------------------
def scan_sources():
    per_src = {}
    for k, (rel, line, ds) in I279.CANDIDATES.items():
        if not (ROOT / rel / "perts.txt").exists():
            continue
        perts, genes, _ = L276.read_bank(rel)
        per_src[k] = dict(line=line, dataset=ds, rel=rel, perts=set(perts), genes=set(genes))
    panel, official = set(L276.PANEL), set(L276.OFFICIAL)
    lines: dict[str, set[str]] = {}
    for k, d in per_src.items():
        for g in d["perts"]:
            if g in panel or g not in official:
                continue
            lines.setdefault(g, set()).add(d["line"])
    genes_held = sorted(g for g, ls in lines.items() if len(ls) >= 2)
    keys = sorted(k for k, d in per_src.items() if len(d["perts"] & set(genes_held)) >= 5)
    axis = L276.common_axis_paths([per_src[k]["rel"] for k in keys])
    return per_src, genes_held, keys, axis


# ---- hallmark indicators (PREREG §5d) ----------------------------------------------------------
def hallmark_sets() -> dict[str, set[str]]:
    out = {}
    for line in GMT.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            out[f[0]] = {x for x in f[2:] if x}
    return out


def programme_directions(axis: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Unit-normalised indicator directions on the axis: (W (9, |axis|), w_prolif, w_stress, sizes)."""
    hs = hallmark_sets()
    aidx = {g: i for i, g in enumerate(axis)}
    W = np.zeros((len(PROGRAMMES), len(axis)))
    sizes = {}
    for i, p in enumerate(PROGRAMMES):
        members = hs[p]
        on = [aidx[g] for g in members if g in aidx]
        W[i, on] = 1.0
        sizes[p] = dict(n_set=len(members), n_on_axis=len(on))
        W[i] /= max(np.linalg.norm(W[i]), 1e-12)

    def combo(names):
        v = np.zeros(len(axis))
        for p in names:
            on = [aidx[g] for g in hs[p] if g in aidx]
            v[on] = 1.0
        return v / max(np.linalg.norm(v), 1e-12)

    return W, combo(PROLIF), combo(STRESS), sizes


def check_prereg() -> str:
    h = hashlib.sha256((HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    want = (WORK / "PREREGISTRATION.sha256").read_text().split()[0]
    if h != want:
        sys.exit(f"preregistration digest mismatch ({h} != {want}) — aborting")
    return h
