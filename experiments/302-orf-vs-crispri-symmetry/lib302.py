"""302 — shared definitions for claim C6.  Frozen by PREREGISTRATION.md sha256 69c2a924….

Nothing here reads the reserve's `X`; it fixes the axis, the cohort and the grouping arithmetic
from metadata + the on-disk H1 CRISPRi bank only.
"""
from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/302-orf-vs-crispri-symmetry"
WORK = HERE / "work"

ATLAS = ROOT / "data/external/scperturb/JoungZhang2023_atlas.h5ad"
ATLAS_SHA = "8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22"
ORFMAP = ROOT / "experiments/301-metadata-acquisitions/work/joung_orf_map.csv"
ORFMAP_SHA = "a6b60664943abc28f3a6a7e9ef257a881adb097ca0bacc363a25f6ce29323211"
H1BANK = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
ZNF279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
PERSRC290 = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"

SEED = 20260907
N_SHARDS = 16
N_DRAWS = 8
MIN_CELLS = 50          # scored set T
MIN_CELLS_PRIMARY = 100  # primary cohort T100
CPM = 5e4
CTRL_SYMBOLS = ("GFP", "mCherry")


def read_orf_map() -> dict[int, str]:
    out = {}
    with open(ORFMAP) as fh:
        for r in csv.DictReader(fh):
            s = (r["gene_symbol"] or "").strip()
            if s:
                out[int(r["perturbation_code"])] = s
    return out


def read_h1_bank() -> tuple[list[str], list[str], np.ndarray]:
    perts = [x.strip() for x in (H1BANK / "perts.txt").read_text().splitlines() if x.strip()]
    genes = [x.strip() for x in (H1BANK / "genes.txt").read_text().splitlines() if x.strip()]
    mask0 = np.asarray(np.load(H1BANK / "measured_mask.npy", mmap_mode="r")[0])
    return perts, genes, mask0


def build_design(atlas: Path) -> dict:
    """The frozen design: axis, scored symbols, per-cell group ids and shards."""
    with h5py.File(atlas, "r") as f:
        codes = f["obs/perturbation/codes"][:].astype(np.int32)
        cats = np.asarray(f["obs/perturbation/categories"][:]).astype(str)
        ncounts = f["obs/ncounts"][:].astype(np.int64)
        var_sym = np.asarray(f["var/gene_symbol"][:]).astype(str)
        n_cells, n_var = [int(x) for x in f["X"].attrs["shape"]]
    assert len(codes) == n_cells and len(var_sym) == n_var
    cat_int = np.asarray([int(c) for c in cats])
    assert np.array_equal(cat_int, np.arange(len(cats))), "atlas categories are not 0..n-1"

    code2sym = read_orf_map()
    perts, h1_genes, mask0 = read_h1_bank()
    h1_row = {p: i for i, p in enumerate(perts)}

    # --- axis: H1 measured genes ∩ atlas var (atlas symbols are unique) --------------------
    vpos = {g: i for i, g in enumerate(var_sym)}
    assert len(vpos) == n_var, "atlas var symbols are not unique"
    axis, axis_var, axis_h1 = [], [], []
    for i, g in enumerate(h1_genes):
        if mask0[i] and g in vpos:
            axis.append(g)
            axis_var.append(vpos[g])
            axis_h1.append(i)
    axis_var = np.asarray(axis_var, np.int64)
    axis_h1 = np.asarray(axis_h1, np.int64)
    apos = {g: i for i, g in enumerate(axis)}

    # --- scored symbols ---------------------------------------------------------------------
    cells_per_code = np.bincount(codes.astype(np.int64), minlength=len(cats))
    codes_of: dict[str, list[int]] = {}
    for c, s in code2sym.items():
        if s in CTRL_SYMBOLS:
            continue
        codes_of.setdefault(s, []).append(c)
    cells_of = {s: int(sum(cells_per_code[c] for c in v)) for s, v in codes_of.items()}
    scored = sorted(s for s, n in cells_of.items() if n >= MIN_CELLS and s in h1_row)
    primary = [s for s in scored if cells_of[s] >= MIN_CELLS_PRIMARY]
    secondary = [s for s in scored if s not in set(primary)]

    ctrl_code = {}
    for c, s in code2sym.items():
        if s in CTRL_SYMBOLS:
            ctrl_code[s] = c

    # --- per-cell group ids -----------------------------------------------------------------
    # symbol groups: 0..len(scored)-1 for TFs, then GFP, mCherry; -1 = not used
    sidx = {s: i for i, s in enumerate(scored)}
    code_to_sym_group = np.full(len(cats), -1, np.int32)
    for s in scored:
        for c in codes_of[s]:
            code_to_sym_group[c] = sidx[s]
    g_gfp, g_mch = len(scored), len(scored) + 1
    code_to_sym_group[ctrl_code["GFP"]] = g_gfp
    code_to_sym_group[ctrl_code["mCherry"]] = g_mch
    n_sym_groups = len(scored) + 2

    # code groups: every ORF code of a scored symbol, then the two controls
    used_codes = sorted(c for s in scored for c in codes_of[s])
    cidx = {c: i for i, c in enumerate(used_codes)}
    code_to_code_group = np.full(len(cats), -1, np.int32)
    for c in used_codes:
        code_to_code_group[c] = cidx[c]
    c_gfp, c_mch = len(used_codes), len(used_codes) + 1
    code_to_code_group[ctrl_code["GFP"]] = c_gfp
    code_to_code_group[ctrl_code["mCherry"]] = c_mch
    n_code_groups = len(used_codes) + 2

    g_sym = code_to_sym_group[codes]
    g_code = code_to_code_group[codes]

    rng = np.random.default_rng(SEED)
    shard = rng.integers(0, N_SHARDS, size=n_cells).astype(np.int32)
    # 8 distinct balanced 8/8 partitions of the 16 shards, same generator
    draws, seen = [], set()
    while len(draws) < N_DRAWS:
        p = rng.permutation(N_SHARDS)
        a = tuple(sorted(p[:N_SHARDS // 2].tolist()))
        b = tuple(sorted(p[N_SHARDS // 2:].tolist()))
        key = min(a, b)
        if key in seen:
            continue
        seen.add(key)
        draws.append((np.asarray(a), np.asarray(b)))

    return dict(n_cells=n_cells, n_var=n_var, codes=codes, ncounts=ncounts,
                axis=axis, axis_var=axis_var, axis_h1=axis_h1, apos=apos,
                scored=scored, primary=primary, secondary=secondary,
                cells_of=cells_of, codes_of=codes_of, used_codes=used_codes,
                g_sym=g_sym, g_code=g_code, n_sym_groups=n_sym_groups,
                n_code_groups=n_code_groups, g_gfp=g_gfp, g_mch=g_mch,
                c_gfp=c_gfp, c_mch=c_mch, shard=shard, draws=draws,
                h1_row=h1_row, h1_genes=h1_genes, ctrl_code=ctrl_code)


def profile(S: np.ndarray, L: float | np.ndarray) -> np.ndarray:
    """log1p(CPM * S / L) — the map's convention, L = total counts over ALL genes."""
    L = np.asarray(L, np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log1p(CPM * S / np.where(L > 0, L, np.nan))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    from scipy.stats import rankdata
    a, b = rankdata(x[ok]), rankdata(y[ok])
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else float("nan")


def partial_spearman(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> float:
    """Spearman of x and y after linearly removing the ranked columns of Z."""
    from scipy.stats import rankdata
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    if ok.sum() < Z.shape[1] + 3:
        return float("nan")
    R = np.column_stack([rankdata(Z[ok, j]) for j in range(Z.shape[1])] + [np.ones(int(ok.sum()))])
    a, b = rankdata(x[ok]), rankdata(y[ok])
    ra = a - R @ np.linalg.lstsq(R, a, rcond=None)[0]
    rb = b - R @ np.linalg.lstsq(R, b, rcond=None)[0]
    d = np.linalg.norm(ra) * np.linalg.norm(rb)
    return float(ra @ rb / d) if d > 0 else float("nan")


def boot_ci(fn, n: int, draws: int, seed: int) -> tuple[float, float, float]:
    """Paired unit bootstrap: fn(idx) -> statistic; returns (lo, hi, P(stat <= 0))."""
    rng = np.random.default_rng(seed)
    vals = np.empty(draws)
    for d in range(draws):
        vals[d] = fn(rng.integers(0, n, n))
    v = vals[np.isfinite(vals)]
    if v.size == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), float((v <= 0).mean())
