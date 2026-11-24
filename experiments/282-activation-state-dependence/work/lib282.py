"""282 -- across-state conservation of a perturbation response within the same donors.

Reuses 280's response object verbatim (`lib280.Bank`, `profile`, `thin`, `residualise`, `cosines`,
`boot_median*`).  Nothing here changes the geometry; it only re-indexes the rows by
(donor pair, state, target) instead of (side, target).
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

ROOT = "~/Code/virtual-cell-challenge-2026"
W280 = f"{ROOT}/experiments/280-gwcd4i-donor-replicates/work"
sys.path.insert(0, W280)
from lib280 import (CONDS, DORDER, PAIRINGS, SEED, Bank, boot_median, boot_median_diff,  # noqa: E402
                    cosines, thin)

# the six donor pairs, as {pairing index: (side a, side b)}
PAIRS = [p for pr in PAIRINGS for p in pr]                       # D1D2, D3D4, D1D3, D2D4, D1D4, D2D3
COMPL = {}
for a, b in PAIRINGS:
    COMPL[a] = b
    COMPL[b] = a
STATE_PAIRS = [("Rest", "Stim8hr"), ("Rest", "Stim48hr"), ("Stim8hr", "Stim48hr")]
PANEL_CSV = f"{ROOT}/data/challenge/pert_counts.csv"


def panel_set() -> set[str]:
    return {l.split(",")[0] for l in open(PANEL_CSV).read().strip().split("\n")[1:]}


def build_index(bank: Bank, filt: str):
    """filt in {'common','keepDE','all'} -> (targets, IDX) where IDX[(cond, D, target)] is the row frame.

    'common'  : (target, donor, guide) rows that pass keep_for_DE in ALL THREE states   [282 PRIMARY]
    'keepDE'  : per-state keep_for_DE                                                   [280's filter]
    'all'     : every guide                                                             [sensitivity]
    """
    P = panel_set()
    o = bank.obs
    tg = o[(o.gtype == "targeting") & o.pert.isin(P)]
    if filt in ("common", "keepDE"):
        tg = tg[tg.keep_for_DE]
    if filt == "common":
        seen = tg.groupby(["pert", "D", "guide"]).cond.apply(set)
        ok = set(seen[seen.apply(lambda s: set(CONDS) <= s)].index)
        keep = np.array([(r, d, g) in ok for r, d, g in zip(tg.pert, tg.D, tg.guide)])
        tg = tg[keep]
    IDX = {}
    for (c, d, t), g in tg.groupby(["cond", "D", "pert"]):
        IDX[(c, d, t)] = g.sort_values("guide")
    targets = sorted({t for (_, _, t) in IDX}
                     - {t for t in {t for (_, _, t) in IDX}
                        if not all((c, d, t) in IDX for c in CONDS for d in DORDER)})
    return targets, IDX


def rows_for(IDX, cond, pair, targets):
    return {t: np.concatenate([IDX[(cond, d, t)].i.values for d in pair]) for t in targets}


def build_rows(bank, IDX, targets, budget=None, seed=SEED):
    """-> dict[(pair_str, state)] = residual matrix (n_targets, n_scored), float64.

    Each (pair, state) matrix uses that state's own control pool for those donors and is centred on its
    OWN across-target mean over `targets`; panel positions are zeroed; only the 17,471 scored columns kept.
    `budget`: dict target -> total UMI to thin every (pair, state) cell to (the depth-matched arm).
    """
    SC, CC = bank.scored, bank.cols_cd4
    rng = np.random.default_rng(seed)
    out = {}
    for cond in CONDS:
        for pair in PAIRS:
            b = bank.ctrl_profile(cond, pair)[CC]
            M = np.zeros((len(targets), len(CC)))
            rw = rows_for(IDX, cond, pair, targets)
            for k, t in enumerate(targets):
                S = bank.counts(rw[t])
                if budget is not None:
                    S = thin(S, budget[t], rng)
                M[k] = bank.profile(S)[CC] - b
            M -= M.mean(0, keepdims=True)
            M[:, ~SC] = 0.0
            out[("".join(pair), cond)] = M[:, SC]
    return out


def depth_budget(bank, IDX, targets):
    """min total UMI over the 18 (donor pair, state) cells, per target."""
    o = bank.obs
    tot = {}
    for t in targets:
        v = []
        for cond in CONDS:
            for pair in PAIRS:
                v.append(sum(float(IDX[(cond, d, t)].total.sum()) for d in pair))
        tot[t] = min(v)
    return tot


def cells_per_cell(IDX, targets):
    """median cells per (donor pair, state) side."""
    out = {}
    for cond in CONDS:
        for pair in PAIRS:
            out[("".join(pair), cond)] = np.array(
                [sum(float(IDX[(cond, d, t)].n_cells.sum()) for d in pair) for t in targets])
    return out


# ------------------------------------------------------------------ the four cells of the 2x2
def stat_S1(R, targets):
    """across state, SAME donor pair. -> (per-target mean over 6 pairs x 3 state pairs, per-cell dict)"""
    cells = {}
    for pair in PAIRS:
        for (si, sj) in STATE_PAIRS:
            cells[(f"{''.join(pair)}", f"{si}|{sj}")] = cosines(R[("".join(pair), si)], R[("".join(pair), sj)])
    agg = np.nanmean(np.vstack(list(cells.values())), axis=0)
    return agg, cells


def stat_S2(R, targets):
    """within state, DIFFERENT donor pair (280's statistic)."""
    cells = {}
    for a, b in PAIRINGS:
        for s in CONDS:
            cells[(f"{''.join(a)}|{''.join(b)}", s)] = cosines(R[("".join(a), s)], R[("".join(b), s)])
    agg = np.nanmean(np.vstack(list(cells.values())), axis=0)
    return agg, cells


def stat_S12(R, targets):
    """across state AND across donor pair."""
    cells = {}
    for a, b in PAIRINGS:
        for (si, sj) in STATE_PAIRS:
            cells[(f"{''.join(a)}|{''.join(b)}", f"{si}|{sj}")] = cosines(R[("".join(a), si)], R[("".join(b), sj)])
            cells[(f"{''.join(b)}|{''.join(a)}", f"{si}|{sj}")] = cosines(R[("".join(b), si)], R[("".join(a), sj)])
    agg = np.nanmean(np.vstack(list(cells.values())), axis=0)
    return agg, cells


def stat_W(bank, IDX, targets, seed=SEED):
    """counting floor: disjoint binomial halves of the SAME (pair, state) counts.  Half depth."""
    SC, CC = bank.scored, bank.cols_cd4
    rng = np.random.default_rng(seed + 1)
    acc = []
    for cond in CONDS:
        for pair in PAIRS:
            b = bank.ctrl_profile(cond, pair)[CC]
            M1 = np.zeros((len(targets), len(CC)))
            M2 = np.zeros((len(targets), len(CC)))
            rw = rows_for(IDX, cond, pair, targets)
            for k, t in enumerate(targets):
                Si = np.rint(bank.counts(rw[t])).astype(np.int64)
                s1 = rng.binomial(Si, 0.5).astype(np.float64)
                M1[k] = bank.profile(s1)[CC] - b
                M2[k] = bank.profile(Si.astype(np.float64) - s1)[CC] - b
            for M in (M1, M2):
                M -= M.mean(0, keepdims=True)
                M[:, ~SC] = 0.0
            acc.append(cosines(M1[:, SC], M2[:, SC]))
    return np.nanmean(np.vstack(acc), axis=0)


def sd_index(s1, s2):
    """SD = 1 - S1/S2 on targets with S2 > 0; NaN elsewhere."""
    out = np.full(len(s1), np.nan)
    ok = np.isfinite(s1) & np.isfinite(s2) & (s2 > 0)
    out[ok] = 1.0 - s1[ok] / s2[ok]
    return out


def boot_median_diff_unpaired(x, y, n=10000, seed=SEED):
    """median(x) - median(y) with INDEPENDENT bootstraps of two disjoint target groups."""
    x = np.asarray(x, float)[np.isfinite(x)]
    y = np.asarray(y, float)[np.isfinite(y)]
    rng = np.random.default_rng(seed)
    d = (np.median(x[rng.integers(0, len(x), size=(n, len(x)))], axis=1)
         - np.median(y[rng.integers(0, len(y), size=(n, len(y)))], axis=1))
    return (float(np.median(x) - np.median(y)), float(np.percentile(d, 2.5)),
            float(np.percentile(d, 97.5)), int(len(x)))


def _ordinal_rank(A):
    """row-wise ordinal ranks (ties broken by position; both inputs here are continuous)."""
    order = np.argsort(A, axis=1, kind="stable")
    r = np.empty_like(order)
    np.put_along_axis(r, order, np.broadcast_to(np.arange(A.shape[1]), A.shape), axis=1)
    return r.astype(np.float64)


def boot_spearman(x, y, n=10000, seed=SEED):
    """Spearman with a paired bootstrap over the resampled unit (targets).  Ranks are recomputed inside
    every resample (ordinal); the point estimate uses scipy's tie-correct spearmanr."""
    from scipy.stats import spearmanr
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x, float)[m], np.asarray(y, float)[m]
    rho = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    rx, ry = _ordinal_rank(x[idx]), _ordinal_rank(y[idx])
    rx -= rx.mean(1, keepdims=True)
    ry -= ry.mean(1, keepdims=True)
    d = (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
    return rho, float(np.nanpercentile(d, 2.5)), float(np.nanpercentile(d, 97.5)), int(len(x))
