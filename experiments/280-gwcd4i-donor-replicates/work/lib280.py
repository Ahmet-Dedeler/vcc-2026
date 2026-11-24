"""280 -- the response object and the donor bookkeeping for GWCD4i's pseudobulk bank.

Effect rows follow `truth267.py` / `build_donor_table.py` arithmetic exactly:
    p     = log1p(TS * S / S.sum()),  TS = 5e4
    b     = the same transform of the summed NON-TARGETING rows of the same donor(s) and condition
    delta = p - b
    r     = delta - mean(delta over the common panel-target set on the same side); panel positions zeroed
Cosines are taken on the mapped, non-panel columns of the official 18,533-gene axis (17,471 columns).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TS = 5e4
DONORS = {"CE0008162": "D1", "CE0010866": "D2", "CE0008678": "D3", "CE0006864": "D4"}
DORDER = ["D1", "D2", "D3", "D4"]
PAIRINGS = [(("D1", "D2"), ("D3", "D4")), (("D1", "D3"), ("D2", "D4")), (("D1", "D4"), ("D2", "D3"))]
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
SEED = 20280906


class Bank:
    """The extracted pseudobulk rows, indexed for (condition, donor, target, guide-rank)."""

    def __init__(self, work: str = "."):
        z = np.load(f"{work}/obs280.npz", allow_pickle=True)
        sel = np.load(f"{work}/bank280_rows.npy")
        self.X = np.load(f"{work}/bank280_X.npy", mmap_mode="r")          # (n_sel, 18129) float32 raw counts
        ax = np.load(f"{work}/axis280.npz", allow_pickle=True)
        self.cols_off = ax["cols_off"]        # official-axis positions of the mapped genes
        self.cols_cd4 = ax["cols_cd4"]        # CD4 columns, same order
        self.is_panel = ax["is_panel"]        # which of those is one of the 300 panel genes
        self.scored = ~self.is_panel          # 17,471 scored columns, in mapped order
        self.panel = [str(x) for x in ax["panel"]]
        obs = pd.DataFrame({k: z[k][sel] for k in ("donor", "cond", "gtype", "run", "pert", "guide",
                                                   "n_cells", "total", "keep_for_DE", "keep_min",
                                                   "keep_eff", "keep_tc", "keep_tg")})
        obs["donor"] = obs.donor.astype(str)
        obs["D"] = obs.donor.map(DONORS)
        for c in ("cond", "gtype", "pert", "guide"):
            obs[c] = obs[c].astype(str)
        obs["i"] = np.arange(len(obs))
        self.obs = obs
        # control pools: (cond, donor) -> summed NTC counts (float64)
        self.ctrl: dict[tuple[str, str], np.ndarray] = {}
        self.ctrl_cells: dict[tuple[str, str], float] = {}
        ntc = obs[obs.gtype != "targeting"]
        for (c, d), g in ntc.groupby(["cond", "D"]):
            self.ctrl[(c, d)] = np.asarray(self.X[np.sort(g.i.values)], dtype=np.float64).sum(0)
            self.ctrl_cells[(c, d)] = float(g.n_cells.sum())
        self.ntc_rows = {(c, d): np.sort(g.i.values) for (c, d), g in ntc.groupby(["cond", "D"])}

    # ------------------------------------------------------------------ selection
    def guides(self, cond: str, donor: str, target: str, keep_for_DE: bool = True) -> pd.DataFrame:
        o = self.obs
        m = (o.cond == cond) & (o.D == donor) & (o.pert == target) & (o.gtype == "targeting")
        if keep_for_DE:
            m &= o.keep_for_DE
        return o[m].sort_values("guide")

    def counts(self, rows) -> np.ndarray:
        rows = np.atleast_1d(np.asarray(rows, dtype=np.int64))
        if len(rows) == 0:
            return np.zeros(self.X.shape[1])
        return np.asarray(self.X[np.sort(rows)], dtype=np.float64).sum(0)

    # ------------------------------------------------------------------ the response
    @staticmethod
    def profile(S: np.ndarray) -> np.ndarray:
        tot = S.sum()
        return np.log1p(TS * S / tot) if tot > 0 else np.zeros_like(S)

    def ctrl_profile(self, cond: str, donors) -> np.ndarray:
        C = np.zeros(self.X.shape[1])
        for d in np.atleast_1d(donors):
            C = C + self.ctrl[(cond, d)]
        return self.profile(C)


def thin(S: np.ndarray, target_total: float, rng) -> np.ndarray:
    """Multinomial (binomial-per-gene) thinning of a UMI count vector to `target_total`."""
    tot = S.sum()
    if tot <= 0 or target_total >= tot:
        return S
    p = target_total / tot
    return rng.binomial(np.rint(S).astype(np.int64), p).astype(np.float64)


def residualise(deltas: np.ndarray) -> np.ndarray:
    """deltas: (n_targets, n_mapped) on the mapped columns, in official order.  Subtract the column mean
    over the rows given (the bank's own centring) -- panel positions are dropped by `scored` afterwards."""
    return deltas - deltas.mean(0, keepdims=True)


def cosines(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    out = np.full(len(A), np.nan)
    ok = (na > 0) & (nb > 0)
    out[ok] = (A[ok] * B[ok]).sum(1) / (na[ok] * nb[ok])
    return out


def boot_median_diff(x: np.ndarray, y: np.ndarray, n: int = 10000, seed: int = SEED):
    """Paired target bootstrap of median(x) - median(y); the resampled unit is the TARGET."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    d = np.median(x[idx], axis=1) - np.median(y[idx], axis=1)
    return float(np.median(x) - np.median(y)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), int(len(x))


def boot_median(x: np.ndarray, n: int = 10000, seed: int = SEED):
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    d = np.median(x[idx], axis=1)
    return float(np.median(x)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), int(len(x))
