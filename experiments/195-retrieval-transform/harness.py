"""195 — the 019 LODO harness, rebuilt on the compact cache.

Every object here is a re-expression of `experiments/019-whiten-interact/interact.py` on the
`feat` axis (7,760 genes = the 3-line gene intersection minus the 259 shared VCC panel
targets).  The rng draw order for the whitener basis subsample is replicated exactly --
`default_rng(0)`, folds in DONORS order, predictors in DONORS order within a fold -- because
the operator depends on which 6,000 perturbations land in the basis.

Phase 0 checks the rebuild against the three numbers AGENTS.md 4.2 publishes.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = HERE / "work" / "cache"
LODO = ["k562_gwps", "orion_hct116", "orion_hek293t"]
BASIS_PERTS = 6000


def sha256_arr(*arrs: np.ndarray) -> str:
    h = hashlib.sha256()
    for a in arrs:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def pds_cosine(P: np.ndarray, T: np.ndarray) -> float:
    """019's scorer verbatim: query = row of P, gallery = all rows of T, cosine distance."""
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


def pds_rows(P: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Per-query pds, so folds can be resampled without recomputing the operator."""
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return 1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)


class Bank:
    """One donor on the feat axis: raw deviations, gene_norm'd deviations, pert index."""

    def __init__(self, name: str):
        self.name = name
        self.raw = np.load(CACHE / f"{name}.s_raw.npy", mmap_mode="r")
        self.scale = np.load(CACHE / f"{name}.scale.npy")
        self.perts = (ROOT / "data" / "donors" / name / "perts.txt").read_text().split("\n")
        self.idx = {p: i for i, p in enumerate(self.perts)}
        self._gn_median_norm: float | None = None

    def gn(self, rows: np.ndarray | None = None) -> np.ndarray:
        """gene_norm'd deviations (float32, as in 019) for the given perturbation rows."""
        src = self.raw if rows is None else self.raw[rows]
        return (np.asarray(src, np.float32) / self.scale).astype(np.float32)

    def gn_median_norm(self) -> float:
        if self._gn_median_norm is None:
            acc = np.empty(self.raw.shape[0])
            for i in range(0, self.raw.shape[0], 2048):
                acc[i:i + 2048] = np.linalg.norm(self.gn(np.arange(i, min(i + 2048, self.raw.shape[0]))), axis=1)
            self._gn_median_norm = float(np.median(acc)) or 1.0
        return self._gn_median_norm


def load_banks() -> dict[str, Bank]:
    return {n: Bank(n) for n in LODO}


def common_perts() -> list[str]:
    return (CACHE / "common.txt").read_text().split("\n")


def fold_basis_rows(banks: dict[str, Bank]) -> dict[str, dict[str, np.ndarray]]:
    """Replicate 019's rng consumption exactly: one default_rng(0) for the whole treat."""
    rng = np.random.default_rng(0)
    out: dict[str, dict[str, np.ndarray]] = {}
    for held in LODO:
        others = [n for n in LODO if n != held]
        sel = {}
        for n in others:
            k = banks[n].raw.shape[0]
            sel[n] = np.sort(rng.choice(k, BASIS_PERTS, False)) if k > BASIS_PERTS else np.arange(k)
        out[held] = sel
    return out


def build_fold(held: str, banks: dict[str, Bank], common: list[str],
               basis_rows: dict[str, np.ndarray]) -> dict:
    """`parts` (per-donor scale-matched panel predictions), `T` (held-line truth), basis."""
    others = [n for n in LODO if n != held]
    sc = {n: banks[n].gn_median_norm() for n in others}
    ref = float(np.median(list(sc.values())))
    crow = {n: np.fromiter((banks[n].idx[p] for p in common), np.int64, len(common))
            for n in LODO}
    T = np.asarray(banks[held].raw[crow[held]], np.float64)
    parts, Bl = {}, []
    for n in others:
        parts[n] = (ref / sc[n]) * banks[n].gn(crow[n]).astype(np.float64)
        Bl.append((ref / sc[n]) * banks[n].gn(basis_rows[n]).astype(np.float64))
    return dict(held=held, others=others, T=T, parts=parts, ref=ref, sc=sc,
                basis=np.vstack(Bl))


def fit_operator(basis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`vcc2026.whiten.Whitener.fit` verbatim (eigh of the Gram, relative eigenvalue floor)."""
    C = (basis.T @ basis) / basis.shape[0]
    ev, V = np.linalg.eigh(C)
    ev = np.maximum(ev, ev.max() * 1e-12)
    return ev, V


def apply_operator(P: np.ndarray, ev: np.ndarray, V: np.ndarray,
                   gamma: float, lam: float = 0.0) -> np.ndarray:
    if gamma == 0.0:
        return P
    w = (ev + lam) ** (-gamma / 2.0)
    return ((P @ V) * w) @ V.T
