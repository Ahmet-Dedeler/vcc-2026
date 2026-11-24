"""480 shared helpers: the scorer's nmae on a fixed reference gate, H1 truth loading, donor-row conversion.

nmae_p = sum_{g in gate_p} |lhat - l| / sum |l| over the reference's significant (BH q < 0.05), finite, non-zero,
non-own genes; perturbations with >= 10 gate genes; mean over perturbations (docs/reference/vcc2026-metrics-BRIEF.md
section 7). A gene missing from the prediction reads as 0.

Donor tables (data/donors/<d>/delta.npy) are natural-log pseudobulk effects: delta = log1p(TS*P_pert/sum) -
log1p(TS*P_ctrl/sum), TS = 5e4 (scripts/build_donor_table.py, scripts/build_orion_donor.py). The scorer's lfc is
log2(mean cell CPM pert / mean cell CPM ctrl). The exact inverse of the pseudobulk transform gives the log2 ratio of
pseudobulk fractions: log2(expm1(b_ctrl + delta) / expm1(b_ctrl)), which is the ratio-of-sums analogue of the scorer's
mean-of-ratios. `donor_lfc2` returns that, with a pseudocount in CPM units to keep near-zero donor genes finite.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
W323 = ROOT / "experiments/323-current-family-validation/work"
TS = 5e4


def gate_from(q: np.ndarray, lfc: np.ndarray, genes: np.ndarray, targets: np.ndarray, alpha: float = 0.05,
              min_gate: int = 10) -> np.ndarray:
    """Truth matrix T (targets x genes): the reference lfc on the gate, NaN elsewhere. Rows with < min_gate are all NaN."""
    T = np.full(lfc.shape, np.nan)
    for k, t in enumerate(targets):
        ok = (q[k] < alpha) & np.isfinite(lfc[k]) & (lfc[k] != 0) & (genes != t)
        if ok.sum() >= min_gate:
            T[k, ok] = lfc[k, ok]
    return T


def nmae_rows(pred: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Per-target nmae (NaN for unscored rows). pred NaN reads as 0."""
    P = np.nan_to_num(np.asarray(pred, float), nan=0.0, posinf=0.0, neginf=0.0)
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    num = np.where(G, np.abs(P - Tz), 0.0).sum(1)
    den = np.abs(Tz).sum(1)
    out = np.full(T.shape[0], np.nan)
    ok = G.any(1) & (den > 0)
    out[ok] = num[ok] / den[ok]
    return out


def nmae(pred: np.ndarray, T: np.ndarray) -> float:
    return float(np.nanmean(nmae_rows(pred, T)))


def boot_delta(rows_a: np.ndarray, rows_b: np.ndarray, n: int = 4000, seed: int = 0, sub: np.ndarray | None = None):
    """Target bootstrap (resample scored targets with replacement) of mean(rows_b) - mean(rows_a).
    Returns (delta, lo, hi, P(delta_boot >= 0))."""
    ok = np.isfinite(rows_a) & np.isfinite(rows_b)
    if sub is not None:
        ok &= sub
    d = (rows_b - rows_a)[ok]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n, len(d)))
    bs = d[idx].mean(1)
    return float(d.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)), float((bs >= 0).mean())


def load_h1():
    """h1-test truth on 323's keep axis. Returns dict(targets, genes, T, lfc, q, ctrl_cpm, cnt)."""
    r = np.load(W323 / "reference.npz", allow_pickle=True)
    keep = r["keep"].astype(bool)
    genes = r["genes"].astype(str)[keep]
    targets = r["targets"].astype(str)
    q = np.asarray(r["q"], float)
    lfc = np.asarray(r["lfc"], float)
    ref_lin = np.asarray(r["ref_lin"], float)
    ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
    lib = float(np.median(np.asarray(r["libraries"], float)))
    T = gate_from(q, lfc, genes, targets)
    return dict(targets=targets, genes=genes, T=T, lfc=lfc, q=q, ctrl_cpm=ref_lin, cnt=ref_lin * lib / 1e6,
                lib=lib)


def load_donor_rows(donor: str, targets, genes, pc_cpm: float = 0.0, clip: float | None = None,
                    root: Path | None = None):
    """Same-target donor rows on (targets x genes), converted to log2 pseudobulk ratios.
    Returns (L2, has_row, on_axis, ctrl_cpm_donor) where L2 is NaN where the donor has no row or no gene.
    pc_cpm: pseudocount in CPM added to both pert and ctrl fractions before the log ratio."""
    base = (root or ROOT / "data/donors") / donor
    P = [l.strip() for l in open(base / "perts.txt")]
    G = [l.strip() for l in open(base / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}
    gcol = {}
    for i, g in enumerate(G):
        if g and g not in gcol:
            gcol[g] = i
    M = np.load(base / "delta.npy", mmap_mode="r")
    b = np.asarray(np.load(base / "b_ctrl.npy"), float)
    cols = np.array([gcol.get(g, -1) for g in genes])
    okc = cols >= 0
    ctrl_frac = np.expm1(b) / TS                       # pseudobulk fraction of the donor control
    ctrl_cpm = np.full(len(genes), np.nan)
    ctrl_cpm[okc] = ctrl_frac[cols[okc]] * 1e6
    L2 = np.full((len(targets), len(genes)), np.nan)
    has = np.zeros(len(targets), bool)
    for k, t in enumerate(targets):
        if t not in pi:
            continue
        has[k] = True
        d = np.asarray(M[pi[t]], float)[cols[okc]]
        bc = b[cols[okc]]
        pert_cpm = np.expm1(bc + d) / TS * 1e6
        c_cpm = np.expm1(bc) / TS * 1e6
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.log2((np.maximum(pert_cpm, 0) + pc_cpm) / (c_cpm + pc_cpm))
        v[~np.isfinite(v)] = np.nan
        if clip is not None:
            v = np.clip(v, -clip, clip)
        L2[k, okc] = v
    return L2, has, okc, ctrl_cpm
