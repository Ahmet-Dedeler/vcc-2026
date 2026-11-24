"""436: the per-(target, gene) donor-sign table for the MDE add on essential-screen-covered targets.

A target is COVERED when any essential screen (K562-essential, RPE1, Nadig Jurkat, Nadig HepG2) has it in its
perturbation list. On a covered target the sign of gene g is sign(mean over every donor table that covers the target
AND measures g of delta[target, g]); the tables are the four essential screens plus the three production donors
(K562-GWPS, Orion HCT116, Orion HEK293T) in `mode="all"`, the four essential screens only in `mode="ess"`; `mode="allx"`
is the all-table mean restricted to genes at least one covering essential screen measures (xcheck436.py).
0 = keep the delivered sign (uncovered target, no table measures g, mean exactly 0, the target's own gene, or fewer
than `min_agree` covering tables share the mean's sign). `h1_2025` is never read (asserted).

Donor deltas are pert-minus-control in log1p-CP50k units (`scripts/build_donor_table.py`); only their sign is used.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ESS_DIRS = {
    "k562_essential": ROOT / "data/donors/k562_essential",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "nadig_hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2",
}
PANEL_DIRS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
ESS3 = ("k562_essential", "rpe1", "nadig_jurkat")      # 434g's set (no HepG2)


def _lines(p):
    return [l.strip() for l in open(p) if l.strip()]


def sha_prefix(p, n=16):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def table_manifest(names):
    out = {}
    for nm in names:
        d = {**ESS_DIRS, **PANEL_DIRS}[nm]
        out[nm] = {"dir": str(d.relative_to(ROOT)), "delta_sha16": sha_prefix(d / "delta.npy"),
                   "perts_sha16": sha_prefix(d / "perts.txt"), "genes_sha16": sha_prefix(d / "genes.txt")}
    return out


def load_rows(names, targets, genes):
    """(len(names), len(targets), len(genes)) float32 donor deltas on the given axes; NaN where not covered/measured."""
    D = np.full((len(names), len(targets), len(genes)), np.nan, dtype=np.float32)
    for di, nm in enumerate(names):
        assert nm != "h1_2025" and "h1_2025" not in str({**ESS_DIRS, **PANEL_DIRS}[nm])
        base = {**ESS_DIRS, **PANEL_DIRS}[nm]
        P = _lines(base / "perts.txt"); G = _lines(base / "genes.txt")
        pi = {p: i for i, p in enumerate(P)}; gcol = {g: i for i, g in enumerate(G)}
        M = np.load(base / "delta.npy", mmap_mode="r")
        assert M.shape == (len(P), len(G)), (nm, M.shape)
        cols = np.array([gcol.get(str(g), -1) for g in genes]); okc = cols >= 0
        for k, t in enumerate(targets):
            if str(t) in pi:
                D[di, k, okc] = np.asarray(M[pi[str(t)]], dtype=np.float32)[cols[okc]]
    return D


def coverage(targets, ess=tuple(ESS_DIRS)):
    """dict screen -> bool per target (membership in the screen's perturbation list), plus 'union'."""
    out = {}
    for nm in ess:
        P = set(_lines(ESS_DIRS[nm] / "perts.txt"))
        out[nm] = np.array([str(t) in P for t in targets])
    out["union"] = np.any([out[nm] for nm in ess], axis=0)
    return out


def donor_signs(targets, genes, mode="all", min_agree=1, ess=tuple(ESS_DIRS)):
    """(len(targets), len(genes)) in {-1, 0, +1}; see the module doc. Returns (S, info)."""
    names = list(ess) + (list(PANEL_DIRS) if mode in ("all", "allx") else [])
    if mode not in ("all", "ess", "allx"):
        raise ValueError(mode)
    cov = coverage(targets, ess)["union"]
    D = load_rows(names, targets, genes)
    with np.errstate(all="ignore"):
        m = np.nanmean(D, 0)
        npos = (np.nan_to_num(D, nan=0.0) > 0).sum(0)
        nneg = (np.nan_to_num(D, nan=0.0) < 0).sum(0)
    s = np.sign(np.nan_to_num(m, nan=0.0))
    agree = np.where(s > 0, npos, np.where(s < 0, nneg, 0))
    s = np.where(agree >= int(min_agree), s, 0.0)
    if mode == "allx":
        # 436 xcheck436: the all-table mean, but only on genes at least one essential screen measures for this target
        # (on genes only a production donor measures, the cross-study panel-only sign is 0.51-0.53 accurate)
        s = np.where(np.isfinite(D[: len(ess)]).any(0), s, 0.0)
    s[~cov] = 0.0
    gidx = {str(g): j for j, g in enumerate(genes)}
    for k, t in enumerate(targets):
        if str(t) in gidx:
            s[k, gidx[str(t)]] = 0.0          # the target's own gene: excluded by the scorer; leave the knockdown alone
    info = {"mode": mode, "min_agree": int(min_agree), "tables": names, "n_covered_targets": int(cov.sum()),
            "covered_targets": [str(t) for t, c in zip(targets, cov) if c],
            "signed_genes_per_covered_target_median": float(np.median((s[cov] != 0).sum(1))) if cov.any() else 0.0,
            "tables_per_covered_target": {str(t): int(np.isfinite(D[:, k]).any(1).sum()) for k, t in enumerate(targets) if cov[k]}}
    return s.astype(np.float64), info
