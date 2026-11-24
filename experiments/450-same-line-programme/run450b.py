#!/usr/bin/env python3
"""450b — six-member H1 price of a SAME-LINE bank's per-target rows in the effect (bulk AND DE design follow the effect).

435's harness (production emitter, 410b key + ballast, coarse MDE), with `run400.native_cap15` monkeypatched: for each
h1-test target Pan 2025 hESC covers (30 of 85), the cross-line effect row (323's `cap15`) is replaced by
    unit(ours) + w * unit(pan)   rescaled to ours' own norm,
so the per-target magnitude (and the 410 key's input scale) is unchanged and only the direction moves. Uncovered targets
are bit-identical to `cap15`. Pan's row: mean over its guides for the target, measured genes only, mapped by symbol onto
the challenge axis (unmeasured genes and the 85 panel genes keep ours' value scaled by the blend's own factor).
Optionally the block prior is Pan's programme instead of CLEAN4 (`--prior pan`), i.e. the full same-line mode.

usage: RUN426_EMITTER=prod uv run python run450b.py <w> <clean4|pan> <run435 args...>
Outputs go to work/h1b_w<w>_<prior>/.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("RUN426_EMITTER", "prod")
W_BLEND = float(sys.argv.pop(1))
CENTRE = os.environ.get("RUN450_CENTRE") == "1"   # rows minus Pan's own programme (its mean over non-h1-test perturbations)
# 450e (2026-09-25 22:00Z): simulate a WEAKER same-line bank. Gaussian noise of RUN450_NOISE x each row's own per-gene std is
# added to every Pan row (seed RUN450_NOISE_SEED), before both the blend and the per-target signs use it.
NOISE = float(os.environ.get("RUN450_NOISE", "0") or 0)
NOISE_SEED = int(os.environ.get("RUN450_NOISE_SEED", "0") or 0)
# Gaussian noise put small non-zeros everywhere and broke the emitter's per-cell totals at w 1 (a harness artefact). The
# realistic weakening is RUN450_FLIP = p: flip the sign of a random fraction p of each row's measured genes (magnitudes and
# sparsity kept), which lowers the bank's sign accuracy directly.
FLIP = float(os.environ.get("RUN450_FLIP", "0") or 0)
PRIOR = sys.argv.pop(1)
assert PRIOR in ("clean4", "pan"), PRIOR
sys.argv.insert(1, PRIOR)          # run450's own argv contract: <prior> <run435 args>
import run450 as P450  # noqa: E402  (sets the prior hook and imports run435 as P450.M)

M = P450.M
WORK = HERE / "work" / (f"h1b_w{W_BLEND:g}_{PRIOR}" + ("_c" if CENTRE else "") + (f"_n{NOISE:g}s{NOISE_SEED}" if NOISE else "")
                        + (f"_f{FLIP:g}s{NOISE_SEED}" if FLIP else ""))
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_native = M.H.native_cap15
T_H1 = sorted(P450.H1_TARGETS)


def pan_rows(targets, genes):
    """Pan's mean row per covered target on the challenge axis `genes` (NaN where Pan does not measure a gene)."""
    db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
    with h5py.File(db, "r") as f:
        ds, pg, vg = P450._cat(f["obs/dataset"]), P450._cat(f["obs/perturbed_gene_name"]), P450._cat(f["var/gene_name"])
        sel = np.sort(np.where((ds == "Pan2025_hESC") & np.isin(pg, list(map(str, targets))))[0])
        X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
        pt = pg[sel]
    vi = {g: i for i, g in enumerate(vg)}
    col = np.array([vi.get(str(g), -1) for g in genes])
    centre = 0.0
    if CENTRE:
        with h5py.File(db, "r") as f:
            rest = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, T_H1))[0])
            Xr = np.where(f["layers/measure_mask"][rest, :] > 0, f["X"][rest, :], np.nan)
        with np.errstate(invalid="ignore"):
            m = np.nanmean(Xr, axis=0)
        centre = np.where(col >= 0, np.nan_to_num(m)[np.maximum(col, 0)], 0.0)
    out = {}
    for t in map(str, targets):
        rows = np.where(pt == t)[0]
        if len(rows):
            with np.errstate(invalid="ignore"):
                r = np.nanmean(X[rows], axis=0)
            out[t] = np.where(col >= 0, r[np.maximum(col, 0)], np.nan) - centre
    if FLIP:
        rng = np.random.default_rng(20260926 + NOISE_SEED)
        for t in sorted(out):
            v = out[t]
            m = np.flatnonzero(np.isfinite(v) & (v != 0))
            flip = rng.random(m.size) < FLIP
            v[m[flip]] = -v[m[flip]]
    if NOISE:
        rng = np.random.default_rng(20260925 + NOISE_SEED)
        for t in sorted(out):
            v = out[t]
            m = np.isfinite(v)
            v[m] = v[m] + rng.normal(0.0, NOISE * float(np.std(v[m])), int(m.sum()))
    return out


def native_blend(b):
    eff, genes, targets = _native(b)
    if W_BLEND == 0:
        return eff, genes, targets
    rows = pan_rows(targets, genes)
    panel = np.isin(genes.astype(str), list(map(str, targets)))
    out = eff.copy()
    n = 0
    for i, t in enumerate(map(str, targets)):
        if t not in rows:
            continue
        o = eff[i]
        p = np.where(np.isfinite(rows[t]) & ~panel, rows[t], 0.0)
        no, np_ = np.linalg.norm(o[~panel]), np.linalg.norm(p)
        if no == 0 or np_ == 0:
            continue
        mix = o / no + W_BLEND * p / np_
        out[i] = np.where(panel, o, mix * (no / np.linalg.norm(mix[~panel])))
        n += 1
    print(f"[450b] blend w {W_BLEND}: {n} of {len(targets)} targets carry a same-line row (norms preserved)", flush=True)
    return out, genes, targets


M.H.native_cap15 = native_blend

if __name__ == "__main__":
    M.main()
