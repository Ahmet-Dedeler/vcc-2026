#!/usr/bin/env python3
"""487 decision 1g: the h1-test emitter price of a same-line bank at the direction quality an IN-MIX bank has.

`whitemix487.py` puts an in-mix bank (Orion HEK293T for HEK293, K562-GWPS for K562) at a whitened direction lift of
x1.16-1.29 at w 2, below 450f's weakest calibrated bank (Pan sign-flipped p 0.3: x2.06, emitter price w 2 +0.0007,
w 0.5 +0.0166). This extends 450f's flip curve to p 0.4 / 0.45 (whitened lift ~x1.3-1.5) on the same harness, so the
w 0.5 vs w 2 choice is priced where the in-mix banks sit. Caveat: flipped rows carry noise independent of ours, an
in-mix bank's noise is partly ours; this is the closest emitter-level proxy the record allows.

It is 450's run450b + run450c (Pan's rows blended at w into the 30 covered targets, norm-preserving; `pt` = the covered
targets' block prior is the row's own sign) re-implemented here so that every output lands in THIS directory
(run450b/c create their output directories under 450's work/ at import). The imported chain is read-only:
run450 (prior 'clean4': its import-time mkdir hits the existing 450/work/h1_clean4) -> run435 -> run400.
KNOWN ANSWER: `inmix487.py 2 0.3 pt base 0 <435 args>` must reproduce 450's h1c_w2_clean4_f0.3s0 metrics exactly.

usage (after 04:30Z, one at a time, ~3 min and ~7.6 GB each):
  RUN426_EMITTER=prod uv run python experiments/487-october-decision-review/inmix487.py <w> <flip p> <pt|nopt> \
      base 0 --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse
outputs: work/h1i_w<w>_f<p>_<pt|nopt>/
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("RUN426_EMITTER", "prod")
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
W_BLEND = float(sys.argv.pop(1))
FLIP = float(sys.argv.pop(1))
PT = sys.argv.pop(1)
assert PT in ("pt", "nopt"), PT
sys.argv.insert(1, "clean4")                     # run450's argv contract: <prior> <run435 args>
sys.path.insert(0, str(ROOT / "experiments/450-same-line-programme"))
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import run450 as P450  # noqa: E402

M = P450.M
WORK = HERE / "work" / f"h1i_w{W_BLEND:g}_f{FLIP:g}_{PT}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_native = M.H.native_cap15


def pan_rows(targets, genes):
    """run450b.pan_rows verbatim for CENTRE off, NOISE off, FLIP = p, NOISE_SEED 0."""
    db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
    with h5py.File(db, "r") as f:
        ds, pg, vg = P450._cat(f["obs/dataset"]), P450._cat(f["obs/perturbed_gene_name"]), P450._cat(f["var/gene_name"])
        sel = np.sort(np.where((ds == "Pan2025_hESC") & np.isin(pg, list(map(str, targets))))[0])
        X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
        pt = pg[sel]
    vi = {g: i for i, g in enumerate(vg)}
    col = np.array([vi.get(str(g), -1) for g in genes])
    out = {}
    for t in map(str, targets):
        rows = np.where(pt == t)[0]
        if len(rows):
            with np.errstate(invalid="ignore"):
                r = np.nanmean(X[rows], axis=0)
            out[t] = np.where(col >= 0, r[np.maximum(col, 0)], np.nan)
    if FLIP:
        rng = np.random.default_rng(20260926 + 0)
        for t in sorted(out):
            v = out[t]
            m = np.flatnonzero(np.isfinite(v) & (v != 0))
            flip = rng.random(m.size) < FLIP
            v[m[flip]] = -v[m[flip]]
    return out


def native_blend(b):
    """run450b.native_blend verbatim."""
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
    print(f"[487i] blend w {W_BLEND} flip {FLIP}: {n} of {len(targets)} targets carry a same-line row", flush=True)
    return out, genes, targets


M.H.native_cap15 = native_blend

if PT == "pt":
    _ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
    GENES = _ref["genes"].astype(str)
    TARGETS = _ref["targets"].astype(str)
    _rows = pan_rows(TARGETS, GENES)          # run450c: per-target SIGNS from the (flipped) raw rows
    PTS = {t: np.where(np.isfinite(r) & (r != 0), np.sign(r), 0.0) for t, r in _rows.items()}
    print(f"[487i] per-target same-line prior for {len(PTS)} covered targets", flush=True)
    _orig = M.carrier_block

    def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
        """run450c.carrier_block verbatim."""
        t = str(TARGETS[i])
        if t not in PTS:
            return _orig(profile, libs, emit_rng, spread, key, i, emit_c)
        saved = M.H.PRIOR_FULL
        assert saved.shape == PTS[t].shape, (saved.shape, PTS[t].shape)
        M.H.PRIOR_FULL = np.where(PTS[t] != 0, PTS[t], saved)
        try:
            return _orig(profile, libs, emit_rng, spread, key, i, emit_c)
        finally:
            M.H.PRIOR_FULL = saved

    M.carrier_block = carrier_block

if __name__ == "__main__":
    M.main()
