#!/usr/bin/env python3
"""450c — run450b (same-line rows blended into the effect) PLUS a per-target DE prior for the covered targets.

In the shipped carriers the block prior flips every unprotected block gene (80 %) to the prior's per-gene sign. For a
target the same-line bank covers, its own row's sign (0.76 on gate genes, 450 §2) beats any programme (0.62), so here the
covered targets' prior is the sign of Pan's row where Pan measures the gene (the default prior elsewhere). This is the
H1 analogue of 402's per-target same-line rows (`psign_by_target`), with the row from the same line.

usage: RUN426_EMITTER=prod uv run python run450c.py <w> <clean4|pan> <run435 args...>   (outputs: work/h1c_w<w>_<prior>/)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
import run450b as B  # noqa: E402  (pops <w> <prior>; installs the effect blend and the programme prior)

M = B.M
WORK = HERE / "work" / (f"h1c_w{B.W_BLEND:g}_{B.PRIOR}" + ("_c" if B.CENTRE else "") + (f"_n{B.NOISE:g}s{B.NOISE_SEED}" if B.NOISE else "")
                        + (f"_f{B.FLIP:g}s{B.NOISE_SEED}" if B.FLIP else ""))
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
GENES = _ref["genes"].astype(str)
TARGETS = _ref["targets"].astype(str)
_centre, B.CENTRE = B.CENTRE, False      # per-target SIGNS always come from Pan's raw rows (450 §2)
_rows = B.pan_rows(TARGETS, GENES)
B.CENTRE = _centre
PT = {t: np.where(np.isfinite(r) & (r != 0), np.sign(r), 0.0) for t, r in _rows.items()}
print(f"[450c] per-target same-line prior for {len(PT)} covered targets "
      f"(median {int(np.median([np.count_nonzero(v) for v in PT.values()]))} signed genes each)", flush=True)
_orig = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    t = str(TARGETS[i])
    if t not in PT:
        return _orig(profile, libs, emit_rng, spread, key, i, emit_c)
    saved = M.H.PRIOR_FULL
    assert saved.shape == PT[t].shape, (saved.shape, PT[t].shape)
    M.H.PRIOR_FULL = np.where(PT[t] != 0, PT[t], saved)
    try:
        return _orig(profile, libs, emit_rng, spread, key, i, emit_c)
    finally:
        M.H.PRIOR_FULL = saved


M.carrier_block = carrier_block

if __name__ == "__main__":
    M.main()
