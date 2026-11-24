"""169 — the one definition of magnitude sharpening, shared by every script here.

    P -> sign(P) * |P|^q          on the gated coordinates, then renormalised.

TWO RENORMALISATIONS, AND WHY THE CHOICE IS NOT COSMETIC
--------------------------------------------------------
`|P|^q` changes the overall scale by orders of magnitude, so something has to put it back.

`median` -- one global factor restoring the panel's MEDIAN row norm.  This is what 030 and
            168 measured.  It leaves the ACROSS-target norm spread free to move: a row whose
            magnitudes are concentrated gains norm relative to a flat one.
`row`    -- one factor per row restoring that row's OWN norm.  The transform is then purely
            the within-row magnitude shape and every row's amplitude is untouched.

`pds_cosine` is invariant to per-row scaling, so the two are IDENTICAL on `pds` -- 168's `pds`
table is a statement about both.  They differ on `nmae` (which reads the row's scale against
the reference's), on `mse` (quadratic in amplitude, and a15 sits on the rho cliff at 99.7% of
the dispersion cap), and on the emitted p-values.  Which one to ship is therefore an
empirical question and is measured in `sharpen169.py`, not assumed.
"""
from __future__ import annotations

import numpy as np


def sharpen(D: np.ndarray, q: float, gate: np.ndarray, mode: str = "median") -> np.ndarray:
    """`sign(D)*|D|^q` on `gate`, renormalised. `D` may be 1-D (one row) or 2-D (a panel)."""
    if mode not in {"median", "row"}:
        raise ValueError(f"unknown renormalisation {mode!r}")
    one = D.ndim == 1
    A = np.atleast_2d(np.asarray(D, dtype=np.float64))
    G = np.where(np.atleast_1d(gate)[None, :], A, 0.0)
    if q == 1.0:
        return G[0].copy() if one else G.copy()
    out = np.sign(G) * np.abs(G) ** q
    n0 = np.linalg.norm(G, axis=1)
    n1 = np.linalg.norm(out, axis=1)
    if mode == "median":
        f = float(np.median(n0)) / max(float(np.median(n1)), 1e-300)
        out = out * f
    else:
        out = out * (n0 / np.maximum(n1, 1e-300))[:, None]
    return out[0] if one else out


def sharpen_lfc(D: np.ndarray, q: float, b: np.ndarray, gate: np.ndarray,
                mode: str = "row", restore: str = "h") -> np.ndarray:
    """Same power transform, but applied in LOG2 FOLD-CHANGE space instead of log-pseudobulk.

    WHY THIS IS A DIFFERENT ARM, NOT A REPARAMETERISATION
    ----------------------------------------------------
    The builder's effect lives in `log1p(TS * p)` space, and `log1p` is nearly LINEAR below
    `TS*p ~ 1` (i.e. under ~20 CPM) and logarithmic above it.  So for one and the same fold
    change, `|h|` is ~2.4x larger on a 100-CPM gene than on a 10-CPM one.  A power transform
    of `h` therefore concentrates on HIGH-EXPRESSION genes; a power transform of `lfc`
    concentrates on LARGE-FOLD-CHANGE genes.

    That distinction is exactly the six-member split: `pds` and `mse` read the effect in the
    `h` space (the scorer forms `bulk_lognorm_means(pred) - bulk_lognorm_means(ctrl)`), while
    `nmae`, `reach`, `fid` and `jac` all read log2 fold changes.  169 measured the
    consequence: in `h` space the `nmae` cost on h1-test at the board-validated gate depth is
    `+0.0416` raw at q = 1.15, and in `lfc` space the same q reads `-0.0140` -- a GAIN.

        u = expm1(b + h) / expm1(b)        linear fold change, pre-normalisation
        L = log2(u)  ->  sharpen  ->  L'
        h' = log1p(expm1(b) * 2**L') - b   on the gate, exactly 0 off it

    `restore` decides WHICH row norm is held fixed, and it is not cosmetic.  Preserving the
    `lfc` norm leaves the h-space amplitude free, and on the challenge panel that reads
    `rA = 0.833` at q = 1.15 -- a 9% amplitude cut riding along with the shape change, i.e.
    TWO variables in one flag (AGENTS.md 11.1 rule 5, and it is the v3 mistake).  `restore='h'`
    rescales the mapped-back effect so the row's L2 norm in the BUILDER's space -- the space
    `VCC_ALPHA` acts in -- is exactly what it was, which makes the arm amplitude-neutral by
    construction and leaves 096's whole `s`-parameterisation of `mse` untouched.

    Exact and invertible; `q = 1` returns the gated input bit-for-bit.
    """
    one = D.ndim == 1
    A = np.atleast_2d(np.asarray(D, dtype=np.float64))
    on = np.asarray(gate, dtype=bool)
    bb = np.asarray(b, dtype=np.float64)
    base = np.expm1(bb)
    live = on & (base > 0)
    G = np.where(live[None, :], A, 0.0)
    if q == 1.0:
        return G[0].copy() if one else G.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(live[None, :], np.expm1(bb[None, :] + G) / np.where(base > 0, base, 1.0),
                     1.0)
    u = np.clip(u, 1e-12, None)
    L = np.where(live[None, :], np.log2(u), 0.0)
    Lq = sharpen(L, q, live, mode=mode)
    out = np.where(live[None, :], np.log1p(base[None, :] * np.exp2(Lq)) - bb[None, :], 0.0)
    if restore == "h":
        n0 = np.linalg.norm(G, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        out = out * (n0 / np.maximum(n1, 1e-300))[:, None]
    elif restore != "lfc":
        raise ValueError(f"unknown restore {restore!r}")
    return out[0] if one else out
