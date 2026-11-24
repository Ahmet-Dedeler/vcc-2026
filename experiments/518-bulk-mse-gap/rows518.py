"""518 — the delivered-row route of arms540 (lib420's own functions, no operator): sharpen q per frame (row mode on the frame's
gate), render (minnorm), realise by 414's ballast emulation at KMAX 16 against PP's unkeyed carriers. dump518's known answer
checks it against arms540's saved rows; score518 uses it on the dumped compositions.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
import numpy as np  # noqa: E402

from lib420 import TS, frames, profile_from_effect, realise_k, sharpen  # noqa: E402


def challenge_frames():
    return dict(frames())


def de_profiles(carrier_design, fr):
    return {ctx: [profile_from_effect(carrier_design[i], b, gate, "minnorm") for i in range(carrier_design.shape[0])]
            for ctx, (b, gate) in fr.items()}


def realised_rows(bulk, fr, P_DE, q=1.75, scale=None, gene_mult=None, return_K=False):
    """{ctx: realised delivered rows}. `scale` = {ctx: factor} after sharpening (VCC_BULK_CTX_SCALE); `gene_mult` = optional
    {ctx: per-gene multiplier on the sharpened effect} (a gene-class allocation arm)."""
    out, Ks = {}, {}
    for ctx, (b, gate) in fr.items():
        F = sharpen(bulk, q, gate, mode="row") if q is not None else bulk.copy()
        if scale is not None:
            F = F * float(scale[ctx])
        if gene_mult is not None:
            F = F * gene_mult[ctx][None, :]
        rea = np.empty_like(F)
        K = np.empty(F.shape[0], int)
        for i in range(F.shape[0]):
            p_bk = profile_from_effect(F[i], b, gate, "minnorm")
            p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
            rea[i] = np.log1p(TS * p_re) - b
        out[ctx] = rea
        Ks[ctx] = K
    return (out, Ks) if return_K else out
