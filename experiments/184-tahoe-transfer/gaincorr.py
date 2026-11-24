"""184 Phase 2f -- the model-free version of the whole finding, on both surfaces.

For every line L and gene j, the empirical per-gene gain g_L[j] is the least-squares coefficient
of the line's own response on the other lines' consensus.  Remove the gene's mean gain across
lines and the line's mean gain across genes (a two-way additive fit); what is left is the
LINE-SPECIFIC part of the gain.  Do the same to the log basal profile.  Their correlation needs
no model, no ridge and no basis, and it is the quantity every arm in phase 2 is trying to exploit.

Reported per surface, per held line, and pooled.  Writes work/gaincorr.json only.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common import WORK, load
from crispri import CORE5, CORE6, build
from phase2 import fit_gain


def twoway(A):
    """Remove row and column means (and the grand mean) from a (lines, genes) matrix."""
    return A - A.mean(0, keepdims=True) - A.mean(1, keepdims=True) + A.mean()


def gains_and_basal(R, B):
    """R (T, L, G) centred responses, B (L, G) log basal -> line-specific gain and basal."""
    L = R.shape[1]
    tot = R.sum(1)
    Gm = np.stack([fit_gain(R[:, l, :], (tot - R[:, l, :]) / (L - 1)) for l in range(L)])
    return twoway(Gm), twoway(B.astype(np.float64)), Gm


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-300))


def report(tag, Rg, Bg, names, out):
    dG, dB, Gm = gains_and_basal(Rg, Bg)
    per = {str(names[l]): corr(dG[l], dB[l]) for l in range(len(names))}
    out[tag] = {
        "n_lines": int(len(names)), "n_genes": int(Rg.shape[2]), "n_perts": int(Rg.shape[0]),
        "pooled_corr": corr(dG.ravel(), dB.ravel()),
        "per_line_corr_mean": float(np.mean(list(per.values()))),
        "per_line_corr_min": float(np.min(list(per.values()))),
        "per_line_corr_max": float(np.max(list(per.values()))),
        "per_line": per,
        "gain_sd_across_lines": float(Gm.std(0).mean()),
        "gain_mean": float(Gm.mean()),
    }
    o = out[tag]
    print(f"{tag:10s} lines {o['n_lines']:2d} genes {o['n_genes']:5d}  "
          f"pooled r {o['pooled_corr']:+.4f}  per-line r mean {o['per_line_corr_mean']:+.4f} "
          f"[{o['per_line_corr_min']:+.4f}, {o['per_line_corr_max']:+.4f}]")


def main():
    out = {}
    R, D, genes, lines, treats, basal, mask = load("plate")
    R = np.nan_to_num(R, nan=0.0)
    report("tahoe", R, np.log10(basal + 1e-8), lines, out)
    del R, D

    for tag, names in [("crispri5", CORE5), ("crispri6", CORE6)]:
        Rc, Bc, g, p = build(names)
        report(tag, Rc, Bc, names, out)
        del Rc, Bc

    with open(os.path.join(WORK, "gaincorr.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
