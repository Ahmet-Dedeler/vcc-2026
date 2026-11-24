"""184 Phase 2c -- what IS the basal-conditioned gain doing?

The obvious deflationary reading of a basal-conditioned gain is a detection floor: a gene the
recipient barely expresses cannot show a large delta, so "predicting" its response is predicting
a zero.  That would still be usable on the challenge (we are given the recipient's control
cells), but it is not line-specific biology and it should not be sold as such.

Test: stratify the gene axis by the MINIMUM basal expression across all 50 lines and re-measure.
On the top stratum every gene is well expressed in every line, so no detection floor can operate;
whatever survives there is genuine recipient-conditional response biology.

Also decomposes the oracle diagonal gain into its global-scalar and per-gene parts, and reports
how much of the recipient's line-specific residual each explains.

Writes work/mechanism.json only.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common import SEED, WORK, load, unit
from phase2 import cos_rows, fit_gain, gain_features, orth, sign_acc_idx, top_index


def main():
    rng = np.random.default_rng(SEED)
    R, D, genes, lines, treats, basal, mask = load("plate")
    T, L, G = R.shape
    ok = np.isfinite(R[:, :, 0])
    R = np.nan_to_num(R, nan=0.0)
    logb = np.log10(basal + 1e-8)

    minb = basal.min(0)                                  # (G,) worst-expressed line per gene
    qs = [0.0, 0.25, 0.5, 0.75, 0.9]
    thr = {q: float(np.quantile(minb, q)) for q in qs}
    out = {"thresholds": thr, "n_lines": L, "strata": {}}

    for q in qs:
        sel = np.flatnonzero(minb >= thr[q])
        acc = {a: [] for a in ["M", "gain-mean", "gain-basal", "gain-basal-clip",
                               "gain-ORACLE", "gain-basal/perm-line", "scale-ORACLE"]}
        sgn = {a: [] for a in acc}
        Rs = np.ascontiguousarray(R[:, :, sel])
        tot = Rs.sum(1)
        for held in range(L):
            S = np.asarray([l for l in range(L) if l != held])
            truth = Rs[:, held, :]
            Mfull = (tot - truth) / (L - 1)
            idx = top_index(truth, 200)
            mu, sd = logb[S][:, sel].mean(0), logb[S][:, sel].std(0) + 1e-6
            Z = (logb[:, sel] - mu) / sd
            # leave-one-out consensus inside the training set, without re-averaging
            gh = np.asarray([fit_gain(Rs[:, l, :], (tot - truth - Rs[:, l, :]) / (L - 2))
                             for l in S])
            gbar = gh.mean(0)
            X = np.concatenate([gain_features(gbar, Z[l]) for l in S])
            y = np.concatenate(list(gh))
            cf = np.linalg.lstsq(X[:, :2], y, rcond=None)[0]
            cx = np.linalg.lstsq(X, y, rcond=None)[0]
            gm = gain_features(gbar, Z[held])[:, :2] @ cf
            gb = gain_features(gbar, Z[held]) @ cx
            zn = Z[rng.choice([l for l in range(L) if l != held])]
            gor = fit_gain(truth, Mfull)
            s_or = float((truth * Mfull).sum() / (Mfull * Mfull).sum())
            preds = {"M": Mfull, "gain-mean": Mfull * gm, "gain-basal": Mfull * gb,
                     "gain-basal-clip": Mfull * np.maximum(gb, 0),
                     "gain-ORACLE": Mfull * gor, "scale-ORACLE": Mfull * s_or,
                     "gain-basal/perm-line": Mfull * (gain_features(gbar, zn) @ cx)}
            for a, P in preds.items():
                acc[a].append(float(np.nanmean(cos_rows(orth(P, Mfull), orth(truth, Mfull)))))
                sgn[a].append(float(np.nanmean(sign_acc_idx(P, truth, idx))))
        out["strata"][str(q)] = {
            "n_genes": int(len(sel)),
            "min_basal_threshold": thr[q],
            "cos_orth": {a: float(np.mean(v)) for a, v in acc.items()},
            "sign": {a: float(np.mean(v)) for a, v in sgn.items()},
        }
        s = out["strata"][str(q)]
        print(f"q={q:<4} genes {s['n_genes']:5d}  " + "  ".join(
            f"{a.split('/')[0][:9]}:{s['cos_orth'][a]:+.3f}" for a in acc))

    with open(os.path.join(WORK, "mechanism.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
