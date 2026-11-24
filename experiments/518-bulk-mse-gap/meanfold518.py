"""518 — the mean-response component on the PLOLO folds against an UNCENTRED truth.

PLOLO's truth is E062.raw_truth = delta_p - mean over the held bank's ~10k perturbations: the NTC-referenced mean response
of the panel is (partly) removed, so the folds cannot see what the board's mse sees there. Here the truth is the held
line's raw delta row (size410.donor_rows: log-pseudobulk vs its own NTC, TS 5e4) on the fold's scored columns; the
prediction rows are score518's exact realised rows.
Reported per fold: the truth's mean-response share n|tbar|^2 / sum|t|^2 (noise-inflated denominator: a lower bound),
cos(our mean row, truth mean), and arms vs base: centre_full / centre_half / mean x1.5 / truth-mean ORACLE added at the
optimal scale (value of a perfect mean) -> r = cos_agg ratio - 1, literal d pds.
KNOWN ANSWER: the centred-truth base reproduces arms540's cos_agg (inside Fold()); raw_truth == raw delta - bank mean on a
spot target (checked to 1e-6 relative on the scored columns, float32 cache).

    uv run python experiments/518-bulk-mse-gap/meanfold518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
import json  # noqa: E402

import numpy as np  # noqa: E402

from score518 import Fold  # noqa: E402

WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def endpoints(rows, T, scored):
    """cos_agg (mean over frames) and literal pds (mean over frames) of realised rows vs a truth T (scored columns)."""
    cos, pds = [], []
    for c in "ABC":
        P = rows[c][:, scored]
        cos.append(float((P * T).sum() / np.sqrt((P * P).sum() * (T * T).sum())))
        Pu = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
        Tu = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
        S = Pu @ Tu.T
        own = np.diag(S).copy()
        n = S.shape[0]
        pds.append(float(np.mean(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1))))
    return float(np.mean(cos)), float(np.mean(pds))


def main():
    from size410 import donor_rows
    out = {}
    for f in FOLDS:
        F = Fold(f)
        R, have = donor_rows(f, F.targets)
        assert have.all()
        Tu = R[:, F.scored]
        Tc = F.truth                                      # PLOLO's centred truth
        # spot check: the difference Tu - Tc is one row vector (the bank mean) for every target
        diff = Tu - Tc
        spread = float(np.max(np.abs(diff - diff[0:1])) / max(np.max(np.abs(diff[0])), 1e-12))
        n = Tu.shape[0]
        tbar = Tu.mean(0)
        share_unc = float(n * (tbar ** 2).sum() / (Tu ** 2).sum())
        share_cen = float(n * (Tc.mean(0) ** 2).sum() / (Tc ** 2).sum())
        base_rows, _ = F.rows()
        hbar = np.mean([base_rows[c][:, F.scored].mean(0) for c in "ABC"], 0)
        cos_m = float((hbar * tbar).sum() / np.sqrt((hbar ** 2).sum() * (tbar ** 2).sum()))
        rec = dict(bank_mean_row_spread=spread, truth_mean_share_uncentred=share_unc, truth_mean_share_centred=share_cen,
                   cos_ourmean_truthmean=cos_m, our_mean_share=float(n * (hbar ** 2).sum() / np.mean(
                       [(base_rows[c][:, F.scored] ** 2).sum() for c in "ABC"])), arms={})
        b_cu, b_pu = endpoints(base_rows, Tu, F.scored)
        b_cc, b_pc = endpoints(base_rows, Tc, F.scored)
        for nm, kw in (("centre_full", dict(centre="effect")), ("centre_half", dict(centre=("shrink", 0.5))),
                       ("mean_x1.5", dict(centre=("shrink", -0.5)))):
            rows, _ = F.rows(**kw)
            cu, pu = endpoints(rows, Tu, F.scored)
            cc, pc = endpoints(rows, Tc, F.scored)
            rec["arms"][nm] = dict(r_uncentred=cu / b_cu - 1, dpds_uncentred=pu - b_pu, r_centred=cc / b_cc - 1, dpds_centred=pc - b_pc)
        # oracle: add alpha * tbar (the truth's own panel mean, NTC-referenced) to every delivered row, alpha optimal
        best = None
        for al in np.linspace(0, 2, 41):
            rows = {c: base_rows[c].copy() for c in "ABC"}
            for c in "ABC":
                rows[c][:, F.scored] += al * tbar[None, :]
            cu, pu = endpoints(rows, Tu, F.scored)
            if best is None or cu > best[1]:
                best = (float(al), cu, pu)
        rec["arms"]["ORACLE_add_truth_mean"] = dict(alpha=best[0], r_uncentred=best[1] / b_cu - 1, dpds_uncentred=best[2] - b_pu)
        out[f] = rec
        print(f"{f:14s} truth mean share (uncentred) {100 * share_unc:.2f} % (centred {100 * share_cen:.3f} %), our mean share "
              f"{100 * rec['our_mean_share']:.2f} %, cos(our mean, truth mean) {cos_m:+.3f}; bank-mean spread {spread:.1e}")
        for nm, a in rec["arms"].items():
            print("   " + nm.ljust(24) + " ".join(f"{k} {v:+.5f}" for k, v in a.items()))
    (WORK / "meanfold518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
