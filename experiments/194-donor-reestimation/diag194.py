"""194 — does each variant do the thing it is supposed to do, before asking whether it helps?

The decisive test for knockdown weighting is MECHANISTIC and does not need a transfer
score.  A cell whose guide failed carries a control profile, so it dilutes the pooled
effect toward zero.  Down-weighting those cells must therefore make the estimated effect
BIGGER along its own direction: writing `d_w` for the reweighted effect and `d_0` for the
pooled one,

    projection  =  <d_w - d_0, d_0> / <d_0, d_0>

should be POSITIVE and the amplification `|d_w| / |d_0|` should exceed 1.  The permuted
control (`kd_perm`, `kd_perm_depth`) reweights by exactly the same weight distribution
with the cell identities scrambled, so it changes the estimate by a comparable amount but
must NOT amplify.  If the real weighting's projection is not clearly above the permuted
one, the per-cell readout is not identifying escapers and no transfer number can rescue it.

Also reported, per the mandate: retained/reweighted cell counts and the effect-norm
distribution of every table.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from tables import SumBank

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"

PANEL = {l.split(",")[0] for l in
         (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()}


def analyse(donor: str, mode: str, variants: list[str]) -> dict:
    bank = SumBank(donor, mode)
    d0, b0, names = bank.table("v0")
    d0 = d0.astype(np.float64)
    scored = np.array([i for i, g in enumerate(bank.genes) if g and g not in PANEL])
    keep_rows = np.array([i for i, t in enumerate(names) if t not in PANEL])
    A = d0[np.ix_(keep_rows, scored)]
    n0 = np.linalg.norm(A, axis=1)
    wsum = np.load(bank.dir / "_wsum.npy")
    cells = np.load(bank.dir / "_cells.npy")
    kd = np.load(bank.dir / "_kd_target.npy")
    out = dict(donor=donor, mode=mode, n_targets=len(names), n_scored_genes=len(scored),
               cells_per_target=dict(min=int(cells[:-1].min()),
                                     median=float(np.median(cells[:-1])),
                                     max=int(cells[:-1].max())),
               n_control_cells=int(cells[bank.ctrl_code]),
               target_knockdown=dict(
                   n_with_readout=int(np.isfinite(kd[0][:-1]).sum()),
                   median_residual=float(np.nanmedian(kd[0][:-1])),
                   p10=float(np.nanpercentile(kd[0][:-1], 10)),
                   p90=float(np.nanpercentile(kd[0][:-1], 90)),
                   frac_below_0p2=float(np.nanmean(kd[0][:-1] < 0.2))),
               effect_norm_v0=dict(median=float(np.median(n0)),
                                   p10=float(np.percentile(n0, 10)),
                                   p90=float(np.percentile(n0, 90))),
               variants={})
    for v in variants:
        dv, bv, nv = bank.table(v)
        assert list(nv) == list(names)
        B = dv.astype(np.float64)[np.ix_(keep_rows, scored)]
        D = B - A
        nb = np.linalg.norm(B, axis=1)
        nd = np.linalg.norm(D, axis=1)
        proj = np.sum(D * A, 1) / np.maximum(np.sum(A * A, 1), 1e-300)
        # Restrict every statistic to the rows the variant ACTUALLY changed.  Only
        # 33% of gwps targets have their own gene measured well enough for a per-cell
        # knockdown readout, so a median over all rows is a median over rows the arm
        # never touched and reads 1.0 for every knockdown arm.
        ch = nd > 0
        nch = int(ch.sum())
        cosv = np.sum(A * B, 1) / np.maximum(n0 * nb, 1e-300)
        base = v.split("@")[0]
        k = bank.schemes.index(base) if base in bank.schemes else bank.schemes.index("v0")
        eff = wsum[k][:-1] / np.maximum(wsum[bank.schemes.index("v0")][:-1], 1e-300)
        rec = dict(n_rows_changed=nch, frac_rows_changed=float(ch.mean()),
                   b_ctrl_shift=float(np.linalg.norm(bv - b0) /
                                      max(float(np.linalg.norm(b0)), 1e-12)))
        if nch:
            p = proj[ch]
            rec.update(
                effect_norm_median=float(np.median(nb[ch])),
                amplification=float(np.median(nb[ch] / n0[ch])),
                rel_change=float(np.median(nd[ch] / n0[ch])),
                projection_median=float(np.median(p)),
                projection_mean=float(np.mean(p)),
                projection_sem=float(np.std(p, ddof=1) / np.sqrt(len(p))),
                cos_to_v0=float(np.median(cosv[ch])),
                effective_cell_ratio=float(np.median(eff[keep_rows][ch])))
            # Second, independent prediction: a target whose guide worked LESS well has
            # more escaper mass to remove, so its amplification must be LARGER.  The
            # correlation between the projection and the target's own measured pooled
            # residual expression should therefore be POSITIVE for a real knockdown
            # weighting and zero for a permuted one.
            res = kd[0][keep_rows][ch]
            good = np.isfinite(res)
            if good.sum() > 10:
                rec["corr_proj_vs_residual"] = float(
                    np.corrcoef(res[good], p[good])[0, 1])
                rec["n_corr"] = int(good.sum())
        out["variants"][v] = rec
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*",
                    default=["k562_gwps", "orion_hct116", "orion_hek293t"])
    ap.add_argument("--variants", nargs="*", required=True)
    ap.add_argument("--out", default=str(WORK / "diag.json"))
    args = ap.parse_args()
    rows = []
    for d in args.donors:
        r = analyse(d, "restricted", args.variants)
        rows.append(r)
        print(json.dumps(r, indent=2), flush=True)
    Path(args.out).write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
