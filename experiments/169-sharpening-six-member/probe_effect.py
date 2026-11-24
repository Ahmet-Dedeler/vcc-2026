"""169 — recover the h1-test artifact's log-space EFFECT from its emitted counts, and prove it.

Sharpening is a transform of the per-gene predicted EFFECT.  Every existing h1-test code path
(093/095/157/164) re-emits from the built prediction's per-gene panel totals, which holds the
profile FIXED -- exactly the axis sharpening has to move.  So the effect has to be recovered.

Production renders
    profile = clip(expm1(b + h), 0)/TS,  then / its own sum        `profile_from_effect`
with `b = log1p(TS * ctrl_totals / ctrl_totals.sum())`             `make_submission_pm_cd4:2056`
and `emit_counts(mode="quasi")` sets `totals = floor(n*L*profile)` `donor_transfer.py:380`.

So `n*L*profile` lies in `[totals, totals+1)`, the unbiased point estimate is `totals + 0.5`
(this is also `realize164.py`'s convention), and

    h = log1p(TS * Z * p) - b,     p = (totals + 0.5) / sum(totals + 0.5)

for one unknown scalar `Z` per target -- the renderer's own normaliser.

TWO THINGS THAT COULD FAIL, AND WHAT THEY WOULD READ IF THE PREMISE WERE WRONG
-----------------------------------------------------------------------------
1.  `profile_from_effect` sets the effect to EXACTLY ZERO outside the CPM gate.  So off-gate
    the recovered `h` must be zero for every target -- and in particular its ACROSS-TARGET
    SPREAD must be at the quantisation floor, because a zeroed gene has the same profile for
    every target.  A wrong `b` (mean-per-cell CPM instead of pseudobulk, or the wrong cell
    set) shows up as a large shared offset; a wrong GATE shows up as across-target spread in
    the mis-assigned band, which no scalar can absorb.
2.  The gate threshold itself is therefore MEASURABLE, not assumed.  This is run over
    `min_cpm in {2, 5, 10}` and the answer is read off the off-gate spread.

Both are reported below.  The measured answer -- the 095 h1-test artifact was built at
`--min-cpm 5.0`, not production's 10.0 -- independently reproduces experiment 164 §3.1, which
established the same fact from a completely different statistic (the between-target spread of
the predicted pseudobulk, 12.14 against a flooring bound of 1.0).

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
TS = 50_000.0
CANDIDATES = (2.0, 5.0, 10.0)


def recover(sums, b, expression, gate):
    """(Z, H) for every target, `Z` fitted on the off-gate block where `h` is known to be 0."""
    off = ~gate
    Zs, H = [], []
    for k in range(sums.shape[0]):
        t = sums[k] + 0.5
        p = t / t.sum()
        m = off & (sums[k] >= 1)        # exclude never-emitted genes: +0.5 there is pure prior
        Z = float(expression[m].sum() / (TS * p[m].sum()))
        Zs.append(Z)
        H.append(np.log1p(TS * Z * p) - b)
    return np.asarray(Zs), np.stack(H)


def main() -> None:
    z = np.load(os.path.join(WORK, "h1test-totals.npz"), allow_pickle=True)
    pg = np.asarray([str(g) for g in z["pred_genes"]])
    hg = np.asarray([str(g) for g in z["harness_genes"]])
    assert np.array_equal(pg, hg), "prediction and harness gene axes differ"
    ct = z["ctrl_totals"].astype(np.float64)
    tar = np.asarray([str(t) for t in z["pred_targets"]])
    allsums = z["pred_sums"].astype(np.float64)
    ncell = z["pred_ncell"]

    b = np.log1p(TS * ct / ct.sum())
    expression = np.expm1(b)
    cpm = 1e6 * expression / expression.sum()
    sel = np.flatnonzero(ncell == 400)
    sums = allsums[sel]
    tsel = tar[sel]
    print(f"genes {pg.size} | sum expm1(b) = {expression.sum():.4f} (TS = {TS:.0f}) | "
          f"{sel.size} perturbation groups, {sums.sum(1).mean() / 400:.0f} counts/cell")

    # ---- 1. which gate did this artifact ship? --------------------------------------
    print(f"\n{'min_cpm':>8}{'gated':>8}{'Z median':>11}{'off-gate |mean|':>17}"
          f"{'off-gate SPREAD':>17}{'on-gate spread':>16}")
    print("  (off-gate spread is the discriminator: a zeroed gene has the SAME profile for")
    print("   every target, so its across-target sd must sit at the quantisation floor)")
    best = None
    for mc in CANDIDATES:
        gate = (cpm > mc) & (b > 0)
        Zs, H = recover(sums, b, expression, gate)
        band = (~gate) & (sums.min(0) >= 20)      # off-gate but well-counted: no floor excuse
        sd_on = float(np.median(H[:, gate].std(0)))
        if int(band.sum()) < 100:
            print(f"{mc:>8.1f}{int(gate.sum()):>8}{np.median(Zs):>11.5f}"
                  f"{'band too small':>17}{'n/a':>17}{sd_on:>16.5f}")
            continue
        sd_off = float(np.median(H[:, band].std(0)))
        mu_off = float(np.median(np.abs(H[:, band].mean(0))))
        print(f"{mc:>8.1f}{int(gate.sum()):>8}{np.median(Zs):>11.5f}{mu_off:>17.5f}"
              f"{sd_off:>17.5f}{sd_on:>16.5f}")
        if best is None or sd_off < best[1]:
            best = (mc, sd_off, Zs, H, gate)
    MIN_CPM, sd_off, Zs, H, gate = best
    print(f"\n  --> the 095 h1-test artifact's gate is `--min-cpm {MIN_CPM}` "
          f"(off-gate across-target sd {sd_off:.5f}, {CANDIDATES}); 164 §3.1 reached the same"
          f" conclusion from a different statistic")

    # ---- 2. quality of the recovery -------------------------------------------------
    off = ~gate
    q, qg = np.abs(H[:, off]), np.abs(H[:, gate])
    qexp_on = 0.5 / np.maximum(sums[:, gate], 1.0)
    print(f"\nZ: median {np.median(Zs):.6f}  min {Zs.min():.6f}  max {Zs.max():.6f}")
    print(f"recovered |h| off-gate: median {np.median(q):.3e}  q99 {np.quantile(q, .99):.3e}")
    print(f"recovered |h| ON-gate : median {np.median(qg):.4f}  q90 "
          f"{np.quantile(qg, .9):.4f}  mean {qg.mean():.4f}")
    print(f"quantisation floor on gated genes (0.5/n): median {np.median(qexp_on):.3e}  "
          f"--> median SNR {np.median(qg) / np.median(qexp_on):.0f}x")
    mg = H[:, gate].mean(0)
    print(f"common mode: on-gate across-target |mean| median {np.median(np.abs(mg)):.4f} "
          f"= {np.median(np.abs(mg)) / np.median(qg) * 100:.1f}% of median |h| "
          f"(096 measured (A-B)/A = 0.84%, i.e. 9.2% of the row norm)")

    np.savez_compressed(os.path.join(WORK, "h1test-effect.npz"),
                        targets=tsel, genes=pg, b=b, gate=gate, cpm=cpm,
                        min_cpm=np.float64(MIN_CPM), H=H.astype(np.float64), Z=Zs,
                        sums=sums.astype(np.int64), ctrl_totals=ct)
    print(f"\nwrote {os.path.join(WORK, 'h1test-effect.npz')}")


if __name__ == "__main__":
    main()
