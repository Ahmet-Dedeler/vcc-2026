# 483 — Testh's "cross-target geometry" lever on H1: negative for our pipeline (analysis only)

Peer-review session, 2026-09-26 ~21:05Z. No build, no board, no production edit.

**Result.** Projecting shared cross-target directions out of our predicted bulk effects lowers h1-test pds at every
setting. Testh's disclosed setting (k 14, β 0.5) gives Δpds −0.0045 with a panel-target donor basis and −0.0024 with a
full-screen donor basis, against random-basis controls at 0 ± 0.0007. Squared error barely moves (−0.4 % … +0.05 %).
Our whitening by the donor signal covariance (γ 0.65) plausibly already removes much of what the projection removes.
That is not the whole story: Testh's 09-07 description says it whitens too ("v6-LOCO transfer (whiten, …)"), so its
+0.0094 board pds (09-15) may come from its basis (LOCO truth deltas of the panel's targets in six contexts, see
Caveats) or from the rest of its pipeline.

## What was tested

- Testh's description (board entry `i2OFyop7G5rc2UFCv1nw`, 09-15): `D <- D − 0.50·(D V14ᵀ) V14`, V14 = top-14
  right-singular vectors of the stacked row-normalised truth deltas of six leave-one-context-out panel folds.
- Here: D = our predicted bulk effects on h1-test, `log1p(TS·sums/Σ) − b_ctrl` from 435's harness sums
  (`stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz`, TS 5e4), truth from 323's `reference.npz`.
- Bases: (a) the same 85 targets' rows in K562-GWPS / Orion HCT116 / Orion HEK293T (82/84/84 targets; no H1 data),
  row-normalised; (b) 3,000 random non-evaluated rows from each full screen; (c) controls: a random orthonormal basis of
  the same rank, and our own predictions' top-k directions (self-deflation).
- pds implemented from the vcc2026 brief §1 (85 target genes removed; ties split).
- **Known answer:** base pds 0.8952380952 reproduces the harness's metrics json exactly.

## Numbers (Δpds vs base, h1-test, seed 0; wins/losses = targets whose PDS rose/fell)

| basis | k | β 0.25 | β 0.5 | β 0.75 | β 1.0 |
|---|---:|---:|---:|---:|---:|
| panel-target donor rows | 5 | −0.0003 | −0.0004 | −0.0041 | −0.0106 |
| | 14 | −0.0020 | **−0.0045** (13/14) | −0.0112 | −0.0217 |
| | 30 | −0.0025 | −0.0060 | −0.0156 | −0.0284 |
| full-screen donor rows | 5 | 0.0000 | 0.0000 | | |
| | 14 | −0.0013 | **−0.0024** (18/12) | | |
| | 30 | −0.0004 | −0.0006 | | |
| random orthonormal (control) | 14 | 0.0000 | +0.0003 | +0.0003 | +0.0006 |
| our own top-k (self) | 14 | −0.0049 | −0.0136 | −0.0263 | −0.0524 |

Full grid in `geom_h1.json` (`geom_h1.py`); the full-screen variant is `geom_h1_full.py`.

## Caveats

- One harness seed, 85 targets; H1 prices direction changes with the right sign in the calibration record (387: H1
  −0.0043, board −0.0098), but this is not a board read.
- The full-screen variant's SVD peaked at 10.8 GB RSS (LAPACK workspace for 9,000 × 18,080); the panel-target variant
  at 0.45 GB. Use the panel-target form, or a randomized SVD, if anyone re-runs it.

Closed at this construction: rank-k projection of donor-shared directions from our whitened bulk. A different basis
(for example LOCO truth deltas of the panel's own targets across six contexts, exactly as Testh built it) is not
tested here, but the full-screen and panel-target donor bases agree in sign.
