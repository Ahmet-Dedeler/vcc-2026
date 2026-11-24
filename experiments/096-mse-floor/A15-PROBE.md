# a14 / a15 — pre-submission `n_pred` probe

Run 2026-08-29 21:14 UTC on the real challenge controls, before either slot was spent.
`probe_de.py` recomputes the scorer's own Wilcoxon call set against each context's real control
cells; no ground truth is used or needed. 12 perturbations per context, median `n_pred`.

| ctx | a12 (board `0.167987`) | a14 | a15 | a15 vs a14 | a15 vs a12 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 4649.5 | 4667.0 | 4640.0 | −0.58% | −0.20% |
| B | 4498.0 | 4481.5 | 4458.5 | −0.51% | −0.88% |
| C | 4171.0 | 4171.5 | 4161.5 | −0.24% | −0.23% |
| **mean** | **4439.5** | **4440.0** | **4420.0** | **−0.45%** | **−0.44%** |

Call set direction balance, % UP: a14 `60.8 / 63.1 / 63.3`, a15 `60.8 / 63.3 / 63.3`.

## What this settles

**099 §8 falsifier 4 does not fire.** It would have fired at `n_pred` moving more than 5%,
which would mean the post-whitening CD4 residual is changing the call set and that `fid`/`reach`
are therefore unpredicted rather than the modelled `~0`. The measured move is **−0.45%**, an
order of magnitude inside the threshold, and the UP fraction is flat to a tenth of a point. The
`fid`/`reach` ≈ 0 assumption in 099 §8's prediction table is sound for this arm.

**a14's own call-set claim is confirmed too.** A14-PREDICTION asserted "`fid`, `reach`, `jac`
~0 — no call-set policy change". a14 moves `n_pred` **+0.01%** against a12. The amplitude
change from `VCC_ALPHA` 1.0108 → 0.982 does not disturb the call set at all.

This is the check A12-BOARD-RESULT demanded ("**Probe its `n_pred` before submitting**") and
that a13 was never given. For contrast, a13 — the arm that lost `−0.0203` — sat at
`4544.5 / 4406.0 / 3971.5`, a **−2.3%** mean move, five times a15's.

Raw records: `work/probe-a14.json`, `work/probe-a15.json`.

## What it does NOT settle

`n_pred` is a property of our own prediction. It shows the call set is stable; it says nothing
about whether the genes in it have the right sign, which is where `nmae`/`reach`/`fid` actually
live and where both harnesses under-called a13 by 8x. Falsifiers 1–3 of 099 §8 — `nmae` raw
rising more than +0.004, `score_pds` falling at all, `score_mse` improving by less than +0.008 —
remain open and are only readable on the board.
