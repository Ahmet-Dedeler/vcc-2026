# 387 — Preregistration: sign repair to the un-whitened consensus (`VCC_SIGN_REPAIR=cons`, experiment 173) on the 319p2 composition. Frozen 2026-09-13 ~06:30Z, before the build.

**The arm.** 319p2's exact configuration (`VCC_NORMREL_EXP=-1.5`, `VCC_NORMREL_MODE=cap`, `VCC_CD4_PARTIAL_G=0.5`,
`VCC_NORMREL_PM` unset) plus `VCC_SIGN_REPAIR=cons`: on the ~16 % of gated (target, gene) entries where whitening
FLIPPED the sign relative to the un-whitened donor consensus, the sign is set back to the consensus and the magnitude
kept (builder `sign_repair_to_consensus`, applied to the whitened matrix before the CD4 residual; inert unless set).
One config key: `sign_repair` absent → `cons`. Builder `d8dc1613…`, driver `27edc486…` (hash-gated).
Package `submissions/387-cap15-partial50-signcons.vcc`, model name "the lint roller has opinions".
Reference for Δ: **319p2, `lqWNXX5YaYArdglu6Stn`, 0.2201159** (same panel, anchors r4).

**Why this slot buys information.** The 2026-09-13 board shows every top-20 team ahead of us on `nmae` (ours 0.935
raw; pack 0.905–0.927; chill 0.835) and the pack ahead on `fid`/`mse`/`reach`, while we lead on `pds`/`jac`: our
stack is tuned for retrieval and pays on the four DE-side members. The scorer computes the DE table (LFC, sign,
p-value) from per-cell arithmetic CPM and `pds`/`mse` from bulk sums; our emitter welds every cell to one library
size, which ties the two channels together. The plan under test today (lab 388/389) is to unweld them so the DE
channel can carry a better sign/magnitude estimator without paying the `pds` tax that killed every DE-side arm of
the last two weeks. The cheapest estimator is the un-whitened consensus (173: +0.014 pooled sign accuracy on the
reference DE pool; H1 six-member net ≈ +0.002; PLOLO says the `pds` cost is ~6× larger on real lines). **This arm
reads, on the real contexts, the two halves separately: the DE-member gain of consensus signs (`score_nmae`,
`score_fid`, `score_reach`) and the `pds` cost (`score_pds`).** The decoupled channel would keep the first and avoid
the second, so the first number is the quantity the plan's expected value depends on.

## Primary reading (member split, in overall points = Δscore_member/6)
| quantity | H1 prediction (173) | reading |
|---|---:|---|
| DE gain = Δnmae + Δfid + Δreach | +0.0012 +0.0018 −0.0015 (reach unresolved) | ≥ +0.002 → consensus signs carry real DE information on A/B/C; the decoupled channel's prize is at least this. 0 ± 0.002 → the sign gain is inside the board's DE noise; the channel's prize is amplitude only. ≤ −0.002 → consensus signs hurt on real lines; the DE channel must not use them. |
| pds cost = Δpds | −0.0016 (H1) … −0.010 (PLOLO 6×) | records the exchange rate of partial un-whitening on real lines; informs how much of 319's whitening the bulk channel should keep. |
| mse | ≈ +0.002 (cos +7.8 % via the cliff mechanics) | secondary. |

## Bins (Δ overall vs 0.2201159), for the record
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | consensus signs net-positive even with the `pds` cost; `cons` joins the composition; re-price the decoupled channel with the measured DE gain. |
| **−0.0005 … +0.0005** | flat; the two halves cancel — read the member split, not the total. |
| **Δ ≤ −0.0005** | the `pds` cost dominates as PLOLO predicted; `cons` does not join the bulk; the DE gain (if positive) is taken only through the decoupled channel. |

**Stated in advance.** The displayed row is expected to FALL (net ≈ −0.006 … +0.003); it is not corrected by
re-submission. Upload 13–17 min. Second slot of the day: not planned (nothing else preregistered; the emitter of
389 cannot be ready before 00:00Z).

## Addendum 2026-09-13 ~06:55Z — written before the 387 entry exists; the bins above are unchanged
388 (H1, analytic, known answers exact) now predicts the consensus-sign DE half NEGATIVE, not positive: the
un-whitened consensus (E2) loses gate sign accuracy −0.009 and reach −0.045 against the delivered signs, and 173's
flip (E3) loses reach −0.045. So the expected member split is DE ≤ 0 with the pds cost on top, i.e. the lower bin.
The arm is still spent because the slot expires unused otherwise and the real-line pds price of partial un-whitening
has never been read. The "decoupled channel with consensus signs" plan it served is closed by 388; the unwelded
channel's live use is the reference's direction base rate on the low-confidence call block (390: +0.0081 raw fid on H1
with exact delivery), which 389's emitter prices next.
