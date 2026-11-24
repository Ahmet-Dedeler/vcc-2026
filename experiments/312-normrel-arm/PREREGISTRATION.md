# 312 — Preregistration, written BEFORE the entry exists

**Frozen 2026-09-07, before `vcc submit`.** A board reading is only evidence if the prediction predates the entry
(AGENTS §11.1; the definition of *board-validated* in FINDINGS §1). This file is that prediction.

## The arm
`VCC_NORMREL_EXP = -1.5` on the A/B/C challenge panel, everything else at the shipped defaults
(`VCC_ALPHA` 0.982, `VCC_CD4_TARGET_F` 0.34, `VCC_H1_UNION_G` 0.28, `VCC_CD4_DEFLATE_K` 5). Per-target donor
weights w_j(t) ∝ (‖row_j(t)‖ / that donor's own panel median)^−1.5, renormalised over the donors carrying t, the
row then rescaled to the uniform-weight row's norm (**shape-only**). Incumbent: `cd4svd5`
`Xv0ayvuoHpsKYWv2iu42`, **0.2160944021991198**.

## What the local instruments say
- PLOLO (308, three cancer folds): Δ`pds` **+0.00259 [+0.00000, +0.00527]**, **3/3 folds positive**; `nmae` leg
  **−0.00095** (a benefit); priced **+0.00104 pessimistic** = 0.55× the +0.0089 gate.
- 097 route (a) a6, an independent earlier measurement on two surfaces: **+0.0008…+0.0009 overall**, "the only arm
  positive on both", direction-only leg 1.046/1.039, magnitude-only leg 1.001/0.999 (nil).

## The prediction, in bins fixed now
Let Δ = (this entry's `score_avg`) − 0.2160944021991198.

| bin | reading |
|---|---|
| **Δ ≥ +0.0040** | the pessimistic pricing is too harsh by ≥ 4×. **Several closed arms reopen** and must be re-priced before Oct 22 — this is the outcome that would change the programme. |
| **+0.0010 ≤ Δ < +0.0040** | pricing is roughly calibrated; the arm is real and small, as both instruments said. Keep it; it is not a route to first place. |
| **−0.0010 < Δ < +0.0010** | null on the board. The local `pds` gain does not survive the six-member composite. |
| **Δ ≤ −0.0010** | the arm LOSES on the board despite 3/3 local folds and an independent 097 replication. That would be the important negative: PLOLO mis-signs this entire class of per-target weighting, and every arm ever priced through it is suspect. |

**Secondary, read from the raw members, no verdict:** the `nmae` leg against its predicted −0.00095; raw
`expr_mse_unbiased_capped_norm` → implied transfer cosine via 197's inversion under §4.1c's validity guard
(ours today 0.09728), to see whether a shape-only per-target reweighting moves direction at all.

## Stated in advance
- Expected magnitude is **small**; this does not approach the +0.0233 gap to `cqawesome` and will not move us off
  3rd. The slot is spent because it is **prize-irrelevant and expires tonight**, and because no per-target
  weighting has ever been read on the real scorer.
- **The board shows the LATEST entry, not the best.** If Δ < 0 the displayed row falls. That is expected, costs
  nothing that is scored, and the entry is NOT to be re-submitted to restore rank (the standing rule).
- `cd4svd5`'s full record is already captured at `data/board-entries/Xv0ayvuoHpsKYWv2iu42.json`, so the incumbent
  row survives this entry replacing it in the API.
