# 446 — q 1.75 levelled (the next dose on the sharpening curve). Preregistration, 2026-09-24 20:55Z
Written before 433h's and 433k's entries exist (09-25 slots) and before any 446 entry exists.

## The arm
446 = 433k with `VCC_BULK_SHARPEN_Q` 1.5 → 1.75 and q 1.75's own re-level. The config diff vs 433k is exactly
`bulk_sharpen_q` + `bulk_ctx_scale`. The DE channel is identical in configuration: carriers, MDE coarse, zero repair.
- **Re-level at matched ρ.** The factors put each context at the same NTC-proxy ρ as 433k's bulk (measured on 417k's
  restored artifact: ρ 0.9909 / 0.9828 / 0.9873, cliffs from `ntclevel417.py`, KA PASS). On 446u, the unlevelled q 1.75
  measurement build, that gives √(ρ_417k / ρ_446u) = **A 1.0927 / B 1.1406 / C 1.1601**. The only difference from 433k is
  the q step at the same position against each context's cliff.
- The built artifact's ρ is re-measured before arming. If it is more than 0.5 % off 433k's in any context, 446 is rebuilt.

## Price (446 − 433k)
- PLOLO on the current bulk (444, realised KMAX 16, stack|K3), q 1.5 → 1.75:
  - cancer folds: rel cos +3.1 … +3.6 %, β-emulated pds +0.005 … +0.012;
  - kolf2: +0.8 % / +0.0018.
- 420's levelled q steps (cancer folds, k = 1): 1.5 → 1.75 ≈ +0.005 overall.
- Board transfer of shape steps is k_sh ≈ 0.42–0.65 of the cancer folds (417a), and kolf2 predicted 417a's board.
- **Prediction: overall +0.001 … +0.004; pds +0.001 … +0.005 raw; mse −0.003 … −0.012 raw.**

## Bins (vs 433k's board row)
- **pds (raw):** Δ ≥ +0.002 means the dose step transfers; +0.0005 ≤ Δ < +0.002 means partial; |Δ| < 0.0005 means flat;
  Δ ≤ −0.001 means q 1.75 overshoots, so q stays at 1.5.
- **mse (raw):** Δ ≤ −0.005 transfers; −0.005 < Δ ≤ −0.001 partial; Δ > 0 means overshoot or mis-level. Re-read
  ntclevel's ρ.
- **DE members:** within the ballast's realisation noise (nmae ≤ 0.001, fid ≤ 0.001, jac ≤ 0.001, reach ≤ 0.01). A DE
  member outside these bins means the render change leaked into the carriers: stop and diagnose.
- **Overall:** ≥ +0.001 with pds ≥ 0 means q 1.75 becomes the bulk default for A/B/C and October (SHIPPED_ENV q, and
  the day-of re-level runs on it). |Δ| < 0.001 means flat: keep q 1.5, which is the smaller extrapolation for October.
  ≤ −0.001 means close the q curve at 1.5.

## Conditions on 433k's read (tonight)
446 goes on the board only if 433k's q step (1.25 → 1.5, levelled) reads pds Δ ≥ +0.001 raw or mse Δ ≤ −0.004 raw
against 433h. If that step was flat, 446 is not submitted and the q curve closes at 1.5.

## Addendum 21:05Z (still before any 09-25/26 entry): the October view, from 438's saved grid
438's panel-agnostic pricer already priced the stack at q 1.75, levelled (the D/E/F default), against the stack at q 1.25
(`work/summary438_abc_real_efp.md`):
| levelled q step vs q 1.25 | abc | real | efp |
|---|---|---|---|
| q 1.5, k = 1 | +0.00847 3/3 | +0.00492 3/3 | +0.00412 3/3 |
| **q 1.75, k = 1** | **+0.01311 3/3** | **+0.00746 3/3** | **+0.00607 3/3** |
| q 2, k = 1 | +0.01486 3/3 | +0.00781 3/3 | +0.00655 3/3 |
| q 1.5, 417a-calibrated (k 0.65 / k_pds 0.42) | +0.00436 3/3 | +0.00221 3/3 | +0.00189 3/3 |
| **q 1.75, calibrated** | **+0.00648 3/3** | **+0.00294 3/3** | **+0.00258 3/3** |
| q 2, calibrated | +0.00680 3/3 | +0.00245 3/3 | +0.00242 3/3 |
Levelled q 1.75 beats q 1.5 on every October-like panel at both transfer rates. The calibrated increment is +0.0021
(A/B/C), +0.0007 (`real`) and +0.0007 (`efp`); q 2 adds nothing over 1.75 at calibrated transfer. Per fold, the levelled
q 1.5 → 1.75 step raises β-emulated pds on all 13 panel×fold cells (+0.005 … +0.012). Cosine falls on `real` (−0.6 …
−1.0 %) and on two `efp` folds (−1.2, −3.2 %), so the October value is carried by pds. If 446 reads in its transfer bins,
the October SHIPPED q moves to 1.75, and the day-of `ntclevel417.py` re-level runs on it.

## Correction 21:30Z (the preregistered re-measurement fired; still before any 09-25/26 entry)
`ntclevel417.py` on 446's restored artifact read ρ **0.9780 / 0.9652 / 0.9685** against the target (433k via 417k)
0.9909 / 0.9828 / 0.9873: 0.65–0.95 % low in scale units, above the 0.5 % rebuild threshold. The render is compressive:
from 446u to 446, ρ grew as factor^1.854 / ^1.862 / ^1.872, not factor². The corrected factors
(ρ_target / ρ_446u)^(1/exponent) are **A 1.1004 / B 1.1517 / C 1.1720**. The arm is rebuilt as **446b** with them, and
its riders become 446bj / 446bz. 446's package (`c90bec90…`) was deleted unsubmitted; its meta stays. Every bin above
applies to 446b.

## 446b verified 22:07Z (still before any 09-25/26 entry)
Package `submissions/446b-…-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc`, sha256 `e171cb39…`.
- Config diff vs 433k is exactly `bulk_ctx_scale` + `bulk_sharpen_q`.
- 0 zero-total gated genes; novelty rc 0.
- Render energy vs 410b 0.9603 / 1.0282 / 1.0785 (433k: 0.9600 / 1.0278 / 1.0778).
- **ntclevel ρ 0.9908 / 0.9827 / 0.9873 against the target 0.9909 / 0.9828 / 0.9873**: matched to 1e-4.
- Riders 446bj (LODO key), 446bz (zs) and 446bd (DLD-1 union 2) are building on 446b's bulk.
