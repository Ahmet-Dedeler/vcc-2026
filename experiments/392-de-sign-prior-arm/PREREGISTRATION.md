# 392 — Preregistration: CPM-block signs from the 4-bank cross-line prior, DE channel only (`VCC_DE_SIGN_PRIOR=clean4`) on 319p2. Frozen 2026-09-13 ~07:20Z, during the build, before any entry exists.

**The arm.** 319p2's exact configuration (`VCC_NORMREL_EXP=-1.5`, `VCC_NORMREL_MODE=cap`, `VCC_CD4_PARTIAL_G=0.5`) plus
`VCC_DE_SIGN_PRIOR=clean4`. Builder `260e03aa…` (knob inert unset), emitter `src/vcc2026/rc_emit_k.py` `a12cf8f9…`
(389's K-class emitter + 391's prior mode), prior `clean4_prior_challenge.npz` `7e532f84…` (mean P_down over K562-GWPS,
K562-essential, Nadig Jurkat, Nadig HepG2 from 371; no challenge or H1 data). Config diff vs 319p2: exactly the five
`de_sign_*` keys. Package `submissions/392-cap15-partial50-designclean4.vcc`, model name "the stapler went on sabbatical".
**Reference for Δ: 319p2, `lqWNXX5YaYArdglu6Stn`, 0.2201159** (raw pds 0.85474 / mse 0.80990 / nmae 0.93480 / fid 0.50789 /
reach 0.27780 / jac 0.04090). Second slot of 2026-09-13 UTC, submitted only after 387's entry is terminal.

**Mechanism.** cell-eval2 computes the DE table from per-cell ARITHMETIC CPM and pds/mse from BULK sums. The package emits
three equal library classes per context (A/B/C control-library quartiles), keeps every per-gene panel total, keeps the
|z| design, and steers the arithmetic mean of each target's 2,200 CPM-block calls to the prior's sign at the delivered
magnitude, down-signed magnitudes shrunk so the block's arithmetic mass is unchanged.

## Evidence before the entry (all H1 development panel unless marked)
| source | fid raw | reach raw | nmae raw (lower better) | pds / mse |
|---|---|---|---|---|
| 390 analytic, exact KA | +0.0116 [+0.0057, +0.0177] | +0.018 [−0.017, +0.056] | +0.0024 [−0.0011, +0.0065] | 0 / 0 |
| **391 end to end, real emitter + scorer** | **+0.0122 [+0.0067, +0.0182]** | +0.0175 [−0.019, +0.055] | +0.0024 [−0.0011, +0.0066] | **0 / 0 exactly** |
| 391 same prior deranged across genes (same dose) | −0.0019 | −0.0267 | +0.0046 (harm) | 0 / 0 |
Prior accuracy by held-out line (371, donor-delta signs): Jurkat +0.020 over all-down, HepG2 +0.016, K562 +0.012, H1
+0.007 (weakest), HCT116 +0.001, debiased HEK293T −0.004. A/B/C audit of 319p2 (20 targets per context): block delivered
down 0.50 / 0.50 / 0.53, prior flips 53 / 52 / 53 %, prior coverage 87 / 79 / 77 % — the same shape as H1.

## Predictions in board units (overall points per raw unit from 920 entries: fid +0.561, reach +0.188, nmae −0.255)
H1 1:1 transfer: fid +0.0069, reach +0.0033, nmae −0.0006, pds 0, mse 0, jac 0 → **Δ ≈ +0.0095 overall.**

## Falsifiers (checked first; any failure voids the reading as a mechanism test)
1. `score_pds` and `score_mse` equal 319p2's to ≤ 1e-9 (per-gene totals identical by construction). A difference means the
   build changed the bulk and the DE reading is confounded.
2. `de_wilcoxon_sig_jaccard` within ±0.0005 raw of 319p2 (call set unchanged by construction).

## Primary reading: Δ score_fid (scaled), predicted +0.041
| bin | reading |
|---|---|
| ≥ +0.020 | the prior's per-gene signs carry A/B/C information and the channel delivers it; build the next variants on this channel. |
| −0.010 … +0.020 | no transfer at this construction (H1-specific or HEK-like lines); the channel stays, the prior does not. |
| ≤ −0.010 | the prior is anti-informative on these lines; close CLEAN4 for A/B/C. |

## Overall bins (Δ score_avg vs 0.2201159)
| bin | reading |
|---|---|
| **≥ +0.004** | transfers at ≥ 40 % of H1; 392 becomes the incumbent. |
| +0.0005 … +0.004 | partial transfer; incumbent if fid is in its top bin, else read the members before any dose step. |
| −0.0005 … +0.0005 | flat. |
| ≤ −0.0005 | harm; close at this construction. |
Secondary, no verdict: reach (board twin noise sd ≈ 0.0045 overall), nmae (expected −0.0006 overall).

**Stated in advance.** The displayed row becomes this entry's score and is not corrected. Upload 13–17 min.

## Addendum 2026-09-13 ~07:47Z — build diagnostics, written after the build and BEFORE the entry exists; bins unchanged
Real A/B/C build (900 blocks, meta `the-stapler-went-on-sabbatical.h5ad.meta.json`, config diff PASS = the five de_sign_* keys):
feasible dose λ median **0.901** (min 0.087, max 1.0) vs H1's median 1.000 (87 % of targets at 1.0); targeted flips per
block median 1,169 (986–1,281) of 2,200; steering error ≤ 1e-15; clips ≤ 5 genes per block; block mass residual ≤ 1.4e-11;
call set delivered median 2,599 (Jaccard vs design median 0.9996); per-cell totals 17,695–22,322.
Production `sign_flip_frac` diagnostic: 319p2 {'max': 0.39339309094571456, 'median': 0.3480954910761901, 'min': 0.22543112403906088} → 392 median 0.360 (0.249–0.412).
Consequence for the reading: a block with λ ≤ 0.5 delivers no sign flips and λ in (0.5, 1) delivers them at magnitude
(2λ − 1)·|E0|, so the delivered fid effect is expected at or somewhat below H1's +0.0122 raw; the primary prediction and all
bins stand as frozen. If score_fid lands in the middle bin, check first whether it is dose rather than transfer (rebuild at a
lower class margin or with 394's q80 variant, which moves less mass per block).
