# 701 — K-UP: the DE call set one step LARGER than shipped (preregistration)

Written 2026-10-14 ~22:40Z by the science session [bf003b]. This is before the arm is built, before its H1 pricing runs, and
before any 10-15 entry exists (10-15 s1 = DN, 671, still building). Slot: 10-15 s2, coordinator [400b6e]'s assignment.

## The arm (one variable vs 446bj2z, entry `be90smZUDcIV7j5MgQ03`, 0.2728260)
- **446bj2zku30** = 446bj2z's exact env + `VCC_RC_CALLSET_CPM=3000`, so every target's call set is `CallSet(400, 3000)`:
  K 3,400 (shipped: `CallSet(400, 2200)`, K 2,600). The knob was patched in 523; its validated range is [100, 5000].
- Built by `run_arm_generic.sh` against 446bj2z's meta, exactly as 523's K16/K20 were (`523/chain523.sh` step 3), on the
  builder as committed after 671's KAu. Expected config diff vs 446bj2z: exactly `rc_callset_cpm` (+ `de_sign_emitter_sha256`
  if the emitter pin differs from 446bj2z's recorded one; that has no output effect, per 499 F2).
- Pass 1 keeps every gene total, so the bulk is untouched: pds and mse must equal 446bj2z's to <= 1e-6 raw. DE-only arms move
  the jackknife mse by O(1e-8), per 496a.
- The base is 446bj2z, not 630arm, on purpose: K16 / K20 / 446bj2z are on it, so this read completes a 4-point dose curve
  on one base.

## Why
- fid is our worst member against the field (scaled 0.0055 vs 0.04–0.09 for #2–#20 on 10-14). F_p = k / max(n_pred, n_real),
  where k counts sign matches on ALL adjudicable called genes. On a target with n_real > n_pred, every added call with a
  right sign pays 1/n_real (coverage). On a target with n_real < n_pred, an added call dilutes F toward its own sign accuracy.
- The board's dose points (vs 446bj2z, Δ(fid + jac) in overall units): K 1,600 (CPM 1,200) **−0.00234**; K 2,000 (CPM 1,600)
  **−0.00086**. Both smaller sizes lost. The fid fall carries the loss; jac is flat (K20's jac was +0.0004 scaled).
- Two shapes fit those three points equally well within the read noise (residuals 0.0002–0.0004 against σ ≈ 0.00035 per
  difference), and they disagree about the next step:
  - linear, 0.00023 per 100 CPM, predicts **+0.0018** at CPM 3,000;
  - quadratic (exact through the 3 points), peak at CPM ≈ 2,216, predicts **−0.0014** at CPM 3,000.

  This read decides between them. Either answer changes the October K: grow it, or confirm 2,200 as the optimum.
- 537's real-emitter H1 (3 seeds) has the same concave shape, larger: K 2,000 −0.0068 / +0.0011 / +0.0030 at f 1 / 0.5 / 0.35.
  The board sits near f 0.87 for SIZE changes (K16 read, FINDINGS §0). I-662's single f 0.61 failed for size.

## Class and noise
- **Scorer mechanics.** No biology changes. It carries to D/E/F only because the final panel was drawn by the same design
  (negligible / subtle / strong bins), from the same lab, at the same depth. The value would be re-read on D/E/F's own
  feedback if members come back.
- **Same emission law, pass 2 re-drawn.** A K change consumes `pass2_rng` in proportion to the designed set. The re-roll sd is
  ≈ 0.00025 overall (537 §2.5, confirmed on 12 board pairs), so σ = 0.00025 below.

## Primary statistic
**Δ(fid + jac) in overall units, vs 446bj2z** = (scaled Δfid + scaled Δjac) / 6, as in 523. Reach and nmae are reported,
and the total Δoverall is a guard on adoption.

## Predictions (written before the H1 runs; an amendment adds the H1 numbers before the package exists)
| source | Δ(fid + jac) overall at CPM 3,000 |
|---|---|
| board curve, linear | +0.0018 |
| board curve, quadratic | −0.0014 |
| 537 H1 at f 1 / 0.87 / 0.61 / 0.5 (701 runs, ×0.72 board convention) | TBD by the amendment below, before the build |

## Bins (σ = 0.00025 overall)
- **POS, ≥ +0.0010 (≥ 4σ):** the board keeps rewarding size above 2,200. Adopt CPM 3,000 for October if Δoverall (all six
  members) is ≥ 0, via a `DEF_DAY_RC_CALLSET_CPM` pass-through with an inert-unset KA. [17adab] offered to write it, 690a
  style. Probe CPM 4,000 next only if the read is ≥ +0.0015.
- **WEAK, +0.0005 … +0.0010 (2–4σ):** adopt only if Δoverall ≥ +0.0005 too; no further probe.
- **NULL, −0.0005 … +0.0005:** the curve is flat above 2,200. Keep 2,200 (`rc_callset_cpm` off).
- **NEG, ≤ −0.0005:** an inverted U with its peak in [2,000, 3,000]. Keep 2,200, and record that the quadratic shape holds.
- **Falsifiers (stop and diagnose before reading):**
  - pds or mse moves by more than 1e-6 raw: a bulk leak;
  - fid falls AND jac rises: the wrong direction for a larger set, i.e. a construction error;
  - nmae moves by more than ±0.002 raw: the delivered lfc on gate genes changed more than a size change implies.
- Deliverability (from the build's `VCC_RC_AUDIT_DUMP`, if one is written, and from the H1 runs): prescribed vs realized calls
  per target. If realized calls grow by less than half the prescribed +800, the arm under-delivers, and NULL/NEG reads as
  "capacity-bound", not as "size does not pay".

## Tries
- Board, on this switch: 3 arms (K16, K20, this one). Best-of-N bias does not apply to a dose curve read as a shape.
- H1: 537 ran 2 sizes × 3 seeds; 701 runs 3 sizes (CPM 2,800 / 3,000 / 3,400) × 3 seeds. Only CPM 3,000 is built.

## What it does not test
Per-target sizing (n_pred(t) = predicted n_real(t)). 192 put its oracle at +0.036 overall, and 521 closed it by payoff
(n_real predicted at r 0.19–0.43 on held-out cancer lines). It reopens only with a better n_real predictor; 700 measures
PIE's.

## AMENDMENT 1 (2026-10-15T21:21Z, [d80a10] science; written before any 701 build or package): the H1 price was not run
- run701.sh v2 waited for R7a's wrapper and was capped at 10-15 12:30Z. R7a never finished: the Mac slept at 05:32Z and
  shut down at 12:55Z. The cap fired with no run started (`work/run701.log`: "R7a not finished by the cap: stopping").
- So the "537 H1 at f 1 / 0.87 / 0.61 / 0.5" row of the prediction table stays EMPTY. The board-curve predictions
  (linear +0.0018, quadratic −0.0014), the bins, the falsifiers and the deliverability check stand unchanged. They never
  depended on the H1 row.
- The H1 runs may still be done later as a post-hoc instrument check. Their numbers will be labelled post hoc and
  decide nothing.
