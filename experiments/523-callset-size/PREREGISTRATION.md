# 523 — the DE call set's size on the board: a dose pair (preregistration, 2026-09-28 15:10Z)

Written before the builder knob exists, before either arm is built, and before any 09-29 entry exists. Slot day: fixed
later by the day's rule (10-06 or after), never before its rule is written.

## The arms (both one variable vs 446bj2z, entry `be90smZUDcIV7j5MgQ03`, 0.2728260)
- **446bj2zk16** = 446bj2z + `VCC_RC_CALLSET_CPM=1200`: every target's call set is `CallSet(400, 1200)`, K 1,600 (today
  `CallSet(400, 2200)`, K 2,600).
- **446bj2zk20** = 446bj2z + `VCC_RC_CALLSET_CPM=1600`: K 2,000.
- Built after `patch_builder523.py` (the knob; unset = 2200, bit-identical, proved by KA0 = 446bj2z's output `ec3a5ae6…`
  rebuilt byte for byte with the knob unset). Expected config diff vs 446bj2z: exactly `rc_callset_cpm` (+
  `de_sign_emitter_sha256` when 496p's emitter pin differs from 446bj2z's recorded `edb5b8dc`: 499 F2, no output effect).
- The bulk is untouched (pass 1 keeps every gene total): pds and mse must equal 446bj2z's to ≤ 1e-6 raw.

## Why (521 §5, §8.1)
The global call-set size is regime-sensitive under today's emission. On H1 (overall units): K 1,600 is −0.0155 at
f 1.0, −0.0003 at f 0.5, **+0.0048 [+0.0009, +0.0100]** at f 0.35 (fid +0.0036, jac +0.0029, reach +0.0077) and +0.0061 in
the compressed regime; K 2,000 is −0.0079 / +0.0007 / +0.0029 / +0.0035. The board's n_real regime is disputed: 192's
fid-coverage fit says f ≈ 0.5 (→ ≈ 0), the board's jac level says f ≈ 0.33 (→ positive). The pair reads the regime
directly, and the answer transfers to D/E/F sizing (same panel design, cancer lines).

## Noise model (why fid + jac is the primary read)
Changing K changes the emitter's pass-2 placement draws (`place_k` consumes `pass2_rng` in proportion to the designed set),
so reach carries the ordering channel's re-roll spread (sd ≈ 0.027 scaled reach ≈ 0.0045 overall between digest-identical
arms, 230). fid and jac are functions of the designed call set and the deterministic sign/magnitude design; their
re-roll spread is small. The primary statistic is **Δ(fid + jac) in overall units** (scaled Δfid + Δjac, / 6).

## Predictions (board, vs 446bj2z; fid/jac at 192's 0.72× retrodiction, reach unpriced)
| arm | f 0.35 regime | f 0.5 regime | f 1.0 regime |
|---|---|---|---|
| K 1,600: Δ(fid + jac) overall | **+0.0024** | ≈ −0.001 … 0 | clearly negative |
| K 2,000: Δ(fid + jac) overall | +0.0014 | ≈ 0 | negative |

## Bins (primary = Δ(fid + jac) overall of K 1,600 vs 446bj2z)
- **≥ +0.0012: the board sits in the small-n_real regime (f ≲ 0.4).** Adopt K 1,600 for A/B/C DE arms; for October,
  write a `DEF_DAY_RC_CALLSET_CPM` pass-through (a 502-style patch with an inert KA) and set 1,200 in the recipe.
- **−0.0005 … +0.0012: f ≈ 0.5 (or the gain is too small to act on).** Keep 2,200; the pair's K 2,000 arm refines only if it
  is outside ±0.0005.
- **≤ −0.0005: the board is closer to H1's regime than either estimate says.** Keep 2,200; record 192's fit as the better
  regime estimate.
- K 2,000 is the dose check: it must sit between 446bj2z and K 1,600 (monotone). A non-monotone pair means the primary read
  is noise-dominated: no adoption on it.
- **Falsifiers:** pds or mse moving by more than 1e-6 raw (a bulk leak: stop); jac moving the wrong way while fid rises by
  more than the prediction (a call-set construction error: diagnose); nmae moving by more than ±0.002 raw (the delivered lfc
  of dropped genes changed more than the design implies: diagnose before reading).
- Reach and nmae are reported, not binned.

## Amendment (2026-09-30T21:10Z, before any 10-06 entry exists): the dose check
The rule "K 2,000 is the dose check: it must sit between 446bj2z and K 1,600 (monotone). A non-monotone pair means the primary
read is noise-dominated: no adoption on it" is **replaced**. Two measurements made after it was written:
- On the real production emitter (worker 537, `experiments/537-arm-pipeline/RESULTS.md` §3.2 P2/P3, H1), Δ(fid + jac)
  overall vs 446bj2z at f 1 / 0.5 / 0.35 / cmp is **K 2,000: −0.0068 / +0.0011 / +0.0030 / +0.0033** (2 seeds) and
  **K 1,600: −0.0132 / +0.0002 / +0.0028 / +0.0033** (3 seeds). K 2,000 ≥ K 1,600 in every regime: an optimum at or above
  K 2,000 is H1's prediction, not noise.
- The board's DE read noise on this emission law is ~0.00025 overall per read (537 §2; 547; 548), far below the bins' gaps.
**New dose rule (applied when both arms are scored):** adopt the K whose Δ(fid + jac) vs 446bj2z is larger, if that value is
≥ +0.0012; if the two are within ±0.0005 of each other, adopt K 2,000 (it loses about half as much in the H1-like regime);
if one is ≥ +0.0012 and the other ≤ −0.0012, the pair is outside H1's shape: diagnose before adopting anything. The primary
bins for K 1,600 alone, the falsifiers and the October consequences above are unchanged (for K 2,000 the October knob value is
1,600, not 1,200). Tries on this switch so far: 0 board arms; H1: 2 sizes × 2–3 seeds (537).
