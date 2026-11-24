# 468 — the strength-gated programme prior on the board: 458g (preregistration, 2026-09-26 17:04Z; committed efce08c2)
Written before 458g is built and before any 09-27…09-30 entry exists.

## The arm
458g = 458d + `VCC_DE_PRIOR_GATE=q33`.
- In each context, a target's design strength is the number of gated genes whose |design row| > 0.2. That is the
  emitter's own `key` (the DE-channel effect row inside the emission gate), and 385's n_real predictor (Spearman 0.48
  with n_real on H1).
- Targets below the context's 33rd percentile of strength keep the incumbent CLEAN4 (q80) block prior. The rest take
  458d's prior: Arc's H1 programme from all public H1 cells (`arch1d_programme_abc.npz`, `6951dbc7…`).
- Nothing else changes: bulk, q80 protection, coarse MDE, zero repair, seeds.
- Builder `693fde15…` = `4feae27d…` + `patch_builder468.py`. Before the arm is built:
  - KA-1: with every new knob unset, it rebuilds 446b byte for byte;
  - KA-2: with VCC_SAMELINE and no gate, it rebuilds 458d byte for byte (`chain468.sh`).
- Expected config diff vs 458d: exactly `de_prior_gate`.
- Name: "the cryostat started a podcast" (renamed 17:15Z from "the petri dish joined a book club", too close to existing
  names). Description: empty.

## Why (467, H1, production emitter, 3 seeds, raw deltas vs CLEAN4 everywhere)
| arm | fid | nmae | reach | overall-equiv |
|---|---|---|---|---|
| Arc's programme on every target (462l) | +0.0198 | −0.0148 | −0.0388 | +0.0067 |
| Arc on the top 67 % by strength, CLEAN4 on the rest (q33) | +0.0161 | −0.0102 | +0.0015 | **+0.0110** |
| q50 | +0.0126 | −0.0076 | +0.0029 | +0.0091 |
| q67 | +0.0067 | −0.0041 | +0.0017 | +0.0051 |
| CLEAN4 on strong, no prior on weak (c4off33) | −0.0021 | −0.0022 | −0.0114 | −0.0032 |

- The gated-minus-ungated gain is +0.0043 overall-equiv, with a paired seed sd of 0.0002.
- The ungated reach cost is four small-gate targets whose few DE genes run against the line's programme (462 §11). The
  gate hands exactly that stratum back to CLEAN4.
- Weak targets want CLEAN4, not the design's own signs (c4off33, and arcoff33 seed 0 at +0.0078 < gated33's +0.0110).
- q33 matches the panel's construction: targets were drawn from an H1 pilot in thirds by DE count (negligible / subtle /
  strong; FINDINGS H.1).

## What H1 cannot tell us
- `H1-BOARD-CALIBRATION.md`: DE sign priors carried 0.07–0.67× of their H1 fid to the board, and a sign prior's reach is
  not priceable on H1. So the transfer of the reach recovery is the open question, and the reason for the slot.
- On A/B/C the design rows come from the donor-transfer model. The strength distribution differs from H1's; the gate uses
  each context's own quantile. The build's `de_prior_gate` audit records the thresholds.

## Prediction (458g − 458d, raw, same day)
- pds and mse identical: the gate touches only the carriers' block prior, and the ballast pins the bulk (on H1 the
  delta was exactly 0.0 for both).
- fid gives back 0–40 % of 458d's Δfid vs 446b (H1: 19 %).
- nmae gives back 0–50 % of 458d's Δnmae vs 446b (H1: 31 %).
- reach recovers 50–100 % of 458d's reach cost vs 446b, if 458d has one. If 458d shows no cost, reach moves ±0.005.
- **Overall +0.0005 … +0.0060**; central about +0.002.

## Noise
The board's ordering channel moves sd 0.0045 overall between digest-identical arms that differ by an emission re-roll.
458g and 458d share every seed and the design. Only the weak third's block signs differ, and their jitter order is the
same. So the pair's difference is mostly the sign change itself, not a re-roll. After the build, `compare458g.py`
checks that strong targets' blocks are byte-identical between the two packages. If they are not, the noise is
re-rolls, and every bin below widens to ±0.009 (2 sd).

## Bins (458g vs 458d; the same day's two entries)
- **Falsifiers:** pds and mse equal to 458d's to ≤ 1e-6 raw, and jac within ±0.0010. Otherwise it is a bug; diagnose it
  before reading anything else.
- **Primary, Δoverall:**
  - **≥ +0.0015:** the gate transfers. October sets `DEF_DAY_PRIOR_GATE=q33` whenever a programme source is set.
  - **+0.0005 … +0.0015:** weakly positive, in the predicted direction. Keep the gate; it is flagged unconfirmed.
  - **−0.0005 … +0.0005:** null. Either A/B/C's programme cost is not concentrated in design-weak targets, or there is no
    cost. October ships the ungated programme (fewer knobs).
  - **≤ −0.0005:** harm. Close the gate.
- **Mechanism (descriptive):** Δreach ≥ +0.003 with Δfid ≤ 0 means the H1 mechanism transfers. Δreach ≤ 0 with Δfid < 0
  means the gate only removes signal.

## October decision, combined with 458d's P (458's bins, vs 446b)
| 458d (P) | 458g vs 458d | October |
|---|---|---|
| positive / partial | ≥ +0.0005 | programme source + gate |
| positive / partial | null or harm | programme source, no gate |
| null / harm | 458g − 446b fid ≥ +0.0015 and overall ≥ 0 | programme source + gate (the gated construction carries, the ungated one does not) |
| null / harm | otherwise | no programme |

## Slot
09-30 s2, beside 458d in s1. 459dr (expected +0.0004, a calibration of 438's H1 axis) moves to a later day. Fallback:
if 458g is not built and gated by the time the 09-30 chain arms, s2 stays 459dr.

## Addendum 18:04Z 09-26 (before 458g exists as a package, before any 09-27…09-30 entry): the three-way split
458g's weak third carries 446b's prior (CLEAN4), and its strong two thirds carry 458d's (Arc).
- `compare458g.py` checks, block by block, whether every 458g block equals 458d's or 446b's. That holds when the
  emitter's draws do not depend on the prior.
- If it does ("exact split"), the three board entries decompose the programme's effect without extra noise:
  - **458g − 446b** is the programme on the design-strong two thirds;
  - **458d − 458g** is the programme on the weak third.
- Both are reported with the read. They are descriptive, not new bins.
- H1's split for comparison (467, raw, mean over targets vs CLEAN4):
  - strong: fid +0.0239, reach +0.0022;
  - weak: fid +0.0113, reach −0.127, from 4 targets whose few DE genes are nearly all down.

## Addendum 18:43Z 09-26 (before any 09-27…09-30 entry): the H1 price on a panel-like subset
Restricted to h1-test's 67 non-essential targets (the 2026 panel is essential-free):
- ungated programme: +0.0015 overall-equiv (fid +0.0138, reach −0.0497);
- gated q33: +0.0078 (fid +0.0103, reach +0.0020). Gated minus ungated: +0.0063 on this subset (+0.0043 on all 85).

This raises the H1-terms expectation for 458g − 458d on A/B/C. It is still subject to the 0.07–0.67× fid transfer and
the unpriced reach transfer. The bins are unchanged. It also makes a near-null 458d (fid up, reach down, overall ≈ 0)
the expected reading, and the combined October table already handles that case: row 3 (458g − 446b carried) applies.

## Addendum 20:13Z 09-26 (before any 09-27…09-30 entry): the H1 evidence's target-sampling interval
A peer review re-bootstrapped 467 over h1-test's targets, and `467/bootstrap467.py` reproduces it:
- the gate's own benefit is +0.0043, 95 % [−0.0034, +0.0132], P(≤ 0) 0.14–0.17;
- without the four reach-loss targets it reverses (−0.0036).

The "paired seed sd 0.0002" in "Why" above is emission noise only. H1 does not establish the gate. The 09-30 pair is the
test: ~250 weak target-contexts with no re-roll noise (exact split). Its generalisation sd to a new panel is about
0.001 overall by the reviewer's scaling. The bins are unchanged.

## Addendum 2026-09-27 07:58Z (before any 09-28 … 09-30 entry exists): a NULL gate read keeps le1 (487 §2)
The peer's robustness review (487 §2) finds `le1` the only gate rule whose H1 value is positive with its CI above 0 on
every panel: +0.0087 A/B/C, +0.0073 Arc-style, +0.0056 efr. An Arc-style October panel carries more design-weak targets
than A/B/C, so the gate matters more there than on the board's panel. The October table changes in one cell:
- G = NULL (|458g − 458d| < 0.0005): **October keeps `DEF_DAY_PRIOR_GATE=le1`** with the programme source. It no longer
  ships the ungated programme.
- G = HARM (≤ −0.0005) still closes the gate. POS is unchanged.
`read468.py` prints the amended NULL line.
