# 433 — the count-balanced MDE on the board. Preregistration, 2026-09-23 09:40Z (build running; no entry exists)

**Arm.** 433a = 417e + `VCC_DE_MDE=0.1,0.1,0.3,coarse`. One variable against 417e: config diff exactly `de_mde`
(`coarse: true`) + `de_sign_emitter_sha256` (edb5b8dc). Builder 2c90d9e1 (commit 5568ee84). Build script
`build433a.sh` (through `experiments/428-zero-gene-audit/run_arm_generic.sh`); the rc dump goes to
`work/rc_audit_433a.jsonl`.

**Why a new arm instead of 426a.** 426a is the same MDE without `coarse`. 429a's per-block dump (431) shows the plain
MDE delivers nothing in context A (`mde_lambda` 0 in 92 % of A blocks; B 0.94, C 0.86). So 426a's board read would be
a B+C read, about 2/3 of the effect, and it is not the version we would ship. 426h's `coarse` delivery gives the same
nmae on H1 (−0.00953 ± 0.00021 vs −0.0096, 3 seeds) with a 17× smaller per-class count shift.

**Switch rule for the 09-24 second slot (fixed now, before the build finishes):**
- 433a replaces 426a if and only if all of these hold:
  - its build gates pass (config diff exactly the two keys, CLIFF equal to 417e's, novelty rc 0, package sha);
  - its dump shows context A's `mde_lambda` median ≥ 0.5;
  - B's and C's medians are not below 429a's (0.936 / 0.855) by more than 0.05.
- Otherwise slot 2 stays 426a, and this file records why.
- The 417e band rule of `arm0924_s2c.sh` is unchanged. If any 417e DE member is outside the band vs 417a, there is no
  automatic slot 2.

**Predictions (433a − 417e, raw board members).** H1 gives nmae −0.00953 at dose 1. The board transfer for pure
magnitude has been between ×1 and ×2.3 of H1.
- **pds: identical** to ≥ 10 significant figures. The bulk is untouched; any difference is a bug.
- **mse:** |Δ| < 2e-4.
- **nmae: −0.009 … −0.022**, central −0.012 (dose ≈ 1 in all three contexts; filled in from the dump below).
- fid, jac: |Δ| ≤ 0.001. reach: −0.004 expected, inside board noise.
- **overall: +0.0016 … +0.0055.**

**Reading bins (nmae), the same as 426a's:**
- **Δnmae ≤ −0.008: the MDE law transfers.** It rides every later arm. Coarse is the delivery for October.
- **−0.008 < Δnmae ≤ −0.003: partial transfer.** A/B/C gate signs are closer to chance than H1's. Keep it if overall
  > 0; no dose increase.
- **Δnmae > −0.003: no transfer.** Close the family and record the mechanism (gate signs too close to chance for
  magnitude to pay).
- **pds not identical: stop.** Diagnose before any further DE arm.

The displayed row may fall: 433a sits on 417e, whose stack has no board read yet.

## Build record (filled in after the build, before any entry)
Written 10:05Z 09-23, before any 09-24 entry exists.
- Built 09:20–09:54Z on builder `2c90d9e1` / emitter `edb5b8dc`.
  - Artifact sha `6396ddee…`, package `433a-…-mde010coarse.vcc` sha `d879cf9f…`.
  - Config diff vs 417e is exactly `de_mde` (`coarse: true`) + `de_sign_emitter_sha256`.
  - CLIFF equals 417e's to 4 decimals (0.9477 / 0.9408 / 0.9447): the bulk is untouched. Novelty rc 0.
- Zero-gene audit: 24 zero-total gated genes, the same defect every package since 410b carries. This arm has no
  repair on purpose: one variable vs 417e.
- **Delivered MDE dose** (`mde_lambda` over blocks, from `work/rc_audit_433a.jsonl`):

  | context | 433a coarse, median / mean | plain MDE (429a), median / mean |
  |---|---|---|
  | A | 0.956 / 0.954 | 0.000 / 0.045 |
  | B | 1.000 / 0.975 | 0.936 / 0.892 |
  | C | 1.000 / 0.991 | 0.855 / 0.834 |

  The mean dose is 0.973, against 0.59 for the plain delivery, so 433a carries ~1.65× 426a's MDE.
- Side effects in the meta:
  - clip_frac median 0.148 (417e 0.032, 426a 0.094);
  - prescribed calls missed median 9 (417e 2, 426a 5);
  - unsuppressible calls median 0 (max 90).

  So the call design is disturbed a little more than 426a's. H1 put that at fid −0.0002, jac −0.0001 and reach −0.004.
- **Switch rule met** (A median 0.956 ≥ 0.5; B/C at 1.0, above 429a's): slot 2 on 09-24 = **433a**. Waiter
  `arm0924_s2d.sh` replaces `arm0924_s2c.sh` (pid 83836, killed while it was still waiting for 417e's capture).
- **Predictions at the delivered dose** (433a − 417e, raw):
  - pds identical; mse |Δ| < 2e-4;
  - **nmae −0.009 … −0.022, central −0.012**;
  - fid/jac |Δ| ≤ 0.001, reach −0.004;
  - **overall +0.0016 … +0.0055**.

  The bins above stand unchanged.
