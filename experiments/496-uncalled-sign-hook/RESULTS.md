# 496 — the uncalled-gene sign hook on H1: nmae exactly as 478 priced it (−0.0015), reach noise-dominated; WEAK (rider only)

2026-09-27 05:56Z, coordinator. Preregistered in `PREREGISTRATION.md` (25bac194) before any result.

**Change.** The emitter's MDE add on complement genes (not head, not block) takes the CLEAN4 per-gene prior sign where
`|P_down − 0.5| ≥ τ` (`mde_add_sign`; `patch_emitter496.py` on a local copy of the production emitter `edb5b8dc`,
`k_emit496.py`). Harness: 435's H1 harness (`run496.py`), coarse MDE (0.1, 0.1, 0.3), vs 435's production-emitter runs at
the same seeds. The production prior file's P_down equals 478's CLEAN4 P on all 6,346 covered genes (max |d| 0); the
harness maps it by gene name (h1-test's axis is not the official order). 2,550 of 10,779 gated genes carry a sign at
τ 0.05 (1,499 down); 810 flips per target (478: 810).

**Known answer: PASS.** The hook emitter without a sign source reproduces 435's production coarse run (seed 0), stats
bit for bit.

| run | nmae | reach | fid | jac | pds, mse | overall-equiv [target 95 % CI] | unsuppressible / target | clip frac | fd |
|---|---:|---:|---:|---:|---|---|---|---|---|
| τ 0.05 seed 0 | −0.00147 | +0.00193 | +0.00007 | −0.00005 | 0 (exact) | +0.00079 [−0.0006, +0.0025] | 5.21 → 3.12 | 0.286 → 0.273 | 0.790 → 0.771 |
| τ 0.05 seed 1 | −0.00143 | +0.00069 | −0.00001 | −0.00005 | 0 | +0.00049 [−0.0010, +0.0023] | 5.31 → 3.07 | 0.286 → 0.273 | 0.790 → 0.771 |
| τ 0.05 seed 2 | −0.00147 | −0.00588 | −0.00003 | −0.00007 | 0 | −0.00075 [−0.0043, +0.0014] | 5.12 → 3.08 | 0.286 → 0.274 | 0.790 → 0.772 |
| **τ 0.05 pooled** | **−0.00146** | **−0.00109** | +0.00001 | −0.00006 | 0 | **+0.00018 [−0.0017, +0.0020], P(≤0) 0.41** | | | |
| τ 0 seed 0 | −0.00240 | −0.00229 | −0.00002 | −0.00009 | 0 | +0.00018 [−0.0033, +0.0035] | 5.21 → 2.02 | 0.286 → 0.264 | 0.790 → 0.763 |

**Against the bins.** Falsifier clean (pds/mse exactly unchanged); no leak (fewer unsuppressible calls, lower clip
fraction). Pooled overall +0.00018 < +0.0005 → **WEAK: a rider only, never its own slot.**

**What the real run adds to 478's analysis.**
- nmae is exactly what the analytic path priced (−0.0015 at τ 0.05, −0.0024 at τ 0), with seed sd 0.00002: the delivery
  is clean and deterministic.
- reach is not. 478 held the call order fixed and found +0.002 … +0.0035; the real |z| realisation moves the order and
  reach lands at +0.0019 / +0.0007 / −0.0059 by seed (τ 0.05) and −0.0023 (τ 0). On H1 the reach term is noise-sized and
  on average slightly negative.
- Priced on nmae alone (the only part H1 prices reliably): +0.0004 overall at H1 slopes, +0.0007 at the MDE's 1.76×
  board transfer. The reach term is ±0.001 on this instrument.

**Status.** The production patches are written and dry-checked, not applied: `patch_emitter496.py` (emitter pin →
`1658ae2c…`) and `patch_builder496.py` (`VCC_DE_MDE_SIGN=clean4:<tau>`, config key `de_mde_sign`; the new emitter pin
changes every later build's `de_sign_emitter_sha256`, so later arms' expected diffs vs a pre-496 reference include it).
It rides only on an arm whose DE members are otherwise pinned (a DE-only change on a bulk-identical base), and only
after the programme-source read (09-30), since the complement source could then be Arc's programme (478b).

## Addendum (06:30Z): where the reach swings come from
Per target on H1 (τ 0.05 vs 435's base, same seeds), reach moves on only 3–4 of 82 targets per seed, and two
small-gate targets carry all of it:
- **AASS** (n_real 17, one of 467/462's four reach-loss targets; base reach 0.75–0.88): −0.125 / −0.062 / −0.625 by seed;
- **STXBP1** (base reach 0.57): +0.286 / +0.286 / +0.143.
- ZNF32 −0.167 on seed 1; the rest move by < 0.003.

The between-seed correlation of per-target deltas is positive (0.85 / 0.60 / 0.37), so the same targets move. But
their size depends on a few sign changes near the 90 % purity edge of tiny gates. On H1's 85 targets that is a two-target
lottery. On the board's ~900 (target, context) pairs it should average toward its mean, which H1 cannot price.
The reading stands: the hook is an nmae lever (−0.0015 raw, seed sd 0.00002) with an unpriced, small-gate reach term.

## Addendum (06:48Z): the hook stacks exactly on the zs rider (496z)
The zs pass (`dn_k250_zs`, 423/429; on the board 09-28 as 446bj2z) was run on top of the hook (τ 0.05) for 3 seeds, against
435's coarse + zs runs (`stack496z.py`, queue `work/queue496z.txt`, ≤ 7.4 GB).
| seed | hook on the zs base: nmae / reach / overall-equiv | hook on the plain base (496) |
|---|---|---|
| 0 | −0.00150 / +0.00193 / +0.00079 | −0.00147 / +0.00193 / +0.00079 |
| 1 | −0.00144 / +0.00069 / +0.00050 | −0.00143 / +0.00069 / +0.00049 |
| 2 | −0.00144 / −0.00588 / −0.00076 | −0.00147 / −0.00588 / −0.00075 |
- nmae −0.00146 (sd 0.00003) on the zs base against −0.00146 on the plain base: additive. The exact-U pass leaves the hook's
  reach untouched to the digit, and pds/mse stay bit-identical.
- Determinism: re-running the hook's pre-pass blocks (needed for the zs tag) reproduced the first run's stats bit for bit on
  all three seeds.
- So in October the hook and zs can be priced as a sum: zs −0.0018 + hook −0.0015 nmae raw on H1.
