# 527 — the programme-sourced uncalled-gene hook on the board: arm 446bj2zmz (preregistration, 2026-09-29 10:44Z)

Adopts 517's draft (`experiments/517-flex-programme/PREREGISTRATION-draft-hook.md`, arm **HK-N**) as written, with the base
fixed to today's best. Written before any 09-30 entry exists (458d/458g unscored), before 496a's read (10-05), and before the
arm is built. Slot: the first free one, 10-07 or later, fixed by that day's rule.

## The arm
**446bj2zmz** = 446bj2z (entry `be90smZUDcIV7j5MgQ03`, 0.2728260) + `VCC_DE_MDE_SIGN=file:<517 meanz_ne hook artifact>:165c6b78…:0.05:0.05:le1`.
The CLEAN4 block prior is untouched. On design-strong targets (strength ≥ 2 by 468's measure) the coarse MDE add on
complement genes (not head, not block) takes the sign of Arc's H1 programme estimated as the mean z over Arc's non-essential
rows (517 §2.6) where |pdown − 0.5| ≥ 0.05; on the rest it takes CLEAN4's (τ 0.05), i.e. 496's hook. DE-only: per-gene totals,
ballast and bulk are unchanged, so no re-level.
- Builder: HEAD f5d5b92f + `patch_builder517.py`; emitter 1658ae2c (496's). `build527.sh`: builder suite 56/0, 517's
  26 parse tests, **KA-1** (unset == 446bj2z's output `ec3a5ae6`), **KA-496** (`clean4:0.05` == 496a's output `e21ff752`).
- **Config diff vs 446bj2z: exactly `de_mde_sign` + `de_sign_emitter_sha256`** (the emitter pin moved with 496p).

## Why
H1, real production emitter, 3 seeds (517's `drive517.sh`, 2026-09-29): nmae **−0.01829 / −0.01826 / −0.01841** raw, fid
+0.0002, jac −0.0001, reach −0.0015 / −0.0026 / −0.0076, pds/mse exact; overall-equivalent +0.0037 … +0.0049. 496's CLEAN4
hook on the same instrument: nmae −0.00147. The hook's bar is the design's own sign accuracy on uncalled genes (~0.53), lower
than the block's CLEAN4 bar, so it can transfer where the block prior (458d) does not (517 §2.4–2.5).

## Prediction (Δ vs 446bj2z, raw)
- pds, mse: **0** (bit-identical). fid, jac: |Δ| ≤ 0.001.
- nmae: H1 −0.0183 × attenuation × transfer. Attenuation from the 09-30 read via
  `517/predict517.py --fid-458d <Δ> --fid-458g-vs-458d <Δ> --nmae-458d <Δ>` (run at the read, recorded in BOARD-RESULT);
  transfer T in [0.35, 1.76] until 496a reads on 10-05, then T = 496a's board Δnmae / −0.00146. Range before either read:
  −0.002 … −0.032.
- reach: unpriced (small-gate lottery; H1 0 … −0.008 by seed).

## Bins
- **Falsifiers (read first):** pds or mse differs at all → a bug, stop. |Δfid| > 0.001 or |Δjac| > 0.001 → a call leak:
  compare `call_set_unsuppressible` with 446bj2z's meta before reading nmae.
- **Primary, Δnmae raw:**
  - **≤ −0.0030: TRANSFERS.** October adds the hook (`DEF_DAY_MDE_SIGN=hookbank:h1_2025_full:0.05:0.05:le1`,
    `DEF_DAY_HOOK_ESTIMATOR=meanz_ne`, via 517's `patch_defday517.py` + `patch_sameline517.py`, an inert KA and a rehearsal).
  - **−0.0030 … −0.0008: PARTIAL.** Adopt only if 496a also transferred (the mechanism confirmed twice); otherwise a second
    read.
  - **> −0.0008: FAILS.** Closed on the board; October keeps 502's CLEAN4 form at most.
- Descriptive: reach; and (446bj2zmz − 446bj2z) − (496a − 446bj2) = the source effect at fixed delivery (zs and the hook act on
  different genes; the additivity is an assumption, stated).
