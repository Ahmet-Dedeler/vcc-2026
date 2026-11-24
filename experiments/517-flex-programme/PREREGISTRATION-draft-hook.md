# DRAFT preregistration — the programme-sourced uncalled-gene hook on the board (517; for the coordinator to adopt)

Written 2026-09-28 by worker 517, BEFORE the 09-30 reads (458d, 458g) and before any hook arm exists. Not armed, not built.
Numbers in [brackets] are filled from `predict517.py` at ~00:30Z 09-30 with the real 458d/458g deltas, before the arm is
built. Everything else is fixed now.

## The arms (one of them, by the 09-30 branch; `RESULTS.md` §6)
- **HK-N** (458d PARTIAL / NULL, or 458g-null with a null 458d): B0 + `VCC_DE_MDE_SIGN=file:$MZ:165c6b78…:0.05:0.05:le1`
  (`$MZ = experiments/517-flex-programme/work/arch1meanz_ne_hook_abc.npz`). The CLEAN4 block is untouched; the MDE add on
  complement genes takes meanz_ne's sign where |pdown − 0.5| ≥ 0.05 on design-strong targets (strength ≥ 2) and CLEAN4's
  (τ 0.05) on the rest.
- **HK-P** (458d POS): B0 + 458g's env (`VCC_SAMELINE=$ARC:6951dbc7…:0`, `VCC_DE_PRIOR_GATE=le1`) + `VCC_DE_MDE_SIGN=file:$MZ:165c6b78…:0.05`
  (the hook follows the block gate).
- Builder: HEAD + `patch_builder496.py` + `patch_builder517.py`; emitter `1658ae2c` (496's). KAs before the arm: KA-1 unset ==
  the base's output sha; KA-496 `clean4:0.05` == 496a's output sha; for HK-P, KA-2 the same env without the hook == 458g's
  output sha (`9929218c…`).
- **Config diff vs its no-hook base: exactly `de_mde_sign` (+ `de_sign_emitter_sha256` if the base predates 496p).** DE-only:
  per-gene totals, ballast and bulk unchanged, so pds and mse must be bit-identical and there is no re-level.
- `[de_mde_sign]` build lines: per context `source file`, n_sign ~6,300–6,600 of ~9,600–10,100 gated genes, weak targets
  ~83–86 (le1).

## Why (H1 + public analogues; `RESULTS.md` §2)
- H1 (analytic, validated path; 3 seeds): the meanz_ne hook on a CLEAN4 block, gated: nmae −0.0183 raw, reach 0, fid 0,
  overall-equiv +0.0051 [+0.0035, +0.0067]; 496's CLEAN4 hook: −0.00146. Not inflated by shared controls (halves517).
- Cross-line within a lab (Xaira; labhook517 / progest517): a same-lab other-line programme beats an other-lab design on
  uncalled genes in every genome-wide case (+0.04 … +0.11 sign accuracy at the hook's coverage), including a pair where the
  same programme HURTS the block (HEK293T → HCT116).
- Real H1 emitter runs (drive517.sh, ≥ 02:00Z 09-29): [nmae, reach, jac, unsuppressible calls per target, vs the analytic].

## Prediction (Δ vs the no-hook base, raw)
- pds, mse: 0 (bit-identical). fid, jac: |Δ| ≤ 0.001 (no calls change; A/B/C's builds leak ~1 call per target).
- nmae: H1-equivalent at the r implied by 458d's read: [r = …, nmae_H1 = …]; board = H1 × T, T in [0.35, 1.76] until 496a
  reads (10-05), then T = 496a's board Δnmae / −0.00146. Central: [...].
- reach: unpriced (small-gate lottery on H1: 0 … −0.009 by r).

## Bins
- **Falsifiers (read first):** pds or mse differs at all → a bug; |Δfid| > 0.001 or |Δjac| > 0.001 → a call leak: diagnose
  before reading nmae (compare `call_set_unsuppressible` with the base's meta).
- **Primary, Δnmae raw:**
  - ≤ −0.0030: **TRANSFERS.** October runs the hook (`DEF_DAY_MDE_SIGN=hookbank:h1_2025_full:0.05:0.05:le1`,
    `DEF_DAY_HOOK_ESTIMATOR=meanz_ne`; for HK-P, without `:le1` since the block gate is on), after its rehearsal.
  - −0.0030 … −0.0008: **PARTIAL.** Ship only if 496a also transferred (the mechanism is then confirmed twice); otherwise
    a second read.
  - > −0.0008: **FAILS.** Closed on the board; October keeps 502's CLEAN4 form at most.
- **Descriptive:** reach reported, not binned (±0.01 is small-gate noise); the HK-N − 496a difference (if both are read on
  the same base) is the source-quality effect at a fixed mechanism.

## Slot
A same-day pair with 496a (10-05: s1 496a CLEAN4 hook, s2 HK-N on the same base) would isolate the source at fixed
delivery; otherwise the first free slot after 09-30. -m from the MODEL-NAMES queue; description empty.
