# 472 — the absolute prior gate: `VCC_DE_PRIOR_GATE=le<N>` (queued 20:16Z 09-26)

**Why (peer review, docs/REVIEW-2026-09-26-second-opinion.md §2).** 468's share form `q33` cuts at `np.quantile` of an
integer strength (count of gated genes with |design row| > 0.2), so the cut jumps between classes.
- H1: threshold 1.72, so weak = strength ≤ 1 (28/85).
- A/B/C: threshold 2, so weak = {0, 1} (86/83/83). This is what 458g tests.
- The efr synthetic October panel: 54–58 % of targets have strength ≤ 1, so the threshold drops to 1 and weak = {0}
  only. The 94–100 strength-1 targets per context keep the programme. That is the class behind all of H1's reach loss:
  CNN3, COL6A1, ZNF664 and AASS all have strength 1.

`le1` names the rule the evidence supports. On A/B/C it selects q33's exact weak sets, so it must rebuild 458g byte for
byte. The 09-30 read of 458g then validates the rule October runs.

**Plan (`chain472.sh`, after 470b):**
1. Patch the committed builder (`patch_builder472.py`). q-mode is untouched.
2. Tests: suite 56/0, 468, 472 (14 parse cases), 451, 417.
3. KA-le1: le1 with 458g's env reproduces 458g's output sha (`9929218c…`).
4. KA-1: unset reproduces 446b's (`9559923a…`).

Pre-checked on a copy of builder 693fde15: the patch applies, compiles, and passes all 14 of 472's parse tests and all
of 468's.

**October:** if 458g reads positive, `DEF_DAY_PRIOR_GATE=le1`, not `q33`. `def_day.sh` passes the string through
unchanged.

**Open:** is ~55 % CLEAN4 right on a panel whose designs are mostly weak (efr)? On H1, strength ≤ 1 targets have median
n_real 9. A design-weak target on an efr-like panel may be weak for a different reason: design shrinkage in a new line.

## Results (22:46Z 09-26)
- Patched: builder `4b8aa34e` → `efb8c9dc`. Suite 56 passed / 0 failed / 1 skipped; 468, 472 (14 cases), 451 and 417
  tests pass.
- **KA-le1 PASS**: `VCC_DE_PRIOR_GATE=le1` with 458g's env rebuilds 458g byte for byte (output `9929218c6f361824…`).
  - Gate lines: "le 1 threshold 2; programme on 214 / 217 / 217 targets, CLEAN4 on 86 / 83 / 83", the same as q33's.
  - So **09-30's read of 458g validates the absolute rule October will run.**
- **KA-1 PASS** (23:16Z): unset rebuilds 446b byte for byte (`9559923a…`). Builder `efb8c9dc` committed.
