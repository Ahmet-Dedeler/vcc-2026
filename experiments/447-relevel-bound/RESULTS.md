# 447 — the re-level bound widened to [0.7, 1.4] (2026-09-25 01:20Z, coordinator)

**Why.** 438's own level file needs these `VCC_BULK_CTX_SCALE` factors on the essential-free, realistic-coverage October
draw (`efr`): 1.145 / 1.172 / **1.208** at q 1.5, **1.212 / 1.213 / 1.261** at q 1.75, 1.29 / 1.26 / 1.32 at q 2. The
builder refused anything outside [0.8, 1.2] at import. So on an efr-like D/E/F panel the day-of re-level could not reach
the cliff at the q the board now favours (433k: q 1.5 transferred; 446b: q 1.75 armed for 09-26). A/B/C needs ≤ 1.19
(446b: 1.1004 / 1.1517 / 1.1720; 446bd2: up to 1.1885).

**Change.** `scripts/make_submission_pm_cd4.py`: the one import-time check moves from [0.8, 1.2] to [0.7, 1.4], with a
comment. Nothing else in the builder moved; every in-range value flows through unchanged, so every existing arm rebuilds
byte-identically by construction. `experiments/417-bulk-sharpen/test_next417.py`: its refusal message string updated.

**Tests (work/):**
- `builder-tests-before.log` and `builder-tests-after.log`: both **56 passed, 0 failed, 1 skipped**.
- 417/410 unit tests pass: `test_next417`, `test_sharpen417`, `test_bulkio417`, `test_carrier417`, `test_ballast410`,
  `test_pert_scale`.
- `bound-accept.log`: A=1.3,B=1.35,C=1.4 accepted; 0.7 / 1.2 / 1.4 accepted; the efr q 1.75 factors
  1.2119 / 1.2130 / 1.2607 accepted; 1.41 and 0.69 refused with the new message; unset → None.
- Builder sha **63dfd7d0 → b6c3f131**.

**Still open for the day:**
- Above 1.2, render compression (ρ ∝ f^1.86 on A/B/C) and ballast capping (KMAX 16) grow. So the day-of loop is:
  unlevelled build → `ntclevel417.py` → levelled build → `ntclevel417.py` again, and a second correction if ρ is off by
  more than 0.5 %.
- The October rehearsal on `ef-bundle-r` must exercise a factor above 1.2 end to end before Oct 22.
- Waiters and chains that pin the builder by hash must use b6c3f131 from now on. The 09-26 waiters submit prebuilt
  packages and are unaffected.
