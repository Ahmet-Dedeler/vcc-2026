# 534 Part C: do the three small base steps hold off the board? Preregistration

Written 2026-09-29 by the session "Public data scope [2e83f3]", at the coordinator's offer (its message of ~14:40Z), before
any Part C number exists. Scope is 534 §5 item 3.

## 1. Question

Three steps sit inside every October build's base and have never been re-priced as separate steps on the October-like
panels:
- **normrel −1.5** (312: +0.0016 on the board, with a 3/3 PLOLO read at the time);
- **the cap form** of it (319: +0.0009, chosen from the board's own 312/312b dose curve, with no held-out read);
- **the CD4 partial rows at g 0.5** (319p / 319p2: +0.0012 + 0.0003).

Together they are +0.0037 on A/B/C. Is any of them A/B/C-specific in AGENTS 11.1c's sense?

## 2. Instrument, panels, designs

- **Instrument:** 438's PLOLO pricer. A copy of `lib438.py` with every cache and `WORK` path redirected under
  `experiments/534-generalization-audit/work/partC/` (as Part B does). Levelled route, k = 1, Δoverall in board units.
- **Panels:** 438's four: `abc`, `real`, `efp`, `efr`. Each has three cancer folds (K562, HCT116, HEK293T held out in
  turn) with its existing held-out operators, read only. **Reserve R6 is never used here** (its access rule allows
  C10a and C10b only).
- **Designs,** all at the October bulk configuration (`h2c1`, LODO key `V1_5_3`, q 1.75):

| id | normrel | partial rows | prices |
|---|---|---|---|
| B0 | −1.5, cap (production) | g 0.5 (production) | known answer |
| B1 | −1.5, plain (no cap) | g 0.5 | the cap step = B0 − B1 |
| B2 | off (exponent 0: equal per-target donor weights) | g 0.5 | the normrel step with its cap = B0 − B2 |
| B3 | −1.5, cap | off (g 0) | the partial-rows step = B0 − B3 |

## 3. Known answers (before any B1–B3 number; abort on failure)

- The toggles set to production values reproduce the unpatched copy bit for bit on one fold of each panel.
- B0 reproduces 438's stored value for the same design on every panel, to |d| ≤ 1e-12. If `h2c1|V1_5_3+sh1.75` is not
  stored, reproduce the nearest stored design (e.g. `h2c1|V1_5_3+sh1.5`) with the toggles at production.
- B2 at exponent 0 equals the uniform-weight composition, where one is stored (for example in 420, 438 or 415);
  otherwise say so.

## 4. Reading rules (frozen now)

1. **Step value** Δ_s = B0 − B_s, per panel as the mean over its three folds, and pooled over the three October-like
   panels (`real`, `efp`, `efr`; 9 folds) by a paired target bootstrap: 4,000 draws, seed 534.
2. **Verdict** on the October-like pool:
   - HOLDS if Δ_s > 0 and the 95 % lower bound > 0;
   - FRAGILE if Δ_s > 0 and the lower bound ≤ 0;
   - FAILS if Δ_s ≤ 0.
3. **The `abc` row** is reported beside the pool. If its sign disagrees with the board's positive read for that step,
   the instrument disagrees with the board for this class. Report it; draw no conclusion.
4. **Consequence.**
   - A step that FAILS on the October-like pool while its board read was positive is A/B/C-specific. It goes to the
     coordinator as a candidate to switch off for October (a `def_day` change: its call, needing a pass-through and a
     rehearsal).
   - A FRAGILE step is flagged.
   - Nothing in production changes here.

## 5. Prior expectation (written before any number)

| step | expected | why |
|---|---|---|
| normrel | HOLDS or FRAGILE | it had a 3/3 PLOLO read before its board read |
| cap | genuinely uncertain, FRAGILE or FAILS | it was chosen from the board's dose curve alone |
| partial rows | FRAGILE (small, either sign) | — |

## 6. Memory, time, execution

- **Gate:** Part B's preregistration §9, unchanged.
- **Heavy scoring:** ≤ 9 GB, strictly serial. No operator fits are needed.
- **Windows:** no heavy process starts after 2026-09-29 23:45Z, during 2026-09-30 00:30–06:00Z, or during
  2026-10-05 02:00–06:00Z. It waits for any running builder.
- **When:** after Part B's `RESULTS-B.md` exists. Target: done before 2026-10-03 00:00Z, so the result can inform the
  recipe freeze after the 10-05 read.
- **Who:** a resumed run of Part B's worker session (it already knows lib438), with a short mandate pointing here.
  A fresh worker if the resume fails.
