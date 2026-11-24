# 453 — the October path, run blind on A/B/C (in progress)

## Setup known answers (2026-09-26 01:05–01:20Z)
Surface: `def_day.sh` (HEAD `b517b6d2`, unmodified) on 180's `def-bundle`, the A/B/C controls and panel relabelled D/E/F.

**Attempt 1 stopped at the driver's whitener guard, as designed.** The October fit on this panel reproduced the pinned
A/B/C operator byte for byte (`e2f65ecf…`, the runbook's known-answer row). On a non-A/B/C bundle the driver takes
that as "the panel did not change" and aborts. On the real day, that guard is what catches a wrong bundle.
- Side finding: the whitener fit is deterministic. 452's note that refits differ by eigensolve noise was wrong (see
  452's correction, 01:25Z). The rehearsal dumps differed only in their recorded env string.

**Attempt 2 stopped at the key stage.** `plolo.Operator` recognised the file as the registered A/B/C operator and
refused the day's metadata json (its hash is not the registry's).

**Attempt 3 (from 01:14Z).** The same arrays, re-saved compressed (`66610be7…`, array-equal to `e2f65ecf`), give the
day's operator its own identity, which is exactly the day-of situation. Its json records `resaved_from_sha256`.
- Setup: every stage rc 0.
  - census: HCT116 290, HEK293T 294, CD4 union 291 / complete 276, all equal to the A/B/C census;
  - identity: D/E/F = Jurkat/HeLa/CAL33.
- **The day's key refit equals the pinned A/B/C key byte for byte** (`pert_scale_def_lo0.3_hi3.npz` sha256 `d7f76a22…`,
  all five arrays equal to 410's `wshare_lo0.3_hi3`).
- So on A/B/C's panel the day-of fitted artifacts (operator, key, CD4 banks, census) are the A/B/C ones exactly.
- What remains different from 446b is the preregistered list:
  - DEF composition rules (carrying CD4 sizing at g 0.344, the V5_abc basis);
  - the October stack H1 ×2 / CD4 ×1;
  - the NTC-proxy re-level;
  - the same-line stage, a no-op here: D = Jurkat, 0 panel rows, `noprog`.

## Builds (attempt 3)
| build | factors D / E / F (= A / B / C) | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 (01:40Z) | 1 / 1 / 1 | 0.8407 | 0.7703 | 0.7346 |
| lvl1 (02:05Z) | 1.0793 / 1.1275 / 1.1546 (ntclevel) | 0.9683 | 0.9622 | 0.9611 |
| **lvl2 (02:30Z)** | **1.0859 / 1.1383 / 1.1662** (secant; local exponents 1.851 / 1.854 / 1.869) | **0.9793** | **0.9794** | **0.9795** |

**CONVERGED after 2 levelled builds** (target ρ* 0.9793). The final build is lvl2. The October path's factors sit 1.3 / 1.2 / 0.5 %
below 446b's tuned A 1.1004 / B 1.1517 / C 1.1720: the NTC-proxy target is lower than 446b's matched ρ, and the CD4 ×1 stack
has less spread.

- The first step asks for D 1.0793 / E 1.1275 / F 1.1546. 446b's tuned factors were A 1.1004 / B 1.1517 / C 1.1720, so
  the October path sits ~2 % lower in scale.
- Two causes, as preregistered: the CD4 ×1 stack has less spread to level, and ρ* 0.9793 is below 446b's matched
  0.99 / 0.98 / 0.99.
