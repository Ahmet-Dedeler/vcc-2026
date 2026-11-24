# 452 — rehearsal of `scripts/def_day.sh`, the October day in one gated script

**Status:** 452 and 452b both PASSED end to end (2026-09-25 20:04Z / 21:49Z). Bundle `experiments/238-essential-free-panel/work/ef-bundle-r`; builder
`1dcb29a3` (VCC_SAMELINE); q 1.75; same-line flags `:2:pt:noprog`.

## Attempt 1 (07:26Z / 08:29Z): stopped at the first re-level step, safely
Both rehearsals (452 real identity, 452b forced) stopped right after lvl0 with "relevel_iter printed no factors". The
builder was restored identical both times.
1. **The first step needs ntclevel's own factors.** `relevel_iter.py` requires two builds per context for a secant. After
   lvl0 alone it refused. 448 had taken the first step from ntclevel's EV scale directly, and `def_day.sh` skipped that.
   Fix: `relevel_iter.py` now takes the first step from the json's own `scale` when a context has one build. Known
   answer on 448's data:
   - lvl0 alone → 1.2084 / 1.2562 / 1.2871 (ntclevel's own);
   - lvl0 + lvl1 → 1.2216 / 1.2718 / 1.2452;
   - all four → CONVERGED 1.2218 / 1.2721 / 1.2514.
   All three reproduce.
2. **bash brace-expanded inline `python3 -c "{…, …}"` one-liners.** The rho log line printed empty three times. Fix: the
   helpers moved to `scripts/def_day_helpers.py` (`rho`, `past_cliff`, `same`, `inbound`, `nctx`), each tested under
   bash.

## Attempt 2 (from 18:17Z), real identity file
- Setup rc 0 (whitener 330 s).
- Identity unanimous: D = Jurkat, E = HeLa, F = CAL33.
- Same-line mode ON for D only: Nadig Jurkat, programme off (`noprog`), 9 covered targets blended at w 2 with the
  per-target prior. E and F unchanged.
- Dump `4c2963da…`.

| build | factors D / E / F | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 | 1 / 1 / 1 | 0.6861 | 0.6206 | 0.5911 |
| lvl1 | 1.1948 / 1.2562 / 1.2871 (ntclevel) | 0.9609 | 0.9566 | 1.0568 → artifact deleted (past its cliff) |
| lvl2 | 1.2069 / 1.2718 / 1.2452 (secant) | **0.9795** | **0.9794** | 0.9656 |
| **lvl3** | 1.2069 / 1.2718 / **1.2517** (bracketing secant) | 0.9795 | 0.9794 | **0.9790** |

E and F reproduce 448's lvl0/lvl1 exactly, as they must (same bundle, no blend there). D's lvl0 ρ is 0.6861 against 448's
0.6706: the same-line blend on 9 targets moves the spread +2.3 %.

**Result (20:04Z): PASS on every criterion.**
- CONVERGED after 3 levelled builds at D 1.2069 / E 1.2718 / F 1.2517 (ρ 0.9795 / 0.9794 / 0.9790, target 0.9793). lvl1's
  artifact was deleted automatically because F was past its cliff.
- Zero-gene audit 0 (7 below −5).
- The meta carries q 1.75, the factors, and `sameline.contexts.D` = {covered 9, programme off, per-target prior on, w 2}.
- `vcc prep --dry-run`: 360,000 cells, targets verified. Builder restored to `1dcb29a3` (identical).
- **Wall time 1 h 47 min** from start to a prepped artifact: setup 6 min, dump 12 s, four builds of ~25 min, audits.
- The per-context secant from lvl0/lvl1 put F at 0.9656 on lvl2, closer than 448's fixed-exponent pass (0.9463). One
  bracketing step then converged it.

**Side observation (452b, 20:11Z):** a fresh whitener fit in a new work dir gives a different dump sha each time:
448 `16cdf3d3`, 452 `4c2963da`, 452b `b721c098` (same bundle, same config). The ρ values agree to 4 decimals (452's
lvl0 E/F = 448's), so the differences are float noise from the eigensolve. Consequence for the day: byte-for-byte
rebuilds of a D/E/F artifact need the SAME fitted operator. The driver reuses `<work>/whitener_<id>.npz` when it exists, so
keep the round's work dir and never refit mid-round.

## Attempt 2b (452b, 20:05–21:49Z): forced identity, heavy blending. PASS
Forced calls:
- D = K562: k562_gwps covers 300/300; programme off under `noprog`.
- E = HCT116: orion_hct116 covers 300/300.
- F = A549: 6/300, rows only.
All at w 2 with the per-target prior.

| build | factors D / E / F | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 | 1 / 1 / 1 | 0.8794 | 0.7199 | 0.5934 |
| lvl1 | 1.0553 / 1.1664 / 1.2847 (ntclevel) | **0.9755** | 0.9680 | 1.0554 → artifact deleted |
| lvl2 | 1.0553 / 1.1735 / 1.2436 | 0.9755 | **0.9794** | 0.9662 |
| **lvl3** | 1.0553 / 1.1735 / **1.2498** | 0.9755 | 0.9794 | **0.9790** |

- **The heavy blend raises the un-levelled spread**: D 0.879 vs 0.686 in 452; E 0.720 vs 0.621. So the first step asks for
  much smaller factors (D 1.055 instead of 1.195), all inside the builder bound. D converged on the first step, E on the
  second, F on the third.
- CONVERGED after 3 levelled builds. Zero-gene audit 0 (5 below −5). Prep dry run: 360,000 cells, targets verified.
  Builder restored to `1dcb29a3`. Wall time 1 h 44 min.

## Verdict
`scripts/def_day.sh` is ready for Oct 22. It covers the same-line mode off (452's E/F), light (452's D, 9 rows) and heavy
(452b's D/E, 300 rows). It converged in 3 levelled builds every time, and it deletes past-cliff artifacts as it goes. The
open inputs are the arguments: q, w/flags, and riders (DEF-RUNBOOK "THE DAY IN ONE PAGE").

## CORRECTION (2026-09-26 01:25Z): the day pipeline IS bit-reproducible; the "side observation" above is wrong
- The refitted operators in every rehearsal work dir are byte-identical: `whitener_union_pruned_def.npz` `75773d59…` in
  448/efr, 452/efr, 452/efr_forced and both attempt1 dirs.
- On A/B/C's panel, the October fit reproduces the pinned operator `e2f65ecf…` exactly (453, 01:11Z).
- The dumps' arrays are identical as well: `bulk` and `carriers` max|d| = 0 between 448's `16cdf3d3` and 452's
  `4c2963da`, and likewise `perts` and `genes`.
- The only difference is the `env` string stored in the dump. It records `VCC_BULK_DUMP`, which is the dump's own
  output path, so a dump written into a different work dir hashes differently.
- So the ρ equality to 4 dp was exact equality. Rebuilding in a new work dir reproduces a D/E/F artifact's content.
  Keeping the round's work dir is still the simplest practice, but refit noise is not a reason for it.
- The byte-level known answer for a dump must compare arrays, not the file sha.

## 452c (from 04:22Z 09-26): the `auto` weight end to end, with forced strong and weak banks
- Setup: `def_day.sh` (b517b6d2) on `ef-bundle-r` with DEF_DAY_FORCE D=K562, E=HCT116, F=A549 and w `auto` (builder
  4feae27d, driver 62b5110d).
- The same-line stage wrote w_D 2.0 (k562_gwps, strong), w_E 2.0 (orion_hct116, strong) and w_F 0.5 (jiang/myllia,
  weak), as designed. D and E cover 300/300 targets with their programmes, and `noprog` keeps them off the prior. F
  covers 6/300 with rows only.
- lvl0 ρ: D 0.8794 / E 0.7199, **identical to 452b** (same w). F reads 0.5918 against 452b's 0.5934: its six rows at w 0.5
  instead of 2 lower the spread slightly. First step: D 1.0553 / E 1.1664 (452b's) / F 1.2865 (452b 1.2847).
- lvl1: D 0.9755 / E 0.9680 (= 452b) / F 1.0565 (past its cliff; artifact deleted automatically).
- lvl2: 0.9755 / 0.9794 / 0.9658.
- **PASS (06:06Z): CONVERGED after 3 levelled builds** at D 1.0553 / E 1.1735 / F 1.2512 (452b: F 1.2498).
  - The final meta carries `sameline.w = "auto"` with w D 2.0 / E 2.0 / F 0.5.
  - Audits clean, prep dry run verified, builder restored.
- **The October default (`auto pt:noprog`) is now rehearsed end to end with mixed bank classes.** Its only difference
  from 452b is F's weight (0.5), and it moves F's factor by 0.6 %.
