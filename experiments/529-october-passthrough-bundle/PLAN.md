# 529 — the October pass-through bundle (plan, 2026-09-29 11:05Z)

One `def_day.sh` / `def_day_helpers.py` / `october_recipe.py` patch that adds every remaining October knob, each default
OFF and byte-identical when unset, proved by ONE inert known answer and ONE rider-on rehearsal. Planned in
NEXT-SLOT-DECISION (09-28 19:45Z) and FINAL-PHASE-FEEDBACK-PLAN §4; this file fixes its contents and schedule.

## Knobs (each: a `DEF_DAY_*` export, a resolver switch in `docs/october-decisions.json`, a `riders_match` row)
| knob | builds | why | decided by |
|---|---|---|---|
| `DEF_DAY_RHO_SCALE=<s>` | ρ* × s in the re-level loop (relevel_iter's target) | FINAL-PHASE rows 2: ρ* is the most D/E/F-specific number; needed to run the ×0.98 / ×1.02 arms if members come back | always installed, default 1 |
| `DEF_DAY_CD4_STACK=x1\|x1.5` | the stack dump's CD4 multiplier | FINAL-PHASE row 4; 457dr read +0.0031 on A/B/C (k 1.43), 438 says ×1.5 loses on October-like panels | installed, default x1 |
| `DEF_DAY_RC_CALLSET_CPM=<n>` | `VCC_RC_CALLSET_CPM` on the level builds | 523's K pair (10-06/07) | the K read |
| `DEF_DAY_MDE_SIGN=hookbank:…` + `DEF_DAY_HOOK_ESTIMATOR=meanz_ne` | 517's `patch_defday517.py` + `patch_sameline517.py` (Arc-sourced hook on a CLEAN4 block) | 527's arm 446bj2zmz (10-01 null branch or 10-07) | the 446bj2zmz read (+ 496a's) |
| `DEF_DAY_LODO_STRENGTH=<θ>` | 522's pct × exp(0.25 z_src) rider in `def_lodo_key.py` | +0.7 … +1.9 % direction on all 9 October-like folds, ≈ 0 on A/B/C (unreadable on the board) | only if 495 (cap 6) transfers on 10-03 |
| `DEF_DAY_CALLSET_BLOCK=prop` | 528's call-propensity block, computed from the day's own controls | 528 (worker, running) | 528's H1 price, then a board read |

## Schedule and safety
- `def_day.sh` must not change while `recipe_dr524.sh` runs it (≥ 10-05 02:00Z, ~1.5 h). The bundle is applied by a chain
  that waits for recipe_dr524's PROCESS to exit, then: apply → `bash -n` → the helpers' unit tests (497/502/517 suites) →
  `repro_fixes` → inert KA (efp relaxed with no new knob == `ac49d319`, as 497 and 502 did) → one rider-on rehearsal on
  relabelled A/B/C with every knob that has been DECIDED by then → `october_recipe.py check` → committer (spec file).
- Before writing it: fold in 526's review findings (cloud, running) so the bundle does not stack on a known defect.
- Each knob stays OFF in `docs/october-decisions.json` until its read decides it (PENDING → DECIDED via
  `october_recipe.py decide`); the resolver maps only DECIDED-on switches to exports.

## Status 2026-09-29 14:15Z
- Cloud 530 delivered the bundle (rho scale, CD4 stack, call-set CPM, Arc hook + meanz_ne, strength rider): tested locally,
  16/16, `--check` clean on today's tree. Not applied; applied after recipe_dr524's 10-05 run (see 530's RESULTS).
- Still to add in the same application (from 531's tabletop and 534's Part C): a `DEF_DAY_LODO_HI` pass-through (fixed cap,
  refused with `--hi-pct`); a documented census override for R0 (531 F16); `def_postbuild.py --require-config-diff` (531
  F17); a default-off knob for any base step 534 Part C marks FAILS: **resolved 09-29, no knob** (the cap form and the partial rows FAIL
  within ±0.0003; both stay on, min-max regret; normrel HOLDS). 528's call-propensity block (`DEF_DAY_CALLSET_BLOCK`) only if its H1 price and a board read
  are positive.

## Status 2026-09-30 11:44Z
- Cloud 538 added the remaining pieces (fixed LODO cap, the call-propensity block re-anchored after 530, census override,
  `--require-config-diff`, REPO from the script path) and one applier that runs 530 + 528 + 538 in order: merged, tested
  locally 13/13, `--check` clean. Apply on 10-05 after recipe_dr524's process exits (538 RESULTS "Commands …").
