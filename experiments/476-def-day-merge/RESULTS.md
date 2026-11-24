# 476 — 485's October day-of fixes merged onto 474 (plan and merge, 2026-09-27 05:35Z)

**Why.** The cloud review 485 found eight October day-of defects in `scripts/def_day.sh` and its helpers (F1 critical: an
in-place retry ships the previous h5ad under the new factor string; F2: a bad secant step aborts before the safe fallback;
F3–F8 medium/low). Its patches were cut against the pre-474 files; 474 (481's patches) has since changed `def_day.sh`,
so 2 of its 9 hunks conflict. The other three patches (helpers, driver dump-assert skip, relevel_iter) apply clean, and
current+patch equals 485's `repro/*.new` byte for byte.

**Merge** (`work/def_day.resolved.sh`, from `git merge-file` ours = 474's 91f350f0, base = 847cced5^, theirs = 485's .new):
- both sides fixed F8 independently. Kept 474's rehearsed forms (contexts from the bundle manifest before setup, `UNIT`,
  `--cd4-thresholds "$CD4T"` on the setup call) and added 485's checks on top: `DEF_DAY_CD4_THRESHOLDS` must be `relaxed`
  or unset; the preflight's contexts must be three bare names and must equal the manifest's;
- conflict 1 (build): 485's body (delete the tag's h5ad and meta first, require driver rc 0, past_cliff rc 2 = stop),
  comparing with 474's `$UNIT`;
- conflict 2: 474's rule-10 postbuild (5b) first, then 485's strict zero-gene audit (6);
- everything else merged without conflict: dump rm + rc 0 + `[bulk-dump] wrote` in the stage log; relevel rc≠0, no
  factor line, or out-of-bound proposal → best_safe; lvlC only when `inband`; meta_matches before the final link;
  nctx must print a count; programme source with 0 matched contexts stops; prep dry-run failure stops; restore trap
  fails closed (KEEP_PINS option).
- `repro/`: 485's `repro_fixes.py` run against the merged script (the harness copy defines `UNIT` next to `FACT`, the
  only harness change): **all fix checks matched the claims** (rc 0).

**Known gaps (for the October switches still pending):** `meta_matches` parses q as a float and compares
`bulk_sharpen_q`, and 5b's allowed-diff list has no `bulk_sharpen_pt`. A per-target sharpening pass-through (479's S_P2
for October) must extend both before it is used.

**Rehearsals** (`chain476.sh`, serial after chain495):
| check | known answer |
|---|---|
| apply | each applied file == repro/*.new; repro_fixes PASS; run_arm_generic's driver pin list + the new driver sha |
| (a) efp strict | stops at setup with the CD4-gate message |
| (b) 465 on efr (LODO, programme h1_2025_full, pt) | final output == 465's `6113937f`; postbuild PASS; the dump stage exits 0 |
| (c) efp relaxed | final output == 474 (c)'s `ac49d319` (D 0.9831 E 1.0324 F 1.0663); postbuild PASS |
| (d) in-place rerun of (c), q 1.5, max_levels 1, (c)'s h5ads on disk | STOP loudly, or DONE with a fresh build (meta q 1.5, output ≠ ac49d319, ρ ≤ 1) |

Any failure restores the five pre-476 files.

## Outcome (2026-09-27, chain476 11:26:56Z → 18:07:00Z): all four rehearsals PASS; committed 5e761657 by its committer
- **apply:** def_day `e3247feb`, helpers `346495d7`, driver `9b2c7eca` (485's dump-assert skip), relevel `bd10f8f2`, each
  equal to its repro/*.new; repro_fixes PASS; run_arm_generic's pin list gained the new driver (→ `c58b9553`).
- **(a)** efp strict stops at setup with the CD4-gate message (186 panel targets map; the gate wants 290). 11:27Z.
- **(b)** 465 retraced on efr byte for byte (`6113937f`); converged at D 1.2198 E 1.2842 F 1.2496 after 3 levelled
  builds; postbuild PASS; the dump stage now exits 0 (the driver logs the arm-assert skip). 14:31Z.
- **(c)** 474's efp relaxed retraced byte for byte (`ac49d319`); converged after 1 levelled build; postbuild PASS. 17:17Z.
- **(d)** in-place rerun of (c) at q 1.5, max_levels 1, with (c)'s two h5ads still on disk: **fail-closed, and both
  guards worked.** It rebuilt fresh instead of reusing (c)'s artifact (lvl0 at q 1.5 had D ρ 1.0721, past the cliff, and
  was deleted; lvl1 was a new build at D 0.9557 E 1.0000 F 1.0354, ρ 0.9820 / 0.9715 / 0.9775, not converged in 1
  level, shipped as the nearest safe build). The postbuild gate then stopped the day because `bulk_sharpen_q` (1.5) is not
  an allowed difference from 446b: the only unexpected key. The pre-485 script would have shipped (c)'s h5ad under q 1.5's
  factor string.
- The builder was HEAD after every rehearsal. chain497 started at 18:08Z on top of these files (def_day `e29e793a`).
