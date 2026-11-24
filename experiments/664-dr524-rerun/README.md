# 664: re-run of 524's October-recipe rehearsal after a memory refusal (coordinator, 2026-10-05)

## What happened
- `recipe_dr524.sh` released at 02:00:41Z. Its gate passed: free+inactive ≥ 25 GiB, FLOOR, no build, HEAD-equal scripts.
- It resolved the recipe with `mde_sign` still PENDING, so the hook was OFF, as 658 assumes. Then it started `def_day.sh`.
- 31 s later, def_pipeline's own pre-build check refused: "free+inactive memory 22.5 GiB < 23.0 GiB; refusing"
  (`work/run.aborted-*/driver-setup.log`). def_day.log: "STOP: setup did not complete".
- The pinned builder was restored (4daae237 == HEAD, "builder after restore").
- Cause: two iOS Simulators (GlowUp, both booted 10-04) keep free+inactive at 19–26 GiB. Setup's identity and key stages
  (2.7 / 4.1 GiB peaks) ate the 2 GiB margin between the 25 GiB release gate and def_pipeline's 23 GiB floor.
- A full def_day holds ~20 GiB for ~24 min (497 efr: repin-panel-support 20.1, build 20.2, prep 22.2 GiB peak). The refusal
  was correct.

## The fix: `rerun_dr524.sh` (armed 02:05Z, pid 18948)
- It re-runs `recipe_dr524.sh` UNCHANGED. The script moves the failed run dir aside itself and re-resolves the recipe, with
  `mde_sign` still OFF.
- It releases only at free+inactive ≥ 27 GiB on two readings 60 s apart, with no build or prep running.
- It retries only on the same memory refusal, after 5 minutes. Any other failure stops it for diagnosis.
- Cap 10-05 19:00Z. chain1006 tracks the rehearsal through `recipe_dr524.pid` and waits for its process until 23:30Z.
- If 524dr is not READY at 20:00Z, chain1006 arms 10-06 s1 = K16, and s2 is empty: K20 was scored 10-04, so the novelty gate
  skips it.
- Ahmet was asked by push to quit the Simulator.

## Ordering kept
- The coordinator decides `mde_sign = clean4:0.05` (496a TRANSFERS) only after this rehearsal exits.
- The 10-05 CHECKLIST v2 also waits for it (step 0).
