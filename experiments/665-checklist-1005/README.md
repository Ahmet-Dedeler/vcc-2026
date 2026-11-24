# 665: the 10-05 CHECKLIST v2 run (coordinator, 2026-10-05)

The procedure is `docs/NEXT-SLOT-DECISION.md` "10-05 CHECKLIST v2". It started after 524dr's re-run exited (07:05Z, gate524 rc 0).

| step | result |
|---|---|
| 0 | procs_clear550 CLEAR; no run_arm_generic.sh running; live == HEAD on the 11 paths |
| 1 | suites: 590 ALL PASS; 538 13/13; 554 5/5; identity554 "all identical: False" (correct; json restored from a byte copy); 641 4/4 |
| 2 | 538 `--check` clean, applied; def_day.sh cb82c34f, def_pipeline 9a718b9e; second apply "already patched" |
| 3 | 554 `--check` clean, applied; def_postbuild 557c2b52, def_day 5851f84b, def_preflight 2d181194 |
| 4 | 590 `--check`, `--apply`; def_identity c98a608f, def_sameline 22c6f0f0; `def_identity.py --verify` VERIFY PASS 0.772 / 0.719 / 0.733 |
| 5 | `bash -n scripts/def_day.sh` OK; `october_recipe.py check --allow-pending` rc 0 |
| 5b | test_october_recipe.py: 4 FAIL after 538, as expected; the six 538 switches decided off in those tests (`OFF538`) → 6/6 OK |
| 6 | post594.py ALL PASS |
| 6b | 641 `--check` clean, applied; def_preflight 17f0836f, def_day 738cb0ae |
| 7 | inert KA (`run_ka1005_inert.sh`) **PASS**: both simulators were shut down by their sessions at 09:08Z, the runner was re-armed at a 25.5 GiB floor, attempt 1 ran 09:39–10:35Z (memory checks 30.3–32.8 GiB); output.sha256 `ac49d319…` == 526/efp_inert; postbuild dose_guard per context; DONE. Next; builder restored == HEAD |
| 8 | run_arm_generic.sh driver list += 9a718b9e; `bash -n` OK; test_gate591 ALL PASS |
| 9 | commit (after 7): the 11 files + run_arm_generic.sh + scripts/test_october_recipe.py, with the KA hash in the message |
| 10 | `october_recipe.py decide mde_sign clean4:0.05` (496a TRANSFERS) and hook_estimator's non-Arc value (652 F1) |

Disk: idall3's 12 GB build h5ad was swept before step 7 (scored 4R0FMEtcLtvZMmg2X0G5; its .vcc sha matched the logged
dad0e10e), taking free space from 33 to 44 GiB.
