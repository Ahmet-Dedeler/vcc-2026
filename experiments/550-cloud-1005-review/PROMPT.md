You are an adversarial reviewer for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
first (§6.2 waiter and process rules, §10 storage, §11 method rules, §14), then docs/DEF-RUNBOOK.md and the bottom of
docs/NEXT-SLOT-DECISION.md (sections dated 2026-09-30). This VM has NO gitignored data (no data/challenge matrices, donors,
harness h5ads, npz artifacts): review code statically and test with synthetic inputs and stubs only.

Experiment 550. Write ONLY under experiments/550-cloud-1005-review/. Commit to your branch only. Never run `vcc`, never contact
anyone, never touch GitHub issues or PRs, never edit files outside your directory. You deliver findings, repro scripts and (if
needed) PATCH SCRIPTS that are not applied.

## Context
On 2026-10-05 the coordinator changes the October (D/E/F) path, the code that builds the prize entry after Oct 22 on a dark
leaderboard. Anything wrong here can zero the final or silently change the recipe. The plan for that day (NEXT-SLOT-DECISION
2026-09-30 11:40Z and 11:44Z items):
1. 02:00-~04:00Z `experiments/524-october-recipe-dr/recipe_dr524.sh` runs automatically (a dress rehearsal through
   `scripts/def_day.sh`, then `gate524.py --write`).
2. After its PROCESS exits: apply the October bundle with ONE applier, `experiments/538-cloud-bundle-extras/apply_october538.py`
   (it runs 530's `patch_bundle530.py`, then 528's `patch_defday528.py`, then 538's extras), using the command list in
   `experiments/538-cloud-bundle-extras/RESULTS.md` ("Commands the coordinator runs on 10-05").
3. The inert known answer: efp relaxed through the patched `def_day.sh` with every new knob unset must reproduce
   `ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3` byte for byte (how it was run before:
   `experiments/502-october-de-riders/chain502.sh`, `experiments/526-cloud-october-review2/fix526.sh`,
   `experiments/476-def-day-merge/chain476.sh`). Then a committer via `scripts/commit_after_chain.sh <spec>`.
4. One rider-on rehearsal with the switches `scripts/october_recipe.py decide` has turned on. The likely October DE form is the
   call-propensity block + the CLEAN4 uncalled-gene hook (the A/B/C arm for it is 446bj2zcph, built by
   `experiments/549-cph-arm/build549.sh`, specified in 537's RESULTS §5.2 item 1 and §3.2).

Snapshots taken for you (2026-09-30 21:16Z; 537 RESULTS as last written at 20:31Z), because the live files are uncommitted or in other agents' dirs:
`experiments/550-cloud-1005-review/snapshots/537-RESULTS-v4-20260930T2031Z.md` (worker 537's current RESULTS; master holds an
older v1 at `experiments/537-arm-pipeline/RESULTS.md`), `snapshots/446bj2zcp.h5ad.meta.json` (the meta of the 10-01 slot-2 arm
built with the call-propensity block) and `snapshots/run446bj2zcp.log` (its `run_arm_generic.sh` log). Cite 537 by that
snapshot's path.

## Review these, and try hard to break each one
A. **Drift since 538 was checked.** 538 was verified on 2026-09-30 11:40Z. Find every commit since then that touched any of the
   nine files the applier patches (list in 538's RESULTS "Applier"), or anything they import or call. On today's master: run
   `apply_october538.py --root <temp copy> --check` and `test_extras538.py`; report pass/fail with output. Any anchor drift is
   a finding.
B. **The 10-05 order of operations.** (1) Which live bash processes could still be executing any of the nine files when the
   applier runs (AGENTS §6.2: bash reads scripts incrementally; `recipe_dr524.sh` runs `def_day.sh` as a child; chains
   `experiments/524-october-recipe-dr/chain1006.sh` at 10-05 20:00Z, `464`'s and `496`'s `arm100x.sh`/`chain100x.sh`, and
   anything else armed). Say exactly which process must have exited, and how the coordinator should test that (process-name
   traps: match interpreter + path, never `pgrep -f` on a bare name). (2) Does anything that runs AFTER the patch (chain1006,
   rule scripts, gate524, readers) depend on the unpatched behaviour, a file sha, or a log line whose text changes? (3) Does
   the inert-KA procedure used before (chain502 / fix526) still work unchanged on the patched tree (pinned shas, greps of log
   text, arguments)? Write the exact KA command sequence for 10-05.
C. **The instructions disagree.** `537-arm-pipeline/RESULTS.md` §5.3 says "apply 530's bundle; re-test `patch_defday528.py` on
   the bundled def_day.sh", while 538 says its applier INCLUDES 530, "do not also run patch_bundle530.py", and 528's patch
   no longer applies after 530 (538 re-anchors it). Establish which is correct by running both sequences on temp copies, and
   write the single correct 10-05 command sequence (check, apply, tests, `bash -n`, `october_recipe.py check`, KA, committer,
   rider-on rehearsal).
D. **Can the October path express exactly the A/B/C arm that will be read?** 549 builds with
   `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ec…` and
   `VCC_DE_MDE_SIGN=clean4:0.05`. On the patched `def_day.sh`, which `DEF_DAY_*` settings produce the same builder env
   (`DEF_DAY_CALLSET_BLOCK=prop` → swap/steer defaults; 502's `DEF_DAY_MDE_SIGN`)? Do the defaults (swap 1.0, steer `kept`)
   equal what built 04fb24ec (read 528's and 535's scripts and `snapshots/446bj2zcp.h5ad.meta.json`)? Does the day-of computation (`callblock528.py` as called by the
   patched def_day.sh) use the same function and parameters as the script that built the A/B/C artifact? Prove equivalence
   on synthetic controls (same inputs → same npz bytes or same arrays), or show the difference.
E. **549's build path vs 537 §5.2 item 1.** Env set/unset lists, REF_META, EXPECT (and whether
   `experiments/428-zero-gene-audit/run_arm_generic.sh` compares EXPECT as a set or as a string, and what happens on order),
   the artifact sha check, the builder/emitter pins, the `busy()` pattern in build549.sh (can it match itself, its children, a
   headless worker's argv, or miss a python builder?), the SIGSTOP/SIGCONT handling of 537's driver (what if the build dies,
   is killed, or the Mac reboots mid-build: is the driver left stopped?). 537 says "pds and mse must equal 446bj2z's to
   ≤ 1e-6 raw in every DE arm (falsifier: a bulk leak)". Is there any LOCAL check before a slot that a DE-only arm left the
   bulk identical to 446bj2z's (e.g. a pseudobulk or census digest in the meta, `postbuild.json`, the CLIFF check)? If none,
   write a cheap check script (reads two prediction h5ads or their `.vcc`, compares per-(context, target) summed counts,
   streams with backed mode, never materialises the matrix; test it on synthetic h5ads) and say where it should run.
F. Anything else on this path that could fail on 10-05 or on Oct 22, or silently change the recipe. Prefer one real defect
   with a repro over ten style notes.

## Deliver
`experiments/550-cloud-1005-review/RESULTS.md`: findings ranked by severity (BLOCKER / HIGH / MEDIUM / LOW), each with
file:line evidence and a repro under `experiments/550-cloud-1005-review/repro/` (synthetic data, stubs for `uv`/python as 530
and 538's tests do), what you ran and its output, what you could not test without data, and the single corrected 10-05
command sequence. Keep it tight. Commit to branch `cloud/550-1005-review`.
