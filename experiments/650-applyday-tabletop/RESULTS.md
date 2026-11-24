# 650 — Tabletop audit of the written 10-05 CHECKLIST and 10-06 APPLY DAY (docs/NEXT-SLOT-DECISION.md)

2026-10-02 16:05–16:35Z, session [7772dc], assigned by the coordinator. Audit base: HEAD `cce5d4b8`. The two procedure
sections are unchanged at `b613fbe2`: the only later edit appends 9 lines elsewhere.

## Method
- A scratch `git worktree` at `cce5d4b8`, outside the live tree. Every step was executed in order as written.
- Worktree-only setup (none of this is a doc issue):
  - gitignored inputs APFS-cloned in (`data/challenge|donors|harness|external`, gitignored `experiments/**/work` files, and
    the bundle h5ads identity554 reads);
  - `.venv` linked to the live venv, with `UV_PROJECT_ENVIRONMENT=<live .venv>`, `UV_NO_SYNC=1` and
    `PYTHONPATH=<worktree>/src`, so nothing was installed into the live venv.
- Heavy tests went through gate602.
- **The 10-05 commit (step 9) was simulated as a detached-HEAD commit inside the worktree, never pushed.** The 10-06 steps
  therefore ran on a faithful post-bundle HEAD. That is how 10-06 failures B2 and B3 were found.
- KA builds were not run. Their command lines, inputs, hash sources and grep targets were checked on paper.
- **The live tree was untouched:** the 15 audited paths are clean and unchanged since `cce5d4b8`, and no patch-installed
  file exists live.
- Runner: `step650.sh`. Per-step logs: `work/logs/`. Index: `work/650log.jsonl` (87 runs).

## Verdict
**The patches themselves are sound.** The full chain applied in the documented order reproduces 606's LAST PINNED block
(14:56Z, `85d4be68`) on all 12 files:

| file | sha |
|---|---|
| def_day | 54fcf040 |
| def_sameline | 5cd3ab07 |
| builder | fe4678c3 |
| def_identity_prior | 428f4a32 |
| rank_up | 7e6b3e52 |
| def_identity | 7073406f |
| helpers | fb715e38 |
| def_pipeline | 76314af6 |
| october_recipe | 9dbcd914 |
| def_preflight | 17f0836f |
| def_postbuild | 557c2b52 |
| gate524 | b16e60c0 |

The builder suite gives 56 passed, 0 failed, 1 skipped.

**The written procedure does NOT run as written.** Four blockers make test steps fail for ordering or tree-state reasons
rather than real defects. One KA rule cannot be met.

### Blockers
- **B1 (10-05 steps 2, 3, 6b): pre-application suites are run AFTER the apply.**
  - test_extras538, test_guard554 / identity554 and test_depth641 build temp trees from the LIVE files and re-apply 538
    (and 554) themselves.
  - Run in the doc's order, each dies with `already patched` and 0 tests run.
  - Run before any apply: test_guard554 5/5, test_depth641 4/4, identity554 as documented (below), test_extras538 12/13
    (see L3).
  - **Fix:** move every bundle suite into step 1, before step 2: test_identity590, test_extras538, test_guard554,
    identity554, test_depth641.
- **B2 (10-06 steps 3b, 4, 4a, 4b, 5, and the compose606 re-pin): suites that apply 641 themselves cannot run once 641 is
  on HEAD.**
  - On the post-10-05 HEAD, `compose606.py` reports `FAIL 641 --check rc=1: already patched` and SOME FAILED, even though
    the rest of the chain composes and the final shas equal the pins.
  - test630dd, test604 and test_skip644 die at setup with the same message. test604f (10-08+) shares compose606's
    DEFAULT_ORDER, so it will fail the same way.
  - test_hek293_642 dies at setup re-applying 538 to the live files.
  - So "every patch's suite green" (step 5) is unattainable on 10-06.
  - **Fix:** treat 641 as part of the bundle in `test598.BUNDLE`, which already skips a step whose `--check` says
    "already". Then drop "641" from `compose606.DEFAULT_ORDER` and `test630dd.PRE`, and have 642's tree builder start
    from `test598.canon_tree`. Alternatively, give every chain builder the same "skip if `--check` says already" rule.
- **B3 (10-06 step 1): test598 after the apply.**
  - In the doc's order (`patch598` apply → test598) the result is **42/43**.
  - T0 asserts canonical builder == LIVE builder, which no longer holds once 598 is applied.
  - Run between `--check` and apply, it is **43/43**.
  - **Fix:** "patch598 --check → test598 → apply → post598".
- **B4 (10-06 step 6): KA630 as specified cannot pass.**
  - The rule requires the def_day log line `630dd: uprn FAILSAFE_TOTAL <n>`.
  - A "run_ka_527.sh style" build goes through `def_pipeline.py` and never runs def_day.sh, so the line is never written.
    A missing line counts as KA FAIL.
  - **Fix:** build KA630 through def_day.sh, in the efp-relaxed form of the inert KA with `DEF_DAY_DE_MAG=dn_k250_zs
    DEF_DAY_DE_UPRN=<variant>`. Or compute FAILSAFE_TOTAL from the run_ka build's rc_audit dump with 630dd's own snippet.
  - Also spell out:
    - the variant: "K 1.0" is ambiguous between b_k100 and the shipping b_k100_q20r3p10 (13:43Z note);
    - the config-diff reference meta for "exactly `de_uprn`".

### Lesser findings, each with a one-line fix
- **L1. Bare script names:**
  - step 2 `apply_october538.py` / `test_extras538.py` → `experiments/538-cloud-bundle-extras/…`;
  - step 3 → `experiments/554-dose-guard/…`;
  - step 8 `run_arm_generic.sh` → `experiments/428-zero-gene-audit/run_arm_generic.sh`;
  - 10-06 3b `patch630dd` → `experiments/630-rank-neutral-up/patch630dd.py`.
- **L2. Missing `--root .`:**
  - 554 (`--root` is required: as written, argparse rc 2);
  - 10-06 598 / 612 / 612dd / 644 (required);
  - 590 and 604 default correctly.
- **L3. test_extras538 is 12/13 even before the apply, not 13/13.**
  - `test_recipe_pending_exports_nothing_new` expects a command without `DEF_DAY_SHARPEN_PT`, but sharpen_pt was DECIDED
    on 10-02.
  - Update its expected line to `DEF_DAY_LODO=1 DEF_DAY_SHARPEN_PT=tercile:2.5,2.0,1.75 DEF_DAY_DE_MAG=dn_k250_zs bash
    scripts/def_day.sh B W N 1.75 0.5 noprog 4`, or record 12/13 as expected.
- **L4. identity554 has no stated expected outcome. Write it down.** It prints `all identical: False`, and that is
  CORRECT per 554 RESULTS:
  - preflight identical on 6/6;
  - postbuild identical on 11/12 (465efr adds its dose-guard FAIL row at prior 0.278);
  - `--rc-audit` on 433a adds no row; 429a adds FAIL, then ALLOWED.
  It also **rewrites the tracked `experiments/554-dose-guard/work/identity554.json`**: add it to step 9's commit, or
  `git checkout` it.
- **L5. Step 5b is out of date.**
  - **4** tests fail, not 5: all_switches_on, riders_need_502, dependency_violations, empty_flags.
  - 538 adds **six** PENDING switches, not two: callset_block, hook_estimator, lodo_hi, lodo_strength, rc_callset_cpm,
    rho_scale.
  - Verified edit, giving 6/6 OK: add `("rho_scale", "off"), ("rc_callset_cpm", "off"), ("hook_estimator", "off"),
    ("lodo_strength", "off"), ("lodo_hi", "off"), ("callset_block", "off")` to those four tests' `self.decide([...])`.
  - The edited `scripts/test_october_recipe.py` is not in step 9's commit list.
- **L6. Step 7's `env -u` list (550's) misses five knobs the post-10-05 scripts read:** `DEF_DAY_ALLOW_CELLS_PER_PERT`,
  `DEF_DAY_ALLOW_DEPTH`, `DEF_DAY_ALLOW_LOW_DOSE`, `DEF_DAY_ALLOW_WIDE_LIBS`, `DEF_DAY_INCUMBENT_META`.
  - Add `-u` for each, or run the KA under `env -i HOME=$HOME PATH=$PATH DEF_DAY_CD4_THRESHOLDS=relaxed …`.
  - Paths, hash source and grep targets are fine: ef-bundle present; `ac49d319…` is in
    526/efp_inert/474efp_final.h5ad.meta.json; `DONE. Next` is a real log line.
- **L7. Step 10: after deciding lodo_hi_pct and mde_sign, `october_recipe.py command` still refuses.**
  - It reports `PENDING: rho_scale, rc_callset_cpm, hook_estimator, lodo_strength, lodo_hi, callset_block`; the doc names
    only the last two.
  - For the efr rehearsal, use `command --allow-pending` (which prints the pending ones as OFF), or decide the other four.
- **L8. Stale counts and shas:**
  - test612 is 23/23 (doc 22/22); test630 is 35/35 (doc 27/27);
  - patch630dd is `dddc5741` with module 7e6b3e52 (doc: patch 4db2a0c5 / 626949ef / module d8c39c3f);
  - the inline "As of f152e4c6" pins (def_day 7075d83d, def_sameline f147f8bb, def_identity_prior 7aeb66b4) are
    superseded. Delete them and point only at the LAST PINNED block.
- **L9. `test_rank_up630` needs arguments.** As written it raises IndexError. With `0 b_k100_q20r3p10 85` it gives PASS
  bit for bit (module 7e6b3e52), ~5 min, peak ~4 GB.
- **L10. 10-06 step 4a:** patch604 prints `WARNING def_sameline.py is ff898053…, not the post-bundle 22c6f0f0`. In this
  order it is benign: 642 ran first, and 604 / 642 commute per 606. Say so, so nobody stops there.
- **L11. 607's items have no commands.**
  - Item (2) measured here: H1 seed-0 dumps, 85 blocks, min(ballast library) − max(carrier library) ≥ 968,667. No block
    violates it, so 630's bug condition does not occur. Moot anyway after the 605 fix.
  - Item (3), counting the 612 targets tied at the key-gate cut: no procedure is given.

## Per-step record (exact command as run; ran Y/N; result vs expected)
| step | command (in the worktree) | ran | result | expected per doc |
|---|---|---|---|---|
| 0 | `bash experiments/550-cloud-1005-review/repro/procs_clear550.sh`; `git diff --quiet -- <11 paths>` | Y | BUSY rc 2 (recipe_dr524 pid 33446 waiting, expected today); 11 paths clean | CLEAR on 10-05 after 524dr's exit |
| 1 | `uv run python experiments/590-identity-validation/test_identity590.py` | Y | ALL PASS (8/8), min_top1 0.772 / 0.719 / 0.733 | pass ✓ |
| 2 | `…/apply_october538.py --root . --check`, then `--root .` | Y | clean; def_pipeline → 9a718b9e | ✓ |
| 2 | `…/test_extras538.py` (after apply, doc order) | Y | ERROR "already patched", 0 run | 13/13 ✗ (B1) |
| 2′ | same, BEFORE the apply | Y | 12/13 (L3) | 13/13 |
| 3 | `patch_guard554.py --check` (as written) | Y | rc 2: `--root` required | ✗ (L2) |
| 3 | `…/patch_guard554.py --root . --check`, then `--root .` | Y | def_postbuild 557c2b52, def_day 5851f84b | ✓ |
| 3 | `test_guard554.py`, `identity554.py` (after apply) | Y | "already patched" | 5/5 ✗ (B1) |
| 3′ | same, BEFORE any apply | Y | 5/5 OK; identity554 = 554's documented result (L4) | ✓ |
| 4 | `patch_identity590.py --check`, `--apply`; `scripts/def_identity.py --verify` | Y | c98a608f / 22c6f0f0; VERIFY PASS 0.772 / 0.719 / 0.733 | ✓ |
| 5 | `bash -n scripts/def_day.sh`; `october_recipe.py check --allow-pending` | Y | OK; rc 0 | ✓ |
| 5b | `scripts/test_october_recipe.py` | Y | 4 FAIL; 6/6 after the L5 edit | "5 fail" (L5) |
| 6 | `experiments/594-def-bundle-composition/post594.py` | Y | ALL PASS; every sha = 594's table except decisions (moves, as the doc says) | ✓ |
| 6b | `patch_depth641.py --root . --check`, then apply | Y | def_preflight 17f0836f, def_day 738cb0ae | ✓ |
| 6b | `test_depth641.py` (after apply) | Y | "already patched" | 4/4 ✗ (B1); 4/4 BEFORE any apply |
| 7 | inert KA via def_day (550 step 4) | N (build) | paths / hash source / grep target OK | `env -u` list incomplete (L6) |
| 8 | driver-list edit (`9a718b9e…`); `bash -n`; `bash experiments/591-prep-gate/test_gate591.sh` | Y | ALL PASS | ✓ |
| 9 | commit the 11 files + run_arm_generic.sh | simulated (detached HEAD in the worktree) | — | add test_october_recipe.py + identity554.json (L4, L5) |
| 10 | `october_recipe.py decide …` (on the worktree copy); `command` | Y | `command` refuses: 6 PENDING | doc names 2 (L7) |
| 11 | chain1006 | checked live | pid 40686 armed since 09-29 for 10-05 20:00Z | ✓ |
| 10-06 1 | `patch598.py --root . --check` → test598 → apply → test598 → `post598.py` | Y | 43/43 before apply; **42/43 after** (T0); post598 ALL PASS, def_pipeline 76314af6 | B3 |
| 2 / 2b | `patch612.py --root .`, `patch612dd.py --root .`; `test612.py` | Y | 23/23 (before and after) | 22/22 (L8) |
| 3 | `patch_builder630.py --root .`; `test630.py`; `test_rank_up630.py` | Y | 35/35; module 7e6b3e52; rank_up KA PASS with args | 27/27 (L8); args (L9) |
| 3b | `patch630dd.py --root .`; `test630dd.py` | Y | applies; test630dd setup FAIL (641) | 39/39 ✗ (B2) |
| 4 | `patch_hek293_642.py --root .`; `test_hek293_642.py` | Y | applies; test setup FAIL (538 re-apply) | ✗ (B2) |
| 4a | `patch604.py`; `test604.py`; `test604_guard.py` | Y | applies with WARNING (L10); test604 setup FAIL (641); guard 5/5 | 27/27 ✗ (B2) |
| 4b | `patch_skip644.py --root .`; `test_skip644.py` | Y | applies; setup FAIL (641) | 4/4 ✗ (B2) |
| pins | end-of-chain shas vs 606's LAST PINNED | Y | **12/12 equal** | ✓ |
| — | `compose606.py` on the post-10-05 HEAD | Y | SOME FAILED (641 already on HEAD); final shas = pins | ✗ (B2) |
| 5 | `experiments/410-pertscale-revisit/run_builder_tests.py` | Y | 56 passed, 0 failed, 1 skipped | ✓ |
| 6 | KA builds | N (paper) | hashes found: ec3a5ae6 (446bj2z meta), dd82740d (479p2g), e21ff752 (496a); ka598 / ka612 present; FAILSAFE line written only by def_day | KA630 unattainable as specified (B4) |
| 7 | driver list += 76314af6; `bash -n`; test_gate591 | Y | ALL PASS | ✓ |
| 9 | `october_recipe.py decide sharpen_pt tercilec:2.5,2.0,1.75` | Y | DECIDED; `check` OK with tercilec | ✓ |
