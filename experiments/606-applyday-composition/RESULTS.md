# 606 — every pending patch composes, in the coordinator's apply order (council 600; [4d0207], 2026-10-02 ~12:40Z)

Several UNAPPLIED patches edit the same files:
- def_day.sh: 641, 598, 612dd, 642;
- def_sameline.py: 642, 604;
- the builder: 598, 612, 630.

Each patch's own test composes it only with its predecessors. `compose606.py` applies ALL of them on a scratch copy of the
canonical post-bundle tree (test598.canon_tree: HEAD + 538 → 554 → 590). It runs each patch's `--check` first, then
`bash -n` on def_day.sh and py_compile on every touched Python file. The live repo is hashed before and after and was
not touched.

| order | result |
|---|---|
| 641 → 598 → 612 → 612dd → 630 → 642 → 604 (the 10-05 / 10-06 plan, + 612dd and 604) | **ALL PASS**: every --check and apply clean; `bash -n` OK; py_compile OK |
| 641 → 598 → 612 → 612dd → 630 → 604 → 642 (604 before 642) | **ALL PASS**; the final def_sameline.py is the SAME (`3fd2dc18…`), so 604 and 642 commute |

Final shas of the composed tree (for the coordinator's post-apply checks):
- def_day.sh `736c71ae90c9…`
- def_sameline.py `3fd2dc18926f…`
- def_identity.py `7073406f6aad…`
- def_day_helpers.py `fb715e3833f3…`
- def_pipeline.py `76314af6b249…` (equals 598's predicted post-598 value)
- october_recipe.py `9dbcd9149020…`
- builder `33dc7a47e1b2…`
- def_preflight.py `17f0836fd1e1…`
- def_postbuild.py `557c2b522b56…`
- gate524.py `b16e60c08e48…`
- def_identity_prior.py `fe3fbeed3c8b…`

Re-run before the apply day (and after any patch changes) with
`uv run python experiments/606-applyday-composition/compose606.py [--order ...]`. This checks anchors and syntax only. Each
patch's own tests and the KA builds remain the real gates.

## Re-run 2026-10-02 ~12:55Z: 644's skip patch v2 as the last link (coordinator request)
`patch_skip644.py` v2 (sha `4c0bb12e…`, replaces v1 `11716c46…`) now edits def_sameline.py too (`--skip-contexts`) and
anchors on 642's lines. compose606 gained `"644"` as the 8th link (log `work/compose606_644.log`).

| order | result |
|---|---|
| 641 → 598 → 612 → 612dd → 630 → 642 → 604 → **644** (the 10-06 order in the live runbook) | **ALL PASS**: every --check and apply clean; `bash -n` OK; py_compile OK; live repo untouched |
| 641 → 598 → 612 → 612dd → 630 → 642 → **644 → 604** | **ALL PASS**, the same final shas: 644 and 604 commute (`work/compose606_644_before604.log`) |

Final shas of the composed tree (these replace the table above for the post-apply checks):
- def_day.sh `3918e0537ee5…` (was 736c71ae: +644), matching [537757]'s expected value
- def_sameline.py `d964aa03f6e2…` (was 3fd2dc18: +644), matching [537757]'s expected value
- builder `16f1fda4fb4d…` (was 33dc7a47). This is NOT from 644, which touches no builder file. 630 v3 (ae39c950,
  12:38Z, after the first 606 run) changed patch_builder630's edits (the v3 knob grammar). The builder is touched only
  by 598, 612 and 630.
- src/vcc2026/rank_up.py `8e1f05a66f4c…` (630's installed module, = `src_rank_up.py`, test630's module sha)
- unchanged: def_identity.py `7073406f…`, def_day_helpers.py `fb715e38…`, def_pipeline.py `76314af6…`,
  october_recipe.py `9dbcd914…`, def_preflight.py `17f0836f…`, def_postbuild.py `557c2b52…`, gate524.py `b16e60c0…`,
  def_identity_prior.py `fe3fbeed…`

644 × 604 semantics, read from the code. 644 marks a skipped context's record `mode = none (644 …)` and `continue`s
before any bank, 642 or 590 logic. It never rewrites `call`. 604's identity branch therefore still passes the skipped
context's TRUE identity call to `identity_prior`, so the same-call guard (never use a reference with the context's own
call) still holds. The skipped context is not `matched`, so it gets its identity programme at w 0 with no rows and no
signs, as 644 documents.

## Re-run 2026-10-02 ~13:35Z: 604 v3 (identity prior: no identity-only extension, side = sign(score))
`work/compose606_604v3.log`: **ALL PASS** in the 10-06 order (641 → 598 → 612 → 612dd → 630 → 642 → 604 → 644). The only
sha that moves is the module patch604 installs: def_identity_prior.py `a2d48da00f01…` (was `fe3fbeed…`). def_sameline.py
stays `d964aa03…`, because patch604's def_sameline edit did not change. All other shas are as in the 12:55Z table.

## Re-run 2026-10-02 ~13:25Z: 612's tie-at-cut fix (740f5bcf; patch612 `f5ecaf2b…`, 612dd unchanged `3fb0e5b1…`)
`work/compose606_612fix.log`: **ALL PASS**, builder `c2ee0e618876…`. This is NOT a pin: 630's `src_rank_up.py` was being
edited (d4cd5f's 605 donor-gain fix, uncommitted) when it ran, and the run picked up that work-in-progress module
(rank_up `d8c39c3f…`). The pins come from the re-run after 630's fix is committed.

## PINNED 2026-10-02 ~13:35Z: every patch committed (630's 605 fix 69972aa1, 612's tie fix 740f5bcf, 604 v3, 644 v2)
`work/compose606_630fix.log`: **ALL PASS** in the 10-06 order 641 → 598 → 612 → 612dd → 630 → 642 → 604 → 644. The working
tree was clean for every patch input when it ran. Post-apply shas for the coordinator's checks:

| file | sha (12 hex) |
|---|---|
| scripts/def_day.sh | `3918e0537ee5` |
| scripts/def_sameline.py | `d964aa03f6e2` |
| scripts/def_identity.py | `7073406f6aad` |
| scripts/def_day_helpers.py | `fb715e3833f3` |
| scripts/def_pipeline.py | `76314af6b249` |
| scripts/october_recipe.py | `9dbcd9149020` |
| scripts/make_submission_pm_cd4.py (builder) | `c2ee0e618876` |
| scripts/def_preflight.py | `17f0836fd1e1` |
| scripts/def_postbuild.py | `557c2b522b56` |
| experiments/524-october-recipe-dr/gate524.py | `b16e60c08e48` |
| scripts/def_identity_prior.py (604 v3) | `a2d48da00f01` |
| src/vcc2026/rank_up.py (630, fixed) | `d8c39c3f4e1a` |

Patch inputs at pin time:
- patch612 `f5ecaf2b…`, patch612dd `3fb0e5b1…`;
- patch_builder630 `d2a2034f…`, which reads the module sha at patch time, so only its DE_UPRN_MODULE_SHA256 line
  follows the module;
- 630's module `d8c39c3f…`;
- patch_skip644 `4c0bb12e…`.

The builder moved from 16f1fda4 to c2ee0e61 through two edits: 612's tie fix and 630's module-sha line.

## 2026-10-02 ~13:40Z: + 630dd (DEF_DAY_DE_UPRN → VCC_DE_UPRN; [4d0207] for [d4cd5f]) — the chain is now 9 links
`experiments/630-rank-neutral-up/patch630dd.py` (sha in the commit) + `test630dd.py` **39/39 PASS** (`630/work/test630dd.log`).
It mirrors 612dd's six anchors (DOC / V / E / C / F / A) on the ORIGINAL DE_MAG lines. V requires DEF_DAY_DE_MAG, the
builder grammar with K in [0.5, 3], and the builder's "630: rank-neutral UP" mark. C requires the final meta's
`de_uprn.variant` == the request AND `module_sha256` == d8c39c3f (the reviewed module; patch time refuses a builder
pinned to any other module). F and A record the knob. Unset = no export, no check, no final.env line, no ALLOW change
(T4, run on the inserted text in a bash harness).

| order | result |
|---|---|
| 641 → 598 → 612 → 612dd → 630 → **630dd** → 642 → 604 → 644 (`work/compose606_630dd.log`) | **ALL PASS** |
| 641 → 598 → 612 → 630 → **630dd → 612dd** → 642 → 604 → 644 (`work/compose606_630dd_first.log`) | **ALL PASS, same shas**: 612dd and 630dd commute (test630dd T2: identical def_day.sh bytes) |

The only post-apply sha that moves from the 13:35Z pin is **def_day.sh `7f78bd047cb1…`** (was 3918e053). All others are
unchanged: builder c2ee0e61, rank_up d8c39c3f, def_sameline d964aa03, def_identity_prior a2d48da0, etc.

## RE-PINNED 2026-10-02 ~13:55Z: 630's module → 7e6b3e52 (c30d64ce: invariant asserts fail safe) and 630dd v2
630dd v2 makes two changes:
- it pins module 7e6b3e52;
- it adds a NON-FATAL per-context report from the final rc_audit dump: mean frac_got / genes moved, guard_failed and
  invariant_failed target counts, the prereg's capacity floor. It is logged and written to `<work>/uprn_audit.json`.

test630dd **43/43 PASS**. `work/compose606_mod7e6b.log` gives **ALL PASS**, and `work/compose606_mod7e6b_ddswap.log`
(630dd before 612dd) gives identical shas.

| file | sha (12 hex) | change vs 13:35Z |
|---|---|---|
| scripts/def_day.sh | `f9552b6bc8d3` | was 3918e053 (+ 630dd) |
| scripts/make_submission_pm_cd4.py (builder) | `2c6cfe7d6ad2` | was c2ee0e61 (630's module-sha line) |
| src/vcc2026/rank_up.py | `7e6b3e52e2ad` | was d8c39c3f |
| everything else | unchanged | def_sameline d964aa03, def_identity 7073406f, helpers fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914, def_preflight 17f0836f, def_postbuild 557c2b52, gate524 b16e60c0, def_identity_prior a2d48da0 |

## RE-PINNED 2026-10-02 ~14:05Z: 630dd v3 (greppable `630dd: uprn FAILSAFE_TOTAL <n>`; coordinator: n > 0 on a KA630 build = KA FAIL)
test630dd **45/45 PASS**; `work/compose606_failsafe.log` **ALL PASS**. **def_day.sh `c2704f3d2e36`** (was f9552b6b). Every other
sha is as at 13:55Z: builder 2c6cfe7d, rank_up 7e6b3e52, def_sameline d964aa03, def_identity_prior a2d48da0, … **Chain
freeze: 10-05 18:00Z** (coordinator). After it, patch changes wait for 10-08, and the 10-06 apply uses the pin current at
the freeze.

## RE-PINNED 2026-10-02 ~14:05Z: 607 follow-ups. 604's fixes (7aeb66b4 / def_sameline f147f8bb), 630's knob grammar (q / r / p ≥ 1, no leading zero; patch_builder630 `5173d789…`), and 630dd mirroring it (test630dd 50/50)
`work/compose606_grammar.log` **ALL PASS**:

| file | sha (12 hex) |
|---|---|
| scripts/def_day.sh | `7075d83de882` |
| scripts/def_sameline.py | `f147f8bbe114` |
| scripts/make_submission_pm_cd4.py (builder) | `fe4678c31b3c` |
| src/vcc2026/rank_up.py | `7e6b3e52e2ad` |
| scripts/def_identity_prior.py | `7aeb66b44f28` |
| unchanged | def_identity 7073406f, helpers fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914, def_preflight 17f0836f, def_postbuild 557c2b52, gate524 b16e60c0 |

## PINNED 2026-10-02 ~14:20Z (authoritative until superseded): all inputs committed, including 644 v3 (d4a3d275), 630dd v4 (bd76246b: FAILSAFE_TOTAL -1 on a missing / bad dump), 630's grammar (52356251) and 604's 607 fixes
`work/compose606_644v3.log`: **ALL PASS** in 641 → 598 → 612 → 612dd → 630 → 630dd → 642 → 604 → 644. Patch inputs: patch612
`f5ecaf2b…`, patch612dd `3fb0e5b1…`, patch_builder630 `5173d789…`, module `7e6b3e52…`, patch_skip644 `65fab4e9…`.

| file | sha (12 hex) |
|---|---|
| scripts/def_day.sh | `4599060940e0` |
| scripts/def_sameline.py | `5cd3ab0726d8` (= [537757]'s after-644-v3 value) |
| scripts/def_identity.py | `7073406f6aad` |
| scripts/def_day_helpers.py | `fb715e3833f3` |
| scripts/def_pipeline.py | `76314af6b249` |
| scripts/october_recipe.py | `9dbcd9149020` |
| scripts/make_submission_pm_cd4.py (builder) | `fe4678c31b3c` |
| scripts/def_preflight.py | `17f0836fd1e1` |
| scripts/def_postbuild.py | `557c2b522b56` |
| experiments/524-october-recipe-dr/gate524.py | `b16e60c08e48` |
| scripts/def_identity_prior.py | `7aeb66b44f28` |
| src/vcc2026/rank_up.py | `7e6b3e52e2ad` |

[537757]'s def_day value (ca84075d) was computed with 630dd v3, before the FAILSAFE -1 hardening. The value above includes
630dd v4. KA630 rule: the line `630dd: uprn FAILSAFE_TOTAL <n>` is present AND n == 0.

## PINNED 2026-10-02 ~14:30Z (authoritative until superseded): 630dd v5 (d4cd5f: a context with > 10 % fail-safe targets is logged "630 EFFECTIVELY OFF in context <c> (<n>/<N>)" and repeated as a WARNING before def_day's DONE; never fatal; test630dd 54/54)
`work/compose606_630dd5.log`: **ALL PASS**. The only change from 14:20Z is **def_day.sh `54fcf04070ca`**. Unchanged:
- def_sameline `5cd3ab07`
- builder `fe4678c3`
- def_identity_prior `7aeb66b4`
- rank_up `7e6b3e52`
- def_identity 7073406f, helpers fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914, def_preflight 17f0836f,
  def_postbuild 557c2b52, gate524 b16e60c0

604f (10-08+, NOT in the frozen chain) is re-anchored on 644 v3's `${SLSKIP:+--skip-contexts "$SLSKIP"}` call line.
The chain + 604f gives ALL PASS (`work/compose606_604f_v3.log`).

## PINNED 2026-10-02 14:53Z (commit fc60a1fa; authoritative until superseded): 604's flip-scale fix (dcef3ca7)
`work/compose606_flipscale.log`: **ALL PASS**. The only change from 14:30Z is **def_identity_prior.py `bb4a99fc5529`**
(was 7aeb66b4). Without a flip it is numerically the same function: test604flip F1 rebuilds every committed artifact
byte for byte. Unchanged:
- def_day `54fcf040`
- def_sameline `5cd3ab07`
- builder `fe4678c3`
- rank_up `7e6b3e52`
- def_identity 7073406f, helpers fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914, def_preflight 17f0836f,
  def_postbuild 557c2b52, gate524 b16e60c0

## PINNED 2026-10-02 14:56Z (commit 85d4be68; authoritative until superseded): 604 module 428f4a32 (flip = base + octile sides only)
`work/compose606_optb.log` **ALL PASS**. Only **def_identity_prior.py `428f4a323964`** moves (was bb4a99fc). Without a
flip it is the same function: test604flip F1 rebuilds every committed artifact byte for byte. Unchanged: def_day
`54fcf040`, def_sameline `5cd3ab07`, builder `fe4678c3`, rank_up `7e6b3e52`, and the rest as at 14:30Z.

## P1 / P2 for the coordinator's "10-06 APPLY DAY v2" ([4d0207], 2026-10-02 ~16:50Z)
**P1 (B2, test infrastructure).** Every chain-tree builder now SKIPS a patch whose `--check` says it is already in the tree
(rc ≠ 0 and "already" in the output; "partially" still fails). After the 10-05 apply HEAD carries 538 / 554 / 590 / 641,
and after 10-06 the rest. Builders changed (26c8eddb, cc4bb2cf):
- compose606.py;
- test630dd.tree;
- test604f.chain_tree;
- test604.variant_samelines;
- test_skip644.chain_tree;
- tree642.make_tree, plus 642's T1, whose pre-590/641 sub-check skips itself when the live tree carries 641.
- verify644.py (the 10-08 runbook re-verify; the coordinator, 3ea6debe): the same rule, 27/27 on today's HEAD, and all six
  steps report "already" on a fully patched tree.

Every patch's re-check was verified to say "already" (641, 598, 612, 612dd, 630, 630dd, 642, 604, 644, 604f, and
538 / 554 / 590). This is the "skip if --check says already" form the coordinator offered as the alternative to moving
641 into test598.BUNDLE. It needs no change to the canonical tree, and today's behaviour is unchanged: on today's HEAD
there are no skips, the shas are identical to the 14:56Z pin, and test630dd is 54/54.

**Proof on a simulated post-10-05 HEAD.**
- Setup: a scratch worktree at 26c8eddb with 538 → 554 → 590 → 641 applied to its live files and committed there
  (472798a2, detached, never pushed). def_sameline is the canonical post-bundle 22c6f0f0. Gitignored data and the main
  venv are linked in.
- compose606: **SKIP 641 (on HEAD), ALL PASS, final shas identical to the PINNED 14:56Z block** (def_day 54fcf040,
  def_sameline 5cd3ab07, builder fe4678c3, def_identity_prior 428f4a32, rank_up 7e6b3e52; the rest as pinned).
- Suites: test630dd **54/54**; test_ka630_audit **6/6**, and **11/11** after the capacity-floor gate (9e239a41, re-run
  in the worktree 17:08Z); test_skip644 **4/4**; test_hek293_642 **6/6** (T1's pre-tree sub-check skipped post-10-05,
  as designed); through 602's gate: **test604 35/35** (17:02Z, peak 6.61 GB) and **test604f 24/24** (17:07Z, peak
  3.85 GB). **P1 is proven: every suite builds and passes on a post-10-05 HEAD.** Logs in the session scratchpad
  (`wt1005_logs/`), not kept.

**P2 (KA630 audit).** `experiments/630-rank-neutral-up/ka630_audit.py <build>.rc_audit.jsonl [--out-dir DIR]`.
- It runs 630dd's OWN report snippet, extracted at run time from `patch630dd.C_NEW`. test_ka630_audit A1 asserts the
  snippet is verbatim in the composed def_day.sh.
- It prints the lines exactly as def_day logs them. Exit 0 ONLY with `630dd: uprn FAILSAFE_TOTAL 0`; exit 1 for n > 0,
  n = −1 (missing, unreadable or no ur_ rows) or no line. test_ka630_audit 6/6.
- **Capacity floor (added 2026-10-02 17:03Z, d4cd5f).** It now also exits 1 unless every context in `--contexts`
  (default A,B,C) is present and meets 630's prereg floor (mean ur_frac_got ≥ 0.127 AND mean ur_n_moved ≥ 500 per
  target). It prints `floor <ctx>: …` for each context and `FLOOR FAIL <ctx>` on a miss. The values come from this run's
  `uprn <ctx>: {…}` lines, never from a uprn_audit.json an earlier run left. test_ka630_audit **11/11**: floor miss on
  frac_got, on moved, a context absent, exactly at the floor (passes), and a stale passing json ignored.
- The KA630 build must set `VCC_RC_AUDIT_DUMP=<path>.jsonl` (remove a stale file first: the builder appends) together
  with `VCC_DE_UPRN`.
