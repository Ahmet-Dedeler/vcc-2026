# 594 — The 10-05 sequence composes: 538 → 554 → 590 apply cleanly in that order, every suite passes at every step, and the final shas are fixed below. One hazard: run_arm_generic.sh pins the pre-538 driver sha

2026-10-01, session "Explore contribution areas [57577c]", at the coordinator's request (due 10-04 12:00Z). Scratch trees
only. The **live files were untouched**: the 11 touched files were hashed before and after, and were equal. Class:
operations. Harness: `compose594.py`. Record: `work/compose594.json`, `work/compose594.log`. The first pass is kept as
`work/compose594.pass1.*`; see "Harness bugs" below.

## How
- **Scratch S.** A git worktree of HEAD `b1a5dffa` (detached, in this session's scratchpad), removed at the end.
- **Overlays.** The gitignored inputs the suites read were symlinked into S read-only: 319 links (control bundles incl.
  the R7b *control* bundle, not its truth file; the CCLE cache; donor arrays; rehearsal metas; the two audit dumps;
  `.venv`).
- **What the patches write.** Only the 11 files below, all real tracked files in S.
- **Order.** In the 10-05 order, with S's own copies of the patch scripts.
- **Suites.** 538's and 554's suites copy the live-side files and apply their own patch(es) themselves. They are
  pre-application suites, so they run at step 0, which is today's HEAD (= live on all 11 files). 590's suite runs at
  steps 0, 1 and 2, i.e. against the tree it will meet on 10-05.
- **Pre-538 `def_day.sh`.** It hard-codes the live repo path (531 #20). It was never executed; only `bash -n` and the
  patches touched it.

## Results
| step | apply | checks (bash -n def_day.sh · py_compile of every touched .py · `october_recipe.py check --allow-pending`) | suites | determinism / extra |
|---|---|---|---|---|
| 0 (HEAD b1a5dffa) | — | PASS · PASS · PASS (pending: `lodo_hi_pct`, `sharpen_pt`, `mde_sign`) | 538 `test_extras538` **PASS** (27 s); 554 `test_guard554` **PASS** (8 s); 554 `identity554` **PASS** (61 s); 590 `test_identity590` **PASS** (21 s) | S == live on all 11 files |
| 1: 538's `apply_october538.py` | `--check` rc 0, apply rc 0 | PASS · PASS · PASS (the new switches land PENDING) | 590 suite against the post-538 tree **PASS** | shas == a fresh reference application of 538 to copies of the live files |
| 2: 554's `patch_guard554.py` (applies by default; `--check` first) | `--check` rc 0, apply rc 0 | PASS · PASS · PASS; `def_preflight.py --verify` **PASS** | 590 suite against the post-538+554 tree **PASS** | shas == `test_guard554.make_tree(with554=True)` on all 10 files it builds |
| 3: 590's `patch_identity590.py` | `--check` rc 0, apply rc 0; a second application refused | PASS · PASS · PASS; `def_identity.py --verify` **PASS**, printing the lowest top-1 correlation 0.772 / 0.719 / 0.733 | gate cases on the applied `def_sameline.py` vs the step-2 (pre-590) file: a (registered, unanimous, 0.554) off ✓; b (0.70) on, npz identical ✓; c (no key) on with WARNING, identical ✓; d (`--min-top1 0`) on, identical ✓; f (not unanimous) off ✓ | — |

## Final sha256 of every file the three patches touch (after 538 → 554 → 590)
| file | HEAD (today) | after 538 → 554 → 590 | changed by |
|---|---|---|---|
| `scripts/def_day.sh` | `0b53ba89bd47…` | `5851f84b9f238c932816d22fbe6cd810c6db68c863db6bdb41d8746d10932f6c` | 538, 554 |
| `scripts/def_day_helpers.py` | `20142142e457…` | `e3ff45bcefc1a6d3af651126903726e51972c61343386752f36be19bfd8215f3` | 538 |
| `experiments/448-def-rehearsal-q175/relevel_iter.py` | `bd10f8f25361…` | `e360a1d9907801a1b5f67832a1e0922a10cf12f605b765eaabbd38739bf71b35` | 538 |
| `scripts/def_sameline.py` | `19596af1c906…` | `22c6f0f08c1f02883887a34690bc5642f172026433d90b19275a578f24d86259` | 538, 590 |
| `scripts/october_recipe.py` | `a8654095d924…` | `503c45a00bf74164a206f13cfcd878c468d6d5d8783a1a04068c32f42bb2c793` | 538 |
| `docs/october-decisions.json` | `e41b165ea00c…` | `8ee9e74dee8ec0c84de2e8d263f4bc2cc05288ca7c112feda57a10c11e82e3db` | 538 |
| `experiments/524-october-recipe-dr/gate524.py` | `277361b87527…` | `7a43d2f5ae0f03c017a6903b9d524e6e7f9352de0f32f4641b0d3a4dba0cf5ae` | 538 |
| `scripts/def_postbuild.py` | `35278c7ad1fd…` | `557c2b522b565636261ea297cc4dcc2e5e1af8f1ea0bbc8c788790388f0b68f6` | 538, 554 |
| `scripts/def_pipeline.py` | `9b2c7eca6c6b…` | `9a718b9eaa24f10cfbf255ce830ea4a748f61252ec503bbb6db69d3ddf5ca9fd` | 538 |
| `scripts/def_preflight.py` | `d95ebbb45f88…` | `2d1811945ebeb506cf2bde739c03c57f5350e84311a4476d11337a2355b1622d` | 554 |
| `scripts/def_identity.py` | `2d71b2595ed8…` | `c98a608fe2838d1ca6b68dbf090874548f286c911f472becf8ef3dcfcc157f59` | 590 |

These hold only if the 11 files are still at the HEAD shas above on 10-05. `recipe_dr524.sh` runs first and requires
`def_day`, `helpers`, the builder and the emitter at HEAD, so check all 11 against this table before applying.

## Ordering (all six permutations on fresh temp trees of the live files)
| order | result | final == canonical |
|---|---|---|
| 538 → 554 → 590 (canonical) | all apply | — |
| 538 → 590 → 554 | all apply | **yes** |
| 590 → 538 → 554 | all apply (590 on the pre-530 `def_sameline.py`, `e49c0c4c…`; 530's edits then apply on top) | **yes** |
| 554 first (554 → 538 → 590, 554 → 590 → 538), or 590 → 554 → 538 | **554 refuses** (rc 1: its `def_day.sh` anchor is 538's text) | — |
- 590 commutes with both 538 and 554. 554 must follow 538, by design, and fails closed otherwise (nothing written).

## Hazard found: `run_arm_generic.sh` pins the pre-538 driver
- `experiments/428-zero-gene-audit/run_arm_generic.sh` accepts the driver only by sha, through its `case "$D" in …` list
  (5463fe28 / 44714b9a / 62b5110d / 9b2c7eca).
- 538's applier changes `scripts/def_pipeline.py` from `9b2c7eca` to **`9a718b9e`**.
- **After 10-05, the first A/B/C arm built with `run_arm_generic.sh` stops with "DRIVER hash … stopping" (exit 4).** This
  covers the 10-07 arms if they are built that way (likely 10-06). It fails closed, but it costs a build window.
- Fix, as chain476 did: add `9a718b9eaa24f10cfbf255ce830ea4a748f61252ec503bbb6db69d3ddf5ca9fd` to the case list (atomic
  rename, no run active). Do it in the 10-05 sequence, after the inert KA shows the post-538 driver rebuilds an A/B/C arm
  unchanged.
- Other pins of the touched files' current shas are historical scripts that must not be re-run anyway: `def_day_helpers`
  `20142142` in `502/chain502.sh` (550: never re-run), `517/patch_defday517.py`, and the 519 / 550 repro scripts.
  Searched: every tracked `*.sh` / `*.py`, by sha prefix.

## Harness bugs (pass 1, `work/compose594.pass1.*`; fixed before the clean run, none a patch defect)
1. **The 554 call was wrong.** `patch_guard554.py` applies by default and has no `--apply`, so my call exited 2 (argparse)
   and 554 was never applied in pass 1.
2. **One overlay path was missing.** identity554 also reads the R7b control bundle (`552/work/r7b/bundle`), which the
   overlay lacked.
3. **590's own T3 had a test defect.** It required `--check`'s post-530 line, which `--check` prints only when the live
   `def_sameline.py` is pre-530. After 538 the file is already post-530, so T3 failed while T1/T2 passed. Fixed in
   `test_identity590.py` (b1a5dffa): T3 accepts a post-530 live file, whose first `--check` line already is the post-530
   proof.

## Commands (coordinator, 10-05)
- Before: compare the 11 live files against the HEAD column.
- Run the sequence (538's applier; 554's `patch_guard554.py --root .` with `--check` first; 590's
  `patch_identity590.py --check`, then `--apply`).
- After: **`uv run python experiments/594-def-bundle-composition/post594.py`** (read-only on the live tree, ~30 s). It
  checks:
  - the 11 shas against the final column;
  - `def_identity.py --verify`, with the three lowest-top-1 lines printed;
  - the 590 gate cases on the live `def_sameline.py` against the pre-590 file, rebuilt from git HEAD + 538's applier in a
    temp tree (sha `03cd591a…`).
- `post594.py` was itself validated on a scratch worktree with all three patches applied: **ALL PASS**, live untouched.
- On today's unapplied live tree it fails by design (no patch is in yet).
- 590's own suite is a pre-application suite (it patches in memory and refuses an already-patched file). Run it before
  the apply, not after.

## Update (2026-10-02): post594.py no longer depends on a frozen sha table
- The coordinator's 10-02 decision (`sharpen_pt` DECIDED) changed `docs/october-decisions.json` on HEAD. Its post-bundle
  sha is now `9da65575…`, not 594's `8ee9e74d…`.
- `october-decisions.json` takes decisions daily until 10-05, so a frozen table goes stale.
- `post594.py` [1] now rebuilds the expected 11 files from **git HEAD's copies + 538 → 554 → 590** in a temp tree and
  compares the live files to that. 594's table is printed for reference.
- This needs HEAD to be the pre-bundle commit when it runs, which holds at checklist step 6 (the commit is step 9).
- Validated: ALL PASS on a scratch worktree of HEAD 6315d037 with the three patches applied, live untouched. The other
  10 files still match 594's table exactly.
- Checklist step 0's "files still at 594's HEAD shas" should read **"live == HEAD on the 11 files"**
  (`git diff --quiet -- <the 11 paths>`).

