# 242 — Hunt every guard that can PASS while the thing it protects is WRONG. Four instances found by accident today.

You are a research subagent. Work in `experiments/242-guard-audit/`. **READ-ONLY over the repo. No
`vcc`, no `cell-eval2`, no build, no submission, no production change, nothing re-measured.** Peak
RSS ≤ 6 GB — this is reading. `uv run python`. **241 is running and shares this repo; do not write
outside your own directory.**

## The class, named from four instances found in one day — all by accident, none on any roadmap

**A check that passes because it is checking the wrong reference.**

| # | guard | how it passed while wrong | cost | status |
|---|---|---|---|---|
| 1 | `solve_cd4_fraction` / `plolo.py:363-369` | median over **all** panel rows → exactly `0` below 50% coverage, then **returns 0 instead of raising** | a real source reads as **identically absent** | **STILL LIVE** (223 diagnosed, never fixed; ~215 call sites) |
| 2 | the whitening operator's panel check | only `n_cols`, which is **17,942 on every panel** | **98.8% of coordinates mis-paired**, −0.021 to −0.0429 on `pds` | fixed by **235**, guard **observed firing** |
| 3 | `def_postbuild.py` | diffed against **`rc_callset`'s** meta — the incumbent from *before* both shipped arms | printed **14/14 PASS** on a **−0.0199671** submission | fixed by **240**, failure **reproduced verbatim** |
| 4 | rule 13's `reach` toll | keyed on a digest whose **own calibration points sit in the zero row** | the charged row contains **zero board points** | corrected by **230** |

**Four in one day, found by accident.** Nobody has ever looked for them on purpose.

## What to look for — four signatures

1. **Stale reference.** A guard comparing against a hardcoded hash, meta file, config or constant
   that names an artefact **superseded since**. `grep` for hardcoded shas, `--incumbent-*` defaults,
   `rc_callset`/`rc_hybrid`/`a15`/`a23`/`191`/`202` in *live* code paths, and any comment of the form
   *"byte-for-byte from …"*.
2. **Silent default instead of a raise.** `if x > 0 else 0.0`, `.get(k, default)`, bare `except`,
   `or 0`, `np.nan_to_num` — anywhere the fallback is **indistinguishable from a valid answer.**
   Signature #1 above is exactly this.
3. **A check that cannot detect its own failure mode.** `n_cols` was the worked example: constant
   across every panel it was meant to discriminate. **Ask of each guard: what value would it print
   if the thing it protects were broken?** If that is the same value it prints when healthy, it is
   not a guard.
4. **Assertions never observed to fire.** 235 and 240 were each required to **demonstrate** the
   raise. **Which existing guards have never been shown firing?** List them; that is the queue.

## Scope, in priority order
`scripts/` **first and in full** — it is the production and Oct-22 path (`make_submission_pm_cd4.py`,
`def_pipeline.py`, `def_postbuild.py`, `build_union_whitener.py`, `cloud_score*.sh`, `aws_reap.sh`,
`board_snapshot.sh`, `agent_supervisor.sh`, `supervisor_loop.sh`). Then the **shared harnesses** many
experiments import — `experiments/162-retrodicting-validator/plolo.py`, `153-de-bench`, and anything
else with ≥ 10 importers (find them). Experiment-local one-off scripts are **out of scope** unless
something imports them.

## Method
- **Rank by blast radius**: how many call sites, and does it touch the **Oct-22 path** or a
  **shipped** artefact? A defect in a one-off script matters far less than one in `scripts/`.
- **For each finding: quote the code with file:line, say what value it prints when healthy, what it
  prints when broken, and why those differ or do not.**
- **Do not fix anything.** 241 is running and the shared harness has ~215 call sites; a mid-run
  change would corrupt a live experiment and break known-answer rows across the archive. **Specify
  the fix and its known-answer requirement instead.**
- **Report the null honestly.** *"`scripts/` is clean apart from the known live defect"* is a real
  and valuable result.

## Deliverable
`RESULTS.md`: a ranked table — **guard, file:line, signature, blast radius, what it would fail to
catch** — then the top three in detail with the fix and the known-answer row each would need. Plus
the **list of guards never observed firing**, which becomes the demonstration queue.

**Fix nothing. Build nothing. Propose no slot.**
