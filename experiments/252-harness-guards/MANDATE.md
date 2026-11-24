# 252 — Close the shared-harness guards 242 found (#4, #5, #6, #13, #14) without moving a single published number

You are a fix-and-prove subagent. Work in `experiments/252-harness-guards/`. **The only files you may edit
outside it:** `experiments/162-retrodicting-validator/plolo.py` and `experiments/064-pm-cd4-composition/run.py`
(the shared harness and its operator loader; 173 and ~9 referencing files). No other agent is running tonight.
**No `vcc`, no build, no `cell-eval2`, no production file.** `uv run python`. Peak RSS ≤ 10 GB; memory gate
before fold loads (`memory_pressure` ≥ 40% free AND free+inactive ≥ 20 GiB); strictly serial; `/usr/bin/time -l`.
Commit each file with its known-answer hashes in the message (`git add <file>`, never `-A`).

## READ FIRST
`docs/READ-FIRST.md`; `experiments/242-guard-audit/RESULTS.md` §1 rows 4, 5, 6, 13, 14 and §2 (the class:
a check that passes because it checks the wrong reference); `AGENTS.md` §11.1 rule 10 (every guard must be
demonstrated firing) and the 2026-09-04 note on the live silent-zero defect in the harness (`docs/LOG.md`,
"Record the still-live silent-zero defect in the shared PLOLO harness, with the fix specified"); 223 §7.1
(the `f_per_unit = 0.0` mirror at `plolo.py:363-369`); 218 §7.1 (the carrying-row rule); 247 §3 (how the
builder now sizes); `experiments/246-direction-hypotheses/lib246.py` and `experiments/250-deflation-october/`
(the most recent heavy users of the harness — their known-answer rows are your regression set).

## The rule that binds this mandate above all others
**Every published number that went through this harness must reproduce bit-for-bit after your change** when
called the way it was called. The fixes are guards and explicit rules, not behaviour changes on the paths
existing experiments used. Where a fix would change a value an experiment depended on (the silent zero), keep
the old path reachable under an explicit argument that names it, make the default RAISE, and show which
published rows (if any) were computed through the zero path — 223 says 0 of 21 were void, so expect none.

## Items
- **#5 / #4** `plolo.Operator(path)` accepts any `(17942, 17942)` file, and the frozen-operator sha check in
  `064/run.py:179-183` compares a value to itself. Make `Operator` record and assert the operator's sha and its
  axis digest against a registry keyed by name (the builder's `WHITENER_REGISTRY` / 235's `axis_sha256`
  convention), and make the 064 check compare against a pinned expectation. Demonstrate: the three
  near-identical operators on disk (`8.036` / `6.257` / `5.895` max eigenvalue) are told apart; a swapped file
  raises; every existing caller that passes the correct file still runs.
- **#6** `build_fold`'s cache-validity test ignores `FOLDS[held]["sources"]`, `["op"]`, `PM_BLEND`, `MIN_CPM`
  and restores a stale `audit`. Key the cache on a digest of those inputs; demonstrate a changed source list
  invalidating the cache and the unchanged one hitting it, with the audit reflecting the live inputs.
- **#13** the on-axis preserve `float(np.median(...)) or 1.0` silently substitutes 1.0 at a zero median (the same
  ≤ 50%-coverage condition as the known defect, one function away). Raise instead, with the coverage in the
  message; demonstrate on a synthetic residual at 40% coverage.
- **#14** `endpoints()` sign accuracy reads 0.0 on an all-zero truth block via `max(live.sum(), 1)`. Return NaN
  (or raise) and make aggregation skip it explicitly; demonstrate on a synthetic all-zero block.
- **The live silent zero** (`plolo.py:363-369`): make the default raise on a residual whose all-rows median is
  zero while carrying rows exist, offer `sizing="carrying"` as the explicit alternative (218's rule), and keep
  `sizing="legacy_all_rows"` reachable by name for reproduction. Demonstrate on `telohaec` / `gse295214_h1` (223's
  rows: fraction exactly 0 today).

## Known-answer rows (abort on failure)
1. 218's `compose218` headline rows (`cd40.34+h1PAIRCARRY0.28`, three cancer folds: macro +7.6803%, `pds`
   +0.004887) and 246's five bridge rows (`STACK`, `B2_cd4_svd5_resolved`, the random control, the shuffled
   projector, `C2_partial`) reproduce at 0.0e+00 through the changed harness, called as those scripts call it.
2. 222's `repro222.py` and 249's PLOLO reproduction still pass unchanged (run them).
3. 233's `BASE` rows on `abc`/`real`/`pess` (its `analyse233.json`) reproduce at 0.0e+00.
4. 216's `f = 0.34` / `f = 0.10` `cos_agg` rows (`work/sweep216.json`) reproduce.

## Deliverable
`RESULTS.md`: for each item the diff summary, the guard shown FIRING against a byte copy of the pre-252 harness
and PASSING on the correct input, the regression table (every known-answer row with its number and 0.0e+00), the
list of call sites whose behaviour changed (expect: none on default paths), commit hashes, and one line:
**HARNESS GUARDED / NOT**. Then stop.
