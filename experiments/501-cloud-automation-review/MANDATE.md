# 501 — adversarial review of the slot and build automation armed on 2026-09-27 (cloud agent, code only)

Repo `Ahmet-Dedeler/virtual-cell-challenge-2026`, branch `master`. Read `AGENTS.md` (§6.2 waiters: macOS `/bin/bash` 3.2,
BSD grep, the pgrep/argv trap, fire-on-process-exit) and `docs/NEXT-SLOT-DECISION.md` from "2026-09-27" to the end. No
gitignored data on this VM; work from code and synthetic inputs. Never run `vcc`, never touch credentials. Commit ONLY to
`cloud/501-automation-review`, ONLY inside `experiments/501-cloud-automation-review/`. No PR.

## Why
Four daily board slots (10-01 … 10-04) and a serial build queue now run unattended for ~48 h. A review this morning (499)
found that every per-target build would have died after building. That was fixed in the driver (`scripts/def_pipeline.py`
`september_rows`), and the queue was re-wired around it. Find anything else that would silently waste a slot, submit the
wrong package, deadlock the queue, or commit the wrong file.

## What to review
1. **Slot chains:** `experiments/479-pertarget-sharpen/chain1002.sh` + `arm1002.sh`, `experiments/495-lodo-cap6/chain1003.sh` +
   `arm1003.sh`, `experiments/498-combined-bulk/chain1004.sh` + `arm1004.sh` + `rule1004.py`. Check the readiness logic
   (package, gate lines, novelty, ρ band), the fallback order, that a package already scored can never be re-submitted,
   that s2 never equals s1, that the -m names (quotes with apostrophes, parentheses, an em dash, "!!") survive every hop
   to `vcc submit -m` and the ledger note, the time caps, and the entry-id extraction (`eid()` in chain1004).
2. **The build queue's gates:** `experiments/479-pertarget-sharpen/chain479b.sh` (+ `relaunch479b.sh`),
   `experiments/495-lodo-cap6/chain495.sh`, `experiments/476-def-day-merge/chain476.sh`, `experiments/497-october-switches/chain497.sh`
   + `chain497e.sh`, `experiments/498-combined-bulk/chain498.sh`, `experiments/496-uncalled-sign-hook/chain496p.sh`. For
   each wait loop: can it deadlock with another chain, a committer, the relauncher or a slot waiter (every long-lived bash
   process's argv is matched by these awk regexes)? Can two builds start at once? Does every failure path end with an
   `END` line and restore what it patched? Are the expected shas and `EXPECT` config diffs consistent with the patches
   they follow (`patch_builder495.py`, `patch_emitter496.py`, `patch_builder496.py`, `patch_defday497.py`)?
3. **The committer:** `scripts/commit_after_chain.sh` and the five `*.spec` files next to the chains (commit479/495/476/497
   and `experiments/496-uncalled-sign-hook/commit496.spec`). Spec parsing under bash 3.2, the process-exit wait, the
   RESTORED refusal, the REQ regexes against the lines each chain actually logs, `git commit --only`, and the index.lock
   retry.
4. **479's builder patch** (`experiments/479-pertarget-sharpen/patch_builder479.py`, function `sharpen_q_by_tercile` in
   `scripts/make_submission_pm_cd4.py`): the floor guard, tie handling, non-finite rows, and what happens on a panel whose
   key or strength has many exact ties or a huge tail (an October-like panel), including `digitize` edge semantics.

## Deliverable
`experiments/501-cloud-automation-review/RESULTS.md`: a severity-ranked table, then one section per finding (file:line,
what goes wrong and when, CONFIRMED with a command and its output / PLAUSIBLE / REFUTED, and a minimal fix as a patch under
`patches/`). Reproductions under `repro/`, runnable without the gitignored data. State which areas you checked and found
clean. Terse and concrete.

## Relaunch note (2026-09-27 07:15Z)
A first run of this mandate ended in ERROR before pushing anything. Its last message said the slot cap, the entry-id
lookup, the committer, the deadlines and the tercile tie rule have problems. Start from current `master`. **Commit and push
`RESULTS.md` after each finding** (a partial report is fine), so an error does not lose the work. Keep reproductions small.
