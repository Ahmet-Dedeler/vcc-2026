You are an adversarial code reviewer for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read
AGENTS.md first (rules), then docs/FINDINGS.md section 0, docs/DEF-RUNBOOK.md ("THE DAY IN ONE PAGE" and section 6), and
the bottom of docs/NEXT-SLOT-DECISION.md (the 2026-09-28 entries). This VM has NO gitignored data (no challenge controls,
donors, harness matrices): review code, and write reproductions that run on synthetic inputs you construct.

Experiment 519. Write ONLY under experiments/519-cloud-october-review/. Commit to your branch only. Never run `vcc`.

Scope — the changes made on 2026-09-28 and the October path they feed (Oct 22 is a one-shot run with a dark board, so a
silent defect costs the final):
1. scripts/def_postbuild.py: the "axis moved off A/B/C's" rows now key on the TARGET set (abc_targets). Is the new
   condition correct for every case (new contexts + new targets; new contexts + A/B/C's targets; A/B/C itself; a panel that
   shares most but not all A/B/C targets)? Can it fail open (pass a build whose operator is for the wrong panel)?
2. experiments/502-october-de-riders/: patch_defday502.py and repro/def_day.sh.new / def_day_helpers.py.new (the
   DEF_DAY_DE_MAG / DEF_DAY_MDE_SIGN pass-throughs), and chain502.sh + commit502.spec (applies them tonight after the build
   queue). Check: exports reach the level builds and never the dump; the allowed-config-diff additions; riders_match's
   logic; the final.env lines; chain502's release condition, failure/restore paths, known-answer checks, and whether any
   path can leave def_day.sh half-patched or commit a wrong file.
3. experiments/497-october-switches/finish497dr.sh (a manual continuation of chain497's step b) and gate497.py: does the
   gate still verify what chain497 would have verified?
4. experiments/456-zs-riders/arm0928.sh (a cap edited 12:00Z -> 23:30Z) and the re-armed slot chains 0929..1005
   (experiments/457-*/chain0929.sh, 459-*/chain0930c.sh, 464-*/chain1001b.sh, 479-*/chain1002.sh, 495-*/chain1003.sh,
   498-*/chain1004.sh, 496-*/chain1005.sh): with the 09-28 slots having fired 12 h late (12:11Z and 12:38Z), can any
   downstream chain misread which entries exist, fire into the wrong UTC day, or double-submit?
5. scripts/october_recipe.py and docs/october-decisions.json with de_mag now DECIDED (dn_k250_zs): does `command` produce
   a def_day invocation that actually passes the rider (it needs 502's pass-through, not yet applied)? What does it do
   today, before 502 lands — refuse, or silently drop the rider?

Known traps to check for (AGENTS 6.2): macOS /bin/bash 3.2 (no associative arrays), BSD grep ($ before \| is ignored;
use grep -E), process matching by pattern (must match interpreter + path), waiters that require a success-only END
line, editing a script a live bash waiter is executing, zsh non-splitting of unquoted variables.

Deliverable: experiments/519-cloud-october-review/RESULTS.md — for each finding: file:line, the failure scenario, a
reproduction script in your dir that demonstrates it on synthetic inputs (or a precise argument why it cannot be
reproduced without data), severity (blocks October / costs score silently / cosmetic), and the minimal fix as a patch
file (NOT applied to production files). Also list what you checked and found correct. Be concrete; no generic advice.
