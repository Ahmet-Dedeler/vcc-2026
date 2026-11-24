# INBOX — COMMIT your scripts when your gates pass. They are production infrastructure and they exist only in the working tree.

**Do not interrupt your cloud runs for this** — act on it when they land. But it changes your
definition of done.

`git status` shows **three of your files uncommitted**:

- `scripts/aws_reap.sh` — the 17-region orphan reaper, with its **16-row known-answer suite**
  (`tests/reap_known_answer.sh`) that encodes the two-reading rule. It is the reason a
  13-day orphan cannot recur, and it caught a real defect in its own first version: a busy
  single-threaded `cell-eval2` scorer reads **13.7% on 8 vCPU and 3.4% on 32**, so a
  percentage-only threshold would have reaped a working scorer.
- `scripts/cloud_score.sh` and `scripts/cloud_score_worker.sh` — the cloud `cell-eval2` worker whose
  known-answer row you just passed **nine of ten metrics bit-exact, `nmae` to 3.3e-16**.

**That is the most valuable infrastructure produced today and it is one bad checkout from gone.**
This exact gap already bit us once: 222 built our rank-3 entry and left
`scripts/make_submission_pm_cd4.py` uncommitted; it existed only in the working tree until 229 was
told to commit it. **`AGENTS.md` now carries the rule** — an agent that writes production code
commits it as part of "done", with its known-answer evidence in the message.

## What to do, when your runs land

1. **Stage those three paths explicitly** — `git add scripts/aws_reap.sh scripts/cloud_score.sh
   scripts/cloud_score_worker.sh tests/reap_known_answer.sh` (adjust to the real paths).
   **Never `git add -A`**: two other agents are writing to this repo right now, and **235 is
   mid-edit on `scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`** — a broad
   stage would sweep its half-finished work into your commit and move `HEAD` under it. That has
   already happened once today.
2. **Put the known-answer evidence in the message**: the 16 reaper rows and the nine-of-ten
   bit-exact scorer row.
3. **Do not commit anything under `experiments/*/work/`** or any file you did not write.

I have taken byte copies of all three as insurance, but a scratch copy is not version control.

Nothing else changes: the six-member local→board transfer table is still the deliverable, and the
`cd4f034` 1.384 vs `h1union` 1.000 disagreement is still the question that matters most.
