#!/bin/bash
# print cited lines: args like R162:9-10  BR222:6  A:49  RB:16-20  BL:8-9  SL:31  DP:84  MB:946
res() { case "$1" in
  A) echo AGENTS.md;; RB) echo docs/DEF-RUNBOOK.md;; BL) echo docs/BOARD-LEDGER.md;; SL) echo submissions/SUBMITTED_LEDGER.json;;
  DP) echo scripts/def_pipeline.py;; MB) echo scripts/make_submission_pm_cd4.py;; CS) echo docs/CURRENT_STATE.md;;
  R*) ls -d experiments/${1#R}-*/RESULTS.md 2>/dev/null | head -1;;
  BR*) ls experiments/${1#BR}-*/BOARD-RESULT*.md 2>/dev/null | head -1;;
esac; }
for c in "$@"; do key=${c%%:*}; rng=${c#*:}; f=$(res "$key"); a=${rng%-*}; b=${rng#*-};
  echo "--- $c  ($f)"; sed -n "${a},${b}p" "$f" | cut -c1-230; done
