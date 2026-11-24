#!/bin/bash
# 819: run the lab's judge on one proposal (experiments/800-checker/JUDGE.md). The judge is a fresh headless run that reads
# and writes only experiments/819-judge/<idea>/VERDICT.md; the lab lead copies the verdict to the idea's dir.
#   bash experiments/819-judge/judge.sh <idea-dir-name> [proposal-file]   e.g. 810-strong-direction PROPOSAL-2.md
# The proposal is copied (frozen) into the judge's dir before the judge starts, so later edits cannot change what is judged.
set -euo pipefail
R=~/Code/virtual-cell-challenge-2026; cd "$R"
I="${1:?idea dir}"; PF="${2:-PROPOSAL.md}"; [ -f "experiments/$I/$PF" ] || { echo "no experiments/$I/$PF" >&2; exit 2; }
N="819-judge/$I-$(date -u +%Y%m%dT%H%M%SZ)"; mkdir -p "experiments/$N"
cp "experiments/$I/$PF" "experiments/$N/PROPOSAL-frozen.md"
cat > "experiments/$N/MANDATE.md" <<EOM
# Judge run for experiments/$I/$PF, frozen copy experiments/$N/PROPOSAL-frozen.md (sha256 $(shasum -a 256 experiments/$N/PROPOSAL-frozen.md | cut -c1-16))
You are the lab's judge. Follow experiments/800-checker/JUDGE.md exactly. Judge the FROZEN copy experiments/$N/PROPOSAL-frozen.md (never the live file), with the idea's MANDATE.md and any earlier
VERDICT-*.md in experiments/$I/ (a resubmission must fix every item of the previous verdict),
experiments/800-checker/PHASEB.md (what the checker can and cannot price), docs/LAB-LOG.md's two "Known-answer report" sections,
docs/FINDINGS.md §1–§2, and every experiment RESULTS.md the proposal or the closures cite; grep experiments/*/RESULTS.md for the
mechanism's key words. Compute nothing heavier than a lookup. Write experiments/$N/VERDICT.md: GO / REVISE / KILL, then each
JUDGE.md check answered with the file and line that settles it, then (for REVISE) exactly what is missing. Be strict: a
plausible idea without a falsifier, or one that re-names a closed construction, does not get compute. Then stop.
EOM
LAB_WORKER_EFFORT=high bash experiments/800-checker/lab_subagent.sh "$N" "experiments/$N/MANDATE.md"
echo "judge for $I launched in experiments/$N"
