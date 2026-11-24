#!/usr/bin/env bash
# Resume 534's Part B worker session for Part C (MANDATE-C.md).
# Writes claude-run-C.{jsonl,err,pid}; never touches Part B's claude-run.jsonl.
# Run it only after Part B's process has exited, and outside the windows in PREREGISTRATION-C §6.
set -euo pipefail
REPO="~/Code/virtual-cell-challenge-2026"
cd "$REPO"
D=experiments/534-generalization-audit
SID="f96b201b-7dfd-43a1-8452-3ae065928cb3"   # Part B's session (claude-run.jsonl's init event)
if [ -f "$D/claude-run.pid" ] && kill -0 "$(cat "$D/claude-run.pid")" 2>/dev/null; then
  echo "Part B worker still running (pid $(cat "$D/claude-run.pid"))" >&2; exit 3
fi
if [ -f "$D/claude-run-C.pid" ] && kill -0 "$(cat "$D/claude-run-C.pid")" 2>/dev/null; then
  echo "Part C already running (pid $(cat "$D/claude-run-C.pid"))" >&2; exit 3
fi
# Same guardrails as scripts/claude_subagent.sh (default: no builds).
SYS="You are a research subagent in the repo at $REPO for the Arc Virtual Cell Challenge 2026. \
Read AGENTS.md first and obey its rules, then docs/FINDINGS.md sections 0-2 (current state, current conclusions, closed families). Absolute rules: NEVER run 'vcc submit'; NEVER build or write a challenge submission package; NEVER write into submissions/; NEVER run 'cell-eval2' (one instance \
only, 25 GB RSS, the machine is shared); use 'uv run python', never bare python; keep peak RSS \
under ~10 GB; write ONLY inside the experiment directory named in your mandate; do not modify \
files outside it. Prefer measured numbers over reasoning; carry a known-answer row in every \
screen; a well-supported negative is a valid result. Finish by writing the results file your mandate names (RESULTS-C.md; never RESULTS.md, which is the owner's) in your \
experiment directory and then a final message summarising the numbers."
nohup command claude -p "$(cat "$D/MANDATE-C.md")" --resume "$SID" \
  --model opus --effort max --dangerously-skip-permissions \
  --append-system-prompt "$SYS" --output-format stream-json --verbose \
  > "$D/claude-run-C.jsonl" 2> "$D/claude-run-C.err" &
echo $! > "$D/claude-run-C.pid"
echo "launched 534 Part C  pid=$(cat "$D/claude-run-C.pid")  log=$D/claude-run-C.jsonl"
