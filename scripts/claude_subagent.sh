#!/usr/bin/env bash
# Launch a research subagent as a HEADLESS Claude Code run (Fable 5.1, max effort by default).
#
# Why this exists (2026-09-02): Ahmet's Claude subscription is the funded runtime for
# subagents; Cursor's third-party-model credits are not. AGENTS.md §11.2 requires the
# strongest model at the highest effort for every subagent here, so this defaults to `--model fable
# --effort max` (override with CLAUDE_SUBAGENT_MODEL / CLAUDE_SUBAGENT_EFFORT). Every run logs its full stream-json event trace next to the experiment it
# owns, so a later session (or a transfer to yet another runtime) can read exactly what it did.
#
# Usage:
#   scripts/claude_subagent.sh <experiment-dir-name> <mandate-file.md>
#   scripts/claude_subagent.sh 164-callset-artefact experiments/164-callset-artefact/MANDATE.md
#
# Monitor:   tail -f experiments/<name>/claude-run.jsonl | grep -o '"name":"[A-Za-z]*"'
# Final:     python3 -c 'import json,sys;[print(json.loads(l).get("result","")) for l in open(sys.argv[1]) if l.startswith("{\"type\":\"result\"")]' experiments/<name>/claude-run.jsonl
# Stop:      kill "$(cat experiments/<name>/claude-run.pid)"
#
# Requires `claude auth status` -> loggedIn: true. If not: `claude auth login` (browser, once).
set -euo pipefail
NAME="${1:?experiment dir name, e.g. 164-callset-artefact}"
MANDATE="${2:?path to the mandate .md}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
[ -f "$MANDATE" ] || { echo "no mandate at $MANDATE" >&2; exit 2; }
mkdir -p "experiments/$NAME"
LOG="experiments/$NAME/claude-run.jsonl"
ERR="experiments/$NAME/claude-run.err"
PIDF="experiments/$NAME/claude-run.pid"
if [ -f "$PIDF" ] && kill -0 "$(cat "$PIDF")" 2>/dev/null; then
  echo "already running: pid $(cat "$PIDF")" >&2; exit 3
fi

# Guardrails repeated at the system level so a long run cannot forget them mid-way.
BUILD_CLAUSE="NEVER build or write a challenge submission package; NEVER write into submissions/;"
if [ "${CLAUDE_SUBAGENT_ALLOW_BUILD:-0}" = "1" ]; then
  BUILD_CLAUSE="You ARE authorised to build a submission package into submissions/ for this mandate \
only, and you MUST verify every pre-committed falsifier before reporting it as buildable;"
fi
SYS="You are a research subagent in the repo at $REPO for the Arc Virtual Cell Challenge 2026. \
Read AGENTS.md first and obey its rules, then docs/FINDINGS.md sections 0-2 (current state, current conclusions, closed families). Absolute rules: NEVER run 'vcc submit'; $BUILD_CLAUSE NEVER run 'cell-eval2' (one instance \
only, 25 GB RSS, the machine is shared); use 'uv run python', never bare python; keep peak RSS \
under ~10 GB; write ONLY inside the experiment directory named in your mandate; do not modify \
files outside it. Prefer measured numbers over reasoning; carry a known-answer row in every \
screen; a well-supported negative is a valid result. Finish by writing RESULTS.md in your \
experiment directory and then a final message summarising the numbers."

nohup command claude -p "$(cat "$MANDATE")" \
  --model "${CLAUDE_SUBAGENT_MODEL:-fable}" --effort "${CLAUDE_SUBAGENT_EFFORT:-max}" \
  --dangerously-skip-permissions \
  --append-system-prompt "$SYS" \
  --output-format stream-json --verbose \
  > "$LOG" 2> "$ERR" &
echo $! > "$PIDF"
echo "launched $NAME  pid=$(cat "$PIDF")  log=$LOG"
