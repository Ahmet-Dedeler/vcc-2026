#!/usr/bin/env bash
# lab_subagent.sh — launch (or resume) a lab worker as a HEADLESS Claude Code run (the lab's copy of the production launcher).
# Differences from the production launcher: the system prompt states the lab's boundaries (docs/LAB-2026-10-16.md incl.
# Amendment 1) instead of production rules that do not apply here; heavy compute (incl. the real scorer) runs ONLY on EC2
# through experiments/800-checker/cloud/labcloud.sh. Model: Opus 5.5 (AGENTS 11.2: tolerates dataset/biology-heavy
# mandates), effort high (charter). Detached (setsid) so an app restart cannot kill it.
# Charter Amendment 1: no production gate pattern may appear in a lab process's argv, so this prompt names none of the
# production chain scripts; the worker reads them from the charter file instead.
#   bash experiments/800-checker/lab_subagent.sh <experiment-dir-name> <mandate.md>
#   bash experiments/800-checker/lab_subagent.sh --resume <experiment-dir-name> <message.md>    (same session, context kept)
set -euo pipefail
RESUME=0
if [ "${1:-}" = "--resume" ]; then RESUME=1; shift; fi
NAME="${1:?experiment dir name}"; MANDATE="${2:?mandate or message path}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
[ -f "$MANDATE" ] || { echo "no file at $MANDATE" >&2; exit 2; }
mkdir -p "experiments/$NAME"
LOG="experiments/$NAME/claude-run.jsonl"; ERR="experiments/$NAME/claude-run.err"; PIDF="experiments/$NAME/claude-run.pid"
if [ -f "$PIDF" ] && [ "$(ps -o comm= -p "$(cat "$PIDF")" 2>/dev/null | tr -d ' ')" = "claude" ]; then echo "already running: pid $(cat "$PIDF")" >&2; exit 3; fi
# reserved windows (window800.sh): at most 2 lab agents, the lead included, so at most ONE worker
if ! W800="$(bash experiments/800-checker/window800.sh)"; then
  live=0; for p in experiments/8[0-9][0-9]-*/claude-run.pid; do [ -f "$p" ] && kill -0 "$(cat "$p")" 2>/dev/null && live=$((live + 1)); done
  [ "$live" -lt 1 ] || { echo "refused: ${W800}; ${live} lab worker(s) already running (max 1 + the lead)" >&2; exit 5; }
fi
RES=()
if [ "$RESUME" -eq 1 ]; then
  SID="$(grep -o '"session_id":"[^"]*"' "$LOG" | tail -1 | cut -d'"' -f4)"
  [ -n "$SID" ] || { echo "no session id in $LOG" >&2; exit 4; }
  RES=(--resume "$SID")
fi
[ -f "$LOG" ] && cp "$LOG" "$LOG.$(date -u +%Y%m%dT%H%M%SZ).bak"
SYS="You are a research worker in the VCC lab (experiments 800-899) of the repo at $REPO (Arc Virtual Cell Challenge 2026). \
Read docs/LAB-2026-10-16.md (the lab charter, INCLUDING its Amendment 1, which binds you), then AGENTS.md, then \
docs/FINDINGS.md sections 0-2, then experiments/800-checker/DESIGN.md. Absolute rules: write ONLY inside the experiment \
directory your mandate names; never edit production code (scripts/, src/, the builder, the October day script, the recipe \
module), shared docs, other experiments' files, anything under data/, waiters or slots; never use the challenge's \
submission CLI; never build a challenge submission; no git commands (the lab lead commits). Never put the names of the \
production chain scripts listed in the charter's Amendment 1 into a command line, a script name or a file you launch \
(production gates match process arguments). The Mac runs production builds: locally you may only read files and run \
sub-second aggregations under 100 MB; any other local job goes through experiments/602-council-gate/gate602.sh (<= 8 GB; it \
refuses to start beside a production build) and never while 'bash experiments/800-checker/window800.sh' prints WINDOW (the \
reserved windows: daily 23:30-01:45Z, and from the D/E/F release until R0 is submitted); in a window also no transfer > 1 GB \
from the Mac; everything heavy (operator fits, dumps, emitter runs, ANY real-scorer run) runs on EC2 \
through experiments/800-checker/cloud/labcloud.sh (stage / launch / wait), on-demand, smallest instance and minutes cap \
that fit, within your mandate's cloud budget and the launcher's caps (EC2 <-> S3 transfers are fine any time); check 'labcloud.sh status' before you finish and leave no instance running. Use 'uv run python', \
never bare python. Carry a known-answer row in every screen and abort on a failed one; a well-supported negative is a \
valid result. Record paths and sha256, never descriptions. Finish by writing RESULTS.md in your experiment directory \
and then a final message summarising the numbers, the known answers and what is left."
nohup python3 -c 'import os, sys; os.setsid(); os.execvp(sys.argv[1], sys.argv[1:])' \
  claude -p "$(cat "$MANDATE")" "${RES[@]+"${RES[@]}"}" --model "${LAB_WORKER_MODEL:-opus}" --effort "${LAB_WORKER_EFFORT:-high}" \
  --dangerously-skip-permissions --append-system-prompt "$SYS" --output-format stream-json --verbose \
  > "$LOG" 2> "$ERR" &
echo $! > "$PIDF"
echo "launched $NAME pid=$(cat "$PIDF") log=$LOG resume=${RESUME}"
