#!/usr/bin/env bash
# Supervise headless research subagents: relaunch any that died on a TRANSIENT API error.
#
# Why this exists (2026-09-03): four of six agents died on `API Error: 529 Overloaded` and
# `API Error: 500 Internal server error` and sat dead for an hour because nobody was watching.
# A transient provider failure must not cost a run. A genuine completion must not be relaunched.
#
# Discriminator: the run's final `{"type":"result"}` event. `is_error: true` AND the result text
# matching an API/overload/timeout pattern => transient, relaunch. Anything else (including a
# clean NO-GO, or the [bio] classifier which is deterministic) => leave it alone.
#
# Usage:  scripts/agent_supervisor.sh <exp-dir> [<exp-dir> ...]
# Emits one stdout line per action, so it can be driven by a Monitor.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
MAXTRIES="${SUPERVISOR_MAX_TRIES:-3}"

classify() {  # $1 = jsonl path -> prints TRANSIENT | DONE | RUNNING_NO_RESULT
  python3 - "$1" <<'PY'
import json,sys,re
p=sys.argv[1]; res=None
try:
    for l in open(p):
        try: o=json.loads(l)
        except Exception: continue
        if o.get("type")=="result": res=o
except FileNotFoundError:
    print("RUNNING_NO_RESULT"); raise SystemExit
if res is None: print("RUNNING_NO_RESULT"); raise SystemExit
txt=str(res.get("result",""))
transient=re.search(r"API Error:\s*(429|5\d\d)|Overloaded|Internal server error|Connection error|timed? ?out|ECONNRESET|socket hang up", txt, re.I)
print("TRANSIENT" if (res.get("is_error") and transient) else "DONE")
PY
}

for NAME in "$@"; do
  D="experiments/$NAME"; PIDF="$D/claude-run.pid"; LOG="$D/claude-run.jsonl"
  [ -d "$D" ] || continue
  # still alive?
  if [ -f "$PIDF" ] && kill -0 "$(cat "$PIDF")" 2>/dev/null; then continue; fi
  # find the mandate this agent was launched with
  MAN="$D/MANDATE-CONTINUE.md"; [ -f "$MAN" ] || MAN="$D/MANDATE.md"
  [ -f "$MAN" ] || { echo "SKIP $NAME (no mandate)"; continue; }
  [ -f "$D/PAUSED.md" ] && continue
  [ -f "$D/SUPERSEDED.md" ] && continue
  # It wrote up: never relaunch. NOTE (2026-09-04): an agent that writes RESULTS.md early and
  # keeps iterating (221 did) is therefore UNGUARDED for the rest of its run. That is deliberate --
  # relaunching a finished agent is worse than losing a tail -- but it means a long-running agent
  # should write RESULTS.md last, and a mandate that says "draft it as you go" trades supervision
  # for safety against context loss. Prefer drafting to a scratch name and renaming at the end.
  [ -f "$D/RESULTS.md" ] && continue
  case "$(classify "$LOG")" in
    TRANSIENT)
      TRIES=$(cat "$D/.supervisor_tries" 2>/dev/null || echo 0); TRIES=$((TRIES+1))
      if [ "$TRIES" -gt "$MAXTRIES" ]; then echo "GIVEUP $NAME after $MAXTRIES relaunches"; continue; fi
      echo "$TRIES" > "$D/.supervisor_tries"
      mv -f "$LOG" "$D/claude-run.attempt$TRIES.jsonl" 2>/dev/null || true
      rm -f "$PIDF"
      echo "RELAUNCH $NAME (transient API error, attempt $TRIES/$MAXTRIES)"
      CLAUDE_SUBAGENT_ALLOW_BUILD="${CLAUDE_SUBAGENT_ALLOW_BUILD:-0}" \
      CLAUDE_SUBAGENT_MODEL="${CLAUDE_SUBAGENT_MODEL:-fable}" \
        bash scripts/claude_subagent.sh "$NAME" "$MAN" >/dev/null 2>&1 || echo "RELAUNCH-FAILED $NAME"
      ;;
    DONE) : ;;
    *) : ;;
  esac
done
