#!/usr/bin/env bash
# Drive agent_supervisor.sh every 3 minutes over every experiment dir that has a mandate
# and no RESULTS.md. Runs as a daemon so agent supervision does not depend on a Monitor
# firing. Logs every action (the supervisor is silent when nothing needs doing).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
LOG="${SUPERVISOR_LOG:-work/supervisor.log}"; mkdir -p "$(dirname "$LOG")"
while :; do
  NAMES=()
  for d in experiments/*/; do
    n="${d%/}"; n="${n##*/}"
    [ -f "$d/MANDATE.md" ] || [ -f "$d/MANDATE-CONTINUE.md" ] || continue
    [ -f "$d/RESULTS.md" ] && continue
    [ -f "$d/PAUSED.md" ] && continue
    [ -f "$d/SUPERSEDED.md" ] && continue
    [ -f "$d/CANCELLED.md" ] && continue
    [ -f "$d/claude-run.pid" ] || continue      # never launched by us: leave alone
    NAMES+=("$n")
  done
  if [ ${#NAMES[@]} -gt 0 ]; then
    OUT="$(bash scripts/agent_supervisor.sh "${NAMES[@]}" 2>&1)"
    [ -n "$OUT" ] && printf '%s %s\n' "$(date -u +%H:%MZ)" "$OUT" >> "$LOG"
  fi
  sleep 180
done
