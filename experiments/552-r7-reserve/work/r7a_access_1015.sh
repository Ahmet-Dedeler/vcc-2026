#!/usr/bin/env bash
# R7a access (science session [bf003b], coordinator-approved 10-14): starts r7a551.py access only when no def_pipeline build
# is in its build stage (def_pipeline.other_builders() empty) AND 551's own gate admits it (it gates every heavy step at
# C10b's floors: 22 GiB free+inactive, pressure 55 %, swap 400 MiB, inside the recorded window). Resumable: a re-run skips
# done steps. Cap: the recorded window's end (10-16 12:00Z). Terminal line: "<ts> END rc=<rc>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "armed (pid $$): waiting for no build in its build stage"
until ! uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610161200 ] && { log "window ended before a builder-free moment"; log "END rc=9"; exit 9; }
  sleep 60
done
log "no build running: starting r7a551.py access"
uv run python experiments/551-r6-c10b/r7a551.py access > experiments/551-r6-c10b/work/r7a/access.log 2>&1
rc=$?
log "access rc=$rc: $(tail -3 experiments/551-r6-c10b/work/r7a/access.log | tr '\n' ' ' | cut -c1-400)"
log "END rc=$rc"
