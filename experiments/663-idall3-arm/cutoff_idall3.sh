#!/usr/bin/env bash
# One-shot (coordinator, 2026-10-04): if idall3's build has not STARTED (no "DRIVER invocation 1" in its run log) by
# 22:30Z, kill it so that recipe_dr524 (10-05 02:00Z; it needs no other build and free+inactive >= 25 GiB) and 10-06's
# 524dr slot are never delayed by a fill.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/663-idall3-arm; LOG=$X/work/build_idall3.log; RL=experiments/428-zero-gene-audit/work/build_idall3/runidall3.log
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610042230 ]; do sleep 60; done
if /usr/bin/grep -q "DRIVER invocation 1" "$RL" 2>/dev/null; then exit 0; fi
for p in $(ps -axo pid,args | awk '/bash/ && /(build_idall3\.sh|run_arm_generic\.sh idall3)/ {print $1}'); do kill "$p" 2>/dev/null; done
printf '%s %s\n' "$(date -u +%FT%TZ)" "NOT READY cut off at 22:30Z: the memory gate never opened (keeps the machine clear for recipe_dr524)" >> "$LOG"
