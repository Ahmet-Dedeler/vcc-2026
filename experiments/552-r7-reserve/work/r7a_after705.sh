#!/usr/bin/env bash
# R7a re-run (R7 owner [d80a10], 10-16): after build705 (lstr705) has written a NEW END and its process is gone, run R7a's
# single access through its v2 wrapper (r7a_access_1015b.sh: no def_pipeline build, 551's gate per heavy step, window end
# 10-16 12:00Z), with the admission-only gate551 deviation recorded in 552 RESULTS + the ledger (a83f3fe2).
# Coordinator's order: build705 -> R7a -> readiness's run690 (the coordinator keeps the 690 HOLD until this wrapper ENDs).
# Release guard: RELEASE-DETECTED before the start -> skip. Terminal line, always: "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/552-r7-reserve/work/r7a_after705.sh > experiments/552-r7-reserve/work/r7a_after705.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/552-r7-reserve/work/r7a_after705.log
T0="$(date -u +%FT%TZ)"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0): waiting for build705 to end"
BL=experiments/705-lodo-strength-safety/work/build_lstr705.log
until newline "$T0" "$BL" " END$" && ! alive experiments/705-lodo-strength-safety/build705.sh; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610161130 ] && { log "build705 still running at 11:30Z: R7a cannot finish inside its window; NOT started"; exit 0; }
  sleep 10
done
log "build705 ended: $(tail -2 "$BL" | head -1 | cut -c1-160)"
[ -e experiments/691-release-watcher/work/RELEASE-DETECTED ] && { log "R7a skipped: release detected"; exit 0; }
log "starting R7a's wrapper"
bash experiments/552-r7-reserve/work/r7a_access_1015b.sh > experiments/552-r7-reserve/work/r7a_access_1015b.nohup 2>&1
log "R7a wrapper exited rc=$?: $(tail -2 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log | tr '\n' ' ' | cut -c1-300)"
