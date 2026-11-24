#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-15): remove [17adab]'s 690 HOLD once build705 has STARTED and then EXITED (its old log carries a
# stale NOT READY from 12:55Z, so a log line alone is not a trigger). Builder order: … → R7a → build705 → run690. Cap 10-18 12:00Z.
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/release_hold690.sh > experiments/670-coord-1015/work/release_hold690.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/670-coord-1015/work/release_hold690.log; HOLD=experiments/690-def-readiness-1015/work/HOLD
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /experiments\/705-lodo-strength-safety\/build705\.sh/' | /usr/bin/grep -q .; }
running(){ ps -axo comm,args | awk '$1 ~ /^(python[0-9.]*)$/ && /make_submission_pm_cd4\.py/' | /usr/bin/grep -q .; }
cap(){ [ "$(date -u +%Y%m%d%H%M)" -ge 202610181200 ]; }
log "armed (pid $$): waiting for build705 to start building, then to exit; then rm $HOLD"
# 1. build705 must be alive AND have passed its start gates (a python builder running while it is alive)
until alive && running; do cap && { log "cap reached before build705 built; HOLD left in place"; log "END"; exit 1; }; sleep 60; done
log "build705 is building"
# 2. build705's process exits
while alive; do cap && { log "cap reached while build705 alive; HOLD left in place"; log "END"; exit 1; }; sleep 60; done
log "build705 exited: $(tail -1 experiments/705-lodo-strength-safety/work/build_lstr705.log | cut -c1-160)"
rm -f "$HOLD" && log "HOLD removed: run690 (R-NEW + R-KA) may start"
log "END"
