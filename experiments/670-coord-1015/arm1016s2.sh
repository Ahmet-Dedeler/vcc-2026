#!/usr/bin/env bash
# 670 (coordinator [400b6e], 2026-10-15): 10-16 slot 2 = K-UP (ku30). Reads the package path from its READY line in
# experiments/701-kup-callset/work/READY-science.log ("<ts> READY ku30 <sha> <pkg>"), then execs arm1015.sh 1016s2 after slot 1.
#   nohup bash experiments/670-coord-1015/arm1016s2.sh "<name>" > experiments/670-coord-1015/work/arm1016s2.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
NAME="${1:?name}"; W=experiments/670-coord-1015/work; LOG=$W/arm1016s2.log; RDY=experiments/701-kup-callset/work/READY-science.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/arm1016s2.pid; log "armed (pid $$): waiting for READY ku30 (cap 202610162100)"
until L="$(/usr/bin/grep -E ' READY ku30 [0-9a-f]{64} ' "$RDY" 2>/dev/null | tail -1)"; [ -n "$L" ]; do
  [ -f "$W/NOGO-1016s2" ] && { log "NOGO-1016s2 present"; log "END"; exit 1; }
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610162100 ] && { log "no READY ku30 by 21:00Z 10-16"; log "END"; exit 1; }
  sleep 60; done
PKG="$(echo "$L" | awk '{print $5}')"; log "READY ku30 seen: pkg $PKG"
exec bash experiments/670-coord-1015/arm1015.sh 1016s2 ku30 "$RDY" experiments/428-zero-gene-audit/work/build_ku30/runku30.log "$PKG" "$NAME" "ku30 slot 2 on 10-16" "$W/submit-1016s1.log" 202610160001 202610162300
