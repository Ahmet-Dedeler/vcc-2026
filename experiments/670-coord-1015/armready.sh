#!/usr/bin/env bash
# 670 (coordinator [400b6e], 2026-10-15): generic READY-line slot waiter. Waits for "<ts> READY <TAG> <sha> <pkg>" in a ready log,
# reads the package path from it, then execs arm1015.sh for slot <SL> behind <prior-log>.
#   nohup bash experiments/670-coord-1015/armready.sh <SL> <TAG> <ready-log> <prior-log> <FLOOR> <CAP> "<name>" "<note>" > …/armready-<SL>.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SL="${1:?slot}"; TAG="${2:?tag}"; RDY="${3:?ready log}"; PRIOR="${4:?prior log}"; FLOOR="${5:?floor}"; CAP="${6:?cap}"; NAME="${7:?name}"; NOTE="${8:?note}"
W=experiments/670-coord-1015/work; LOG=$W/armready-$SL.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/armready-$SL.pid; log "armed (pid $$): waiting for READY $TAG in $RDY (cap $CAP)"
until L="$(/usr/bin/grep -E " READY $TAG [0-9a-f]{64} " "$RDY" 2>/dev/null | tail -1)"; [ -n "$L" ]; do
  [ -f "$W/NOGO-$SL" ] && { log "NOGO-$SL present"; log "END"; exit 1; }
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "no READY $TAG by CAP"; log "END"; exit 1; }
  sleep 60; done
PKG="$(echo "$L" | awk '{print $5}')"; log "READY $TAG seen: pkg $PKG"
exec bash experiments/670-coord-1015/arm1015.sh "$SL" "$TAG" "$RDY" "experiments/428-zero-gene-audit/work/build_$TAG/run$TAG.log" "$PKG" "$NAME" "$NOTE" "$PRIOR" "$FLOOR" "$CAP"
