#!/usr/bin/env bash
# Second slot of the day: wait until the FIRST armed submission's log shows a successful submit and an
# entry id, then hand off to arm_submit_after.sh, which waits for that entry to be terminal before firing.
#   scripts/arm_second_slot.sh <first-slot-log> <YYYYmmddHHMM-UTC floor> <package.vcc> <sha256> "<joke name>" <logfile>
set -uo pipefail
FIRSTLOG="${1:?first-slot log}"; FLOOR="${2:?floor}"; PKG="${3:?package}"; SHA="${4:?sha}"; NAME="${5:?name}"; LOG="${6:?log}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "ARMED-SECOND waiting for a successful first submit in $FIRSTLOG; pkg=$PKG name='$NAME' pid=$$"
while :; do
  if grep -q "SUBMIT rc=0" "$FIRSTLOG" 2>/dev/null; then
    ID="$(sed -n '/SUBMITTING/,$p' "$FIRSTLOG" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | head -1)"
    [ -n "$ID" ] && { log "first slot entry id: $ID"; break; }
    log "first submit succeeded but no entry id found yet; retrying"; sleep 60
  elif grep -q "ERROR" "$FIRSTLOG" 2>/dev/null; then
    log "first slot logged an ERROR -- the second slot is NOT armed automatically; lead decides"; exit 1
  fi
  sleep 60
done
exec bash scripts/arm_submit_after.sh "$ID" "$FLOOR" "$PKG" "$SHA" "$NAME" "$LOG"
