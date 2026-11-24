#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16): prior guard for slot <SL>. arm1015 makes the NEXT slot wait until submit-<SL>.log has a terminal
# line. If <SL>'s armready ends without a package (no READY by its cap, or a NOGO), it never execs arm1015, so submit-<SL>.log never gets
# one and the next slot would sit until ITS cap. This guard watches armready-<SL>.log for an END written after the guard started, with no
# arm1015/armready for <SL> alive and no terminal line in submit-<SL>.log. Then it appends "NOT SUBMITTING (prior_guard) …" + "END" there,
# with no pid and no entry id, so the next slot "goes alone".
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/prior_guard.sh <SL> <CAP YYYYmmddHHMM> > …/prior_guard-<SL>.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SL="${1:?slot}"; CAP="${2:?cap}"; W="${PG_W:-experiments/670-coord-1015/work}"
A="$W/armready-$SL.log"; S="$W/submit-$SL.log"; LOG="$W/prior_guard-$SL.log"; T0="$(date -u +%FT%TZ)"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
newend(){ awk -v t0="$T0" 'substr($0, 1, 20) > t0 && $0 ~ / END$/ {f = 1} END {exit f ? 0 : 1}' "$A" 2>/dev/null; }
live(){ ps -axo comm,args | awk -v s=" $SL " '$1 ~ /^bash$/ && /(arm1015|armready)\.sh/ && index($0, s) {n++} END {exit n ? 0 : 1}'; }
terminal(){ /usr/bin/grep -qE " END$|SUBMIT rc=|NOT SUBMITTING|ERROR " "$S" 2>/dev/null; }
echo $$ > "$W/prior_guard-$SL.pid"; log "armed (pid $$, T0 $T0): watching $A; writes a terminal line to $S if armready ends without a package"
while :; do
  terminal && { log "$S is terminal already: nothing to do"; log "END"; exit 0; }
  if newend && ! live; then
    printf '%s NOT SUBMITTING (prior_guard): armready-%s ended without a package: %s\n' "$(date -u +%FT%TZ)" "$SL" \
      "$(/usr/bin/grep -v ' END$' "$A" | tail -1 | cut -c22-140)" >> "$S"
    printf '%s END\n' "$(date -u +%FT%TZ)" >> "$S"
    log "wrote the terminal line to $S"; log "END"; exit 0
  fi
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "cap $CAP reached"; log "END"; exit 0; }
  sleep 30
done
