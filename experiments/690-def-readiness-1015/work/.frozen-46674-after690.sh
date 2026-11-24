#!/usr/bin/env bash
# 690 launcher (10-15 21:35Z, readiness [4592ee]): starts run690 (R-NEW + R-KA) after the slot builds, without itself
# looking like run690. Science's build705 refuses to start while any bash process names run690.sh, and run690 must start
# after build705, so arming run690 directly would deadlock both. This waits until BOTH:
#   - the coordinator's 690 HOLD file ($W/HOLD) is gone (it stays until build705 exits), and
#   - none of the queued chains is alive under bash: ka693b, rabc694, build701, R7a's r7a_access_1015b, build705
# then execs run690.sh (CAP690 default 10-19 00:00Z; run690 re-checks HOLD, builders, memory and disk itself).
# Waiting stops at AFTER690_CAP (default 10-18 18:00Z). Log: work/after690.log. Terminal line: "<ts> exec run690 ..." (then
# run690.log takes over) or "<ts> END" on a stop.
# usage: nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/after690.sh >> experiments/690-def-readiness-1015/work/after690.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work; LOG=$W/after690.log; HOLD=$W/HOLD
# Run from a frozen per-launch copy, so an edit of this file never corrupts a live run (bash reads scripts by offset)
if [ -z "${AFTER690_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-after690.sh"; cp "$REPO/$E/after690.sh" "$FZ" || exit 1
  AFTER690_FROZEN=1 exec bash "$FZ" "$@"
fi
CAP="${AFTER690_CAP:-202610181800}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"' EXIT
# bash + a path substring in the args, never pgrep -f (AGENTS 6.2); substrings also match the chains' frozen copies
alive(){ ps -axo comm=,args= | awk -v p="$1" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
CHAINS="experiments/693-dnmag-speed/ka693b.sh rabc694.sh experiments/701-kup-callset/build701.sh r7a_access_1015b.sh experiments/705-lodo-strength-safety/build705.sh"
echo $$ > "$W/after690.pid"
log "after690 armed (pid $$): run690 starts once the 690 HOLD is gone and none of [$CHAINS] is alive; cap $CAP"
last=""
while :; do
  LIVE=""; for c in $CHAINS; do alive "$c" && LIVE="$LIVE ${c##*/}"; done
  [ -e "$HOLD" ] && LIVE="$LIVE HOLD"
  [ -z "$LIVE" ] && break
  [ "$LIVE" != "$last" ] && log "waiting for:$LIVE"; last="$LIVE"
  if [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; then
    log "STOP: still waiting at cap $CAP (for:$LIVE)"
    # science's sci1016 (build707) waits for a NEW "END" in run690.log: give it one, so a run690 that never starts cannot
    # hold build707 until its own cap. Only on this terminal stop (a kill for a re-arm must not release build707).
    printf '%s %s\n%s END\n' "$(date -u +%FT%TZ)" "after690 STOP at cap $CAP: run690 never started (for:$LIVE)" \
      "$(date -u +%FT%TZ)" >> "$W/run690.log"
    exit 3
  fi
  sleep 60
done
trap - EXIT
log "exec run690 (CAP690 ${CAP690:-202610190000}); its log: $W/run690.log"
CAP690="${CAP690:-202610190000}" exec bash "$REPO/$E/run690.sh"
