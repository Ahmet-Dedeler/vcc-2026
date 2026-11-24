#!/bin/bash
# power602.sh -- power guard for gate602 ([4d0207] for [c941b4], 2026-10-03 10:45Z: the Mac ran on battery at 55 %).
# While `pmset -g batt` says "Battery Power", it keeps work/quiet_hold 15 min in the future. Every gate602 >= v6 waiter
# refuses to START while that hold is active (jobs already running are untouched). When AC returns, it removes the
# hold it set, so starts resume on the next 10 s check. If this watcher dies, its last hold lapses within 15 min
# (fail-open, like a normal quiet hold).
#   usage: nohup bash experiments/602-council-gate/power602.sh > /dev/null 2>&1 &     (one instance; pid in work/power602.pid)
#   test:  POWER602_W=<dir> POWER602_PMSET=<script printing pmset-like output> POWER602_ONCE=1 bash power602.sh
set -u
W=${POWER602_W:-~/Code/virtual-cell-challenge-2026/experiments/602-council-gate/work}
PM=${POWER602_PMSET:-pmset}
LOG=$W/gate602.log; echo $$ > "$W/power602.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
STATE=unknown
while :; do
  B=$("$PM" -g batt 2>/dev/null)
  if printf '%s\n' "$B" | head -1 | grep -q "Battery Power"; then
    echo $(( $(date +%s) + 900 )) > "$W/quiet_hold"
    [ "$STATE" = battery ] || log "POWER battery ($(printf '%s\n' "$B" | sed -n 2p | awk '{print $3,$4,$5}')): gate602 starts held (power602 pid $$)"
    STATE=battery
  elif printf '%s\n' "$B" | head -1 | grep -q "AC Power"; then
    if [ "$STATE" != ac ]; then
      [ "$STATE" = battery ] && { rm -f "$W/quiet_hold"; log "POWER AC: hold lifted, gate602 starts resume"; }
      [ "$STATE" = unknown ] && log "POWER AC at start (power602 pid $$)"
      STATE=ac
    fi
  fi
  [ "${POWER602_ONCE:-0}" = 1 ] && break
  sleep 30
done
