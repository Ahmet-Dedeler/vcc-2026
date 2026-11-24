#!/usr/bin/env bash
# R7a resume (R7 owner [d80a10], 10-16 ~07:30Z; coordinator [2d6f46] placed R7a after readiness's R-NEW and before build707).
# The 06:34Z resume failed before any panel; the fix and the ledger's driver pin are 27780ccb / 45e51ed2 (test_access552 PASS).
# Trigger, by env (the coordinator's choice):
#   (a) default: run690.log has a NEW " END" line after this waiter started, and neither run690 nor after690 is alive;
#   (b) R7A_TRIGGER_RE=<ERE>: run690.log has a NEW line matching it (e.g. readiness's "R-NEW done" line, with run690 paused
#       at the coordinator's HOLD before R-KA), and no bash def_day.sh / run_arm_generic.sh chain is alive.
# Then R7a's v2 wrapper (r7a_access_1015b.sh: no def_pipeline build in its build stage, 551's gate per heavy step).
# build707.sh lists THIS waiter in its CHAINS, so sw707 cannot start beside R7a.
# Start cap: the ledger R7 section's window extension minus 30 min (an edited extension line moves it; fallback 10-16 17:30Z).
# Release guard: RELEASE-DETECTED before the start -> skip. Terminal line, always: "<ts> END".
#   [R7A_TRIGGER_RE='...'] nohup bash scripts/detach.sh bash experiments/552-r7-reserve/work/r7a_after690.sh > experiments/552-r7-reserve/work/r7a_after690.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/552-r7-reserve/work/r7a_after690.log
T0="$(date -u +%FT%TZ)"; TRE="${R7A_TRIGGER_RE:-}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
r7cap(){ local w c
  w="$(awk '/^## Reserve R7/{s=1;next} s&&/^## /{s=0} s' docs/CONFIRMATION-LEDGER.md \
       | /usr/bin/grep -oE 'R7a WINDOW END EXTENDED TO [0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}Z' | tail -1 | awk '{print $NF}')"
  c="$(date -j -u -v-30M -f '%Y-%m-%dT%H:%MZ' "$w" +%Y%m%d%H%M 2>/dev/null)"; echo "${c:-202610161730}"; }
PL=experiments/690-def-readiness-1015/work/run690.log
ready(){ if [ -n "$TRE" ]; then newline "$T0" "$PL" "$TRE" && ! alive scripts/def_day.sh && ! alive run_arm_generic.sh
         else newline "$T0" "$PL" " END$" && ! alive run690.sh && ! alive after690.sh; fi; }
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0): trigger ${TRE:+(b) /$TRE/ in run690.log + no def_day/run_arm_generic}${TRE:-(a) run690 END + exit}; start cap $(r7cap)"
# the trigger must hold on two polls 60 s apart: a kill + relaunch of run690 (readiness did two at 07:18-07:21Z) leaves a
# gap in which neither process is alive
until ready && sleep 60 && ready; do
  CAP="$(r7cap)"
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "trigger not met by the start cap $CAP: R7a cannot finish inside its window; NOT started (coordinator decides)"; exit 0; }
  sleep 10
done
log "trigger met: $(tail -1 "$PL" | cut -c1-160)"
[ -e experiments/691-release-watcher/work/RELEASE-DETECTED ] && { log "R7a skipped: release detected"; exit 0; }
log "starting R7a's wrapper"
bash experiments/552-r7-reserve/work/r7a_access_1015b.sh > experiments/552-r7-reserve/work/r7a_access_1015b.nohup 2>&1
log "R7a wrapper exited rc=$?: $(tail -2 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log | tr '\n' ' ' | cut -c1-300)"
