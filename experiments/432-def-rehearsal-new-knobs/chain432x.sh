#!/usr/bin/env bash
# 432 EXTRAS chain (2026-09-23 ~12:25Z). Waits for the main chain (chain432.sh, pid in work/chain432.pid) to exit, then
# runs build432x.sh's extra phases strictly in order, each behind build432x.sh's own gate:
#   r:  dump_b0env                         the LOAD known answer (b0's own-env dump; 13 s)
#   ef: dump_stack, b1, level_b1, audit_b1 one low-coverage build (238's essential-free pessimistic bundle) with every
#                                          knob but the re-level: key, ballast, sharpen, stack LOAD, coarse MDE, zero
#                                          repair, zs
# RESUMABLE like chain432.sh (a phase with a "[<name>] phase done" line in state432.log is skipped). No pgrep -f.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=$REPO/experiments/432-def-rehearsal-new-knobs; W=$X/work; LOG=$W/chain432x.log; STATE=$W/state432.log
HEADB=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
bsha(){ shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1; }
MAIN="$(cat "$W/chain432.pid")"
log "START extras (pid $$); waiting for the main chain pid $MAIN to exit"
while kill -0 "$MAIN" 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609231800 ] && { log "main chain still running at 18:00Z -- giving up"; exit 1; }; sleep 30; done
log "main chain exited; its last line: $(tail -1 "$W/chain432.log")"
[ "$(bsha)" = "$HEADB" ] || { log "builder $(bsha) != HEAD $HEADB -- STOP"; exit 1; }
done_already(){ grep -qE " \[$1\] phase done$" "$STATE" 2>/dev/null; }
run(){ local crit="$1" set="$2" phase="$3" name="$3"; [ "$set" = ef ] && name="ef_$3"
       case "$phase" in level_*|audit_*) [ "$set" = ef ] && name="ef_$phase";; esac
       if done_already "$name"; then log "phase $set:$phase already done -- skipped"; return 0; fi
       log "phase $set:$phase START"; SET="$set" bash "$X/build432x.sh" "$phase"; local rc=$?; log "phase $set:$phase rc=$rc"
       if [ $rc -ne 0 ] && [ "$crit" = 1 ]; then log "critical phase $set:$phase failed -- STOP"; exit 1; fi; }
run 0 r  dump_b0env
run 1 ef dump_stack
run 1 ef b1
run 0 ef level_b1
run 0 ef audit_b1
log "builder at END: $(bsha)"
log "END"
