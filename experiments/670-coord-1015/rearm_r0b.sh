#!/usr/bin/env bash
# 670 (coordinator [3c7a7e], 2026-10-22): the coordinator's half of readiness's plan B for R0 (experiments/690-def-readiness-1015/planb_r0.sh).
# Run it only AFTER planb_r0.sh has launched def_day in experiments/692-def-r0/work/r0b and r0ready692 on OUT692=…/r0b-ready.
# It re-arms R0's submit waiter on plan B's fresh ready / run logs (the first arm stopped on r0's NOT READY and never submitted):
#   1. refuses if a bash `arm1015.sh R0` is alive, or if plan B's r0ready692 is not running;
#   2. moves the first attempt's experiments/670-coord-1015/work/submit-R0.log aside (submit-R0.attempt1.log), so the new arm and
#      anything reading that path start from a clean log;
#   3. launches arm1015 R0 with the SAME package path, name, note, prior stub, floor and cap as the first arm (5181), on
#      ready-log …/r0b-ready/r0def.ready.log and run-log …/r0b-ready/r0def.run.log, and name_guard R0 on a matching spec.
# wait708: in the die case it has already run day708 on the first arm's END (its decisions do not depend on R0's package), so it
# is NOT re-armed here. If day708 did not run, re-arm wait708 by hand with planb_r0.sh's printed command after this script.
#   DRY=1 bash experiments/670-coord-1015/rearm_r0b.sh     (prints, touches nothing)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
DRY="${DRY:-0}"; W=experiments/670-coord-1015/work; OUTB="${R0B_READY_DIR:-experiments/692-def-r0/work/r0b-ready}"
# R0B_READY_DIR: the dir whose r0def.ready.log / r0def.run.log the new arm reads (plan B r0b-ready, or a re-gate's fresh ready dir).
R0B_BUILD_DIR="${R0B_BUILD_DIR:-experiments/692-def-r0/work/r0b}"   # the r0ready692 to require alive (skip with R0B_NO_692=1)
RDY=$OUTB/r0def.ready.log; RUNL=$OUTB/r0def.run.log; PKG=submissions/r0def.vcc
NAME="I am persistent but open-minded."; NOTE="${R0B_NOTE:-R0 D/E/F (plan B r0b)}"; PRIOR=experiments/692-def-r0/work/prior-stub.log
FLOOR=202610190000; CAP=202610270000
say(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*"; }
live(){ ps -axo comm=,args= | awk -v p="$1" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
if live "arm1015.sh R0 "; then
  [ "$DRY" = 1 ] && say "DRY: an arm1015 R0 is alive now (a real run would STOP here)" || { say "STOP: an arm1015 R0 is alive"; exit 1; }
fi
[ "${R0B_NO_692:-0}" = 1 ] || live "r0ready692.sh data/challenge-def $R0B_BUILD_DIR" || [ "$DRY" = 1 ] \
  || { say "STOP: no r0ready692 on $R0B_BUILD_DIR (set R0B_NO_692=1 if the READY line is written another way)"; exit 1; }
[ -s "$PRIOR" ] || { say "STOP: prior stub $PRIOR missing"; exit 1; }
SPEC=$W/name_guard-R0b.spec
if [ "$DRY" = 1 ]; then
  say "DRY: mv $W/submit-R0.log $W/submit-R0.attempt1.log"
  say "DRY: arm1015.sh R0 r0def $RDY $RUNL $PKG \"$NAME\" \"$NOTE\" $PRIOR $FLOOR $CAP"
  say "DRY: name_guard.sh $SPEC $CAP  (spec: R0 r0def $RDY $RUNL $PKG \"$NOTE\" $PRIOR $FLOOR $CAP)"
  exit 0
fi
[ -e "$W/submit-R0.log" ] && mv "$W/submit-R0.log" "$W/submit-R0.attempt1.log" && say "moved the first attempt's log to submit-R0.attempt1.log"
printf '%s\n' R0 r0def "$RDY" "$RUNL" "$PKG" "$NOTE" "$PRIOR" "$FLOOR" "$CAP" > "$SPEC"
nohup bash scripts/detach.sh bash experiments/670-coord-1015/arm1015.sh R0 r0def "$RDY" "$RUNL" "$PKG" "$NAME" "$NOTE" "$PRIOR" "$FLOOR" "$CAP" \
  > experiments/692-def-r0/work/arm-R0b.nohup 2>&1 &
sleep 3
nohup bash scripts/detach.sh bash experiments/670-coord-1015/name_guard.sh "$SPEC" "$CAP" > $W/name_guard-R0b.nohup 2>&1 &
sleep 3
say "armed: $(tail -1 $W/submit-R0.log 2>/dev/null | cut -c1-160)"
say "name_guard: $(tail -1 $W/name_guard-R0.log 2>/dev/null | cut -c1-120)"
