#!/usr/bin/env bash
# 670 (coordinator [223d39], 2026-10-19): one-command recovery of the coordinator's standing waiters after a reboot (AGENTS §6.2 step 2).
# IDEMPOTENT: each item is launched only if no copy is running, so running this twice (or beside another session) never doubles anything.
# Covers the final-phase set: arm1015 R0, name_guard R0, guard_release, board_snapshot, pisces_watch, caffeinate. Readiness's lane
# (watch691, r0ready692, wait708) is in OCT22-CHECKLIST §B. A/B/C slot waiters are re-armed by hand only while that slot is still open.
#   bash experiments/670-coord-1015/rearm_coord.sh          (prints what it launched and what was already running)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/670-coord-1015/work
running(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^(bash|caffeinate)$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
launch(){ local what="$1"; shift
  if running "$what"; then echo "running already: $what"; return; fi
  nohup bash scripts/detach.sh "$@" > /dev/null 2>&1 & sleep 1; echo "launched: $what"; }
# R0's waiter, re-armed 10-22 05:28Z on the re-gate's logs (r0g-ready; the first arm stopped on r0's NOT READY). Floor 10-19 00:00Z, cap 10-27.
launch "arm1015.sh R0" bash experiments/670-coord-1015/arm1015.sh R0 r0def experiments/692-def-r0/work/r0g-ready/r0def.ready.log \
  experiments/692-def-r0/work/r0g-ready/r0def.run.log submissions/r0def.vcc "I am persistent but open-minded." \
  "R0 D/E/F (re-gated r0g, DEF_DAY_ALLOW_LOW_DOSE=1)" experiments/692-def-r0/work/prior-stub.log 202610190000 202610270000
launch "name_guard-R0b.spec" bash experiments/670-coord-1015/name_guard.sh $W/name_guard-R0b.spec 202610270000
launch "guard_release.sh" bash experiments/670-coord-1015/guard_release.sh
launch "scripts/board_snapshot.sh" bash scripts/board_snapshot.sh
launch "scripts/pisces_watch.sh" bash scripts/pisces_watch.sh
running "caffeinate -is" && echo "running already: caffeinate" || { nohup bash scripts/detach.sh caffeinate -is -t 1600000 > /dev/null 2>&1 & sleep 1; echo "launched: caffeinate"; }
echo "--- then: readiness's lane (OCT22-CHECKLIST §B: watch691, r0ready692, wait708) if that session is not back."
