#!/usr/bin/env bash
# One-shot (coordinator, 2026-10-03; docs/NEXT-SLOT-DECISION.md "10-04 slot 2 FILL"): when chain1004 (pid 24456) exits after
# 20:00Z, and it armed 10-04 slot 1 but left slot 2 empty (rule1004's fill list is exhausted: 460dr was scored on 10-01),
# arm 498/arm1004.sh slot 2 with 446bj2zk20 (K20, READY, its own prereg). arm1004.sh slot 2 waits for slot 1's capture.
# It does nothing if chain1004 armed its own slot 2, or armed nothing at all. Cap: 10-04 00:30Z.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/660-slot-1007/work; LOG=$W/fill1004s2.log; C=experiments/498-combined-bulk/work/chain1004.log
A=experiments/498-combined-bulk
PKG=submissions/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc
RL=experiments/428-zero-gene-audit/work/build_446bj2zk20/run446bj2zk20.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/fill1004s2.pid
log "armed (pid $$): waiting for chain1004 (pid 24456) to exit"
while kill -0 24456 2>/dev/null; do
  if [ "$(date -u +%Y%m%d%H%M)" -ge 202610040030 ]; then log "cap 10-04 00:30Z reached with chain1004 still running -- nothing armed"; log "END"; exit 0; fi
  sleep 20
done
LINE="$(/usr/bin/grep -E 'armed 10-04: s1' "$C" | tail -1)"
log "chain1004 exited; its arm line: ${LINE:-<none>}"
case "$LINE" in
  "") log "chain1004 armed nothing (slot 1 empty): a slot-2 fill would wait for a slot 1 that never comes -- not arming"; log "END"; exit 0;;
  *"s2 <none>"*) ;;
  *) log "chain1004 armed its own slot 2 -- no fill needed"; log "END"; exit 0;;
esac
nohup bash $A/arm1004.sh 2 446bj2zk20 "$RL" "$PKG" "You're behind—very behind." "446bj2zk20 slot 2 on 10-04 (coordinator fill)" \
  > $A/work/arm1004_s2.nohup 2>&1 & echo $! > $A/work/arm1004_s2.pid
sleep 3
log "armed the K20 fill as 10-04 slot 2 (pid $(cat $A/work/arm1004_s2.pid))"
log "END"
