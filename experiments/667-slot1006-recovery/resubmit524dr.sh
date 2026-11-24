#!/usr/bin/env bash
# 667 (coordinator, 2026-10-06 ~05:35Z): resubmit 524dr after its first entry was stuck.
# Entry VqBgraiBxlaK0WP5efOF (submitted 00:01:43Z) sat in "launching" for 5.5 h: its Batch scoring job never got a machine, while
# other teams' later entries were scored. The coordinator SIGTERMed the waiting vcc (the id was in its output), then ran
# `vcc cancel --yes` ("stopped the queued scoring job before it started running"; in_flight 0, spent 0). K16 (slot 2) took the
# freed slot at 05:31:20Z. This waits for K16's arm (pid 71275) to exit, then runs chain1006's own slot-1 command for 524dr
# (same package, gate log, -m name). arm1006.sh re-verifies the gates and the package sha and carries the 536b retries.
# 524dr was never scored or ledgered, so the novelty check passes. Cap 10-06 11:30Z (arm1006's own cap is 12:00Z).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/667-slot1006-recovery; LOG=$X/work/resubmit524dr.log; A=experiments/524-october-recipe-dr; K16PID=71275
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/resubmit524dr.pid
log "armed (pid $$): waiting for K16's arm (pid $K16PID) to exit; cap 10-06 11:30Z"
while kill -0 $K16PID 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610061130 ] && { log "cap 11:30Z with K16's arm still running -- not resubmitting 524dr"; log "END"; exit 0; }
  sleep 30
done
log "K16's arm exited; its log tail: $(tail -2 $A/work/submit1006s2.log | tr '\n' ' ' | cut -c1-300)"
nohup bash $A/arm1006.sh 1 524dr $A/work/gates524.log submissions/524dr-october-recipe-abc.vcc \
  "The most dangerous worldviews are the worldviews of those who have never viewed the world." "524dr slot 1 on 10-06 (resubmitted after the stuck entry VqBgraiBxlaK0WP5efOF was cancelled)" \
  > $A/work/arm1006_s1b.nohup 2>&1 & echo $! > $A/work/arm1006_s1b.pid
sleep 3; log "armed arm1006 slot 1 for 524dr (pid $(cat $A/work/arm1006_s1b.pid))"; log "END"
