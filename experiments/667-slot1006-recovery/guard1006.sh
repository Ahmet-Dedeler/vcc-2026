#!/usr/bin/env bash
# 667 guard (coordinator, 2026-10-06 ~10:25Z), modelled on 498/guard1005.sh. If 10-06's K16 entry is STILL queued at 23:00Z, it
# would block 10-07's 00:01Z submits ("one in flight"). Our high-nnz scoring jobs have been getting no machine today (524dr's entry
# VqBgraiBxlaK0WP5efOF queued 5.5 h until cancelled; K16 queued since 05:31Z).
# At 23:00Z, ONLY IF K16's arm (pid 71275) is alive AND submit1006s2.log has no capture AND limits show in_flight 1:
#   SIGTERM K16's waiting `vcc submit` (click.echo already flushed "entry: <id>" into the arm's capture), read the id from
#   submit1006s2.log, `vcc status <id>`, and `vcc cancel <id> --yes` only if the status is still "launching" (queued; a running
#   job cannot be cancelled anyway). Cancelling a queued entry spends no slot. Never submits anything.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/667-slot1006-recovery; LOG=$X/work/guard1006.log; S2=experiments/524-october-recipe-dr/work/submit1006s2.log
AT=202610062300; ARM=71275
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/guard1006.pid
log "armed (pid $$): at 23:00Z cancel K16's entry if it is still queued"
until [ "$(date -u +%Y%m%d%H%M)" -ge $AT ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610070001 ] && { log "woke after 10-07 00:01Z -- nothing to do"; log "END"; exit 0; }
  kill -0 $ARM 2>/dev/null || { log "K16's arm exited before 23:00Z -- nothing to do"; log "END"; exit 0; }
  sleep 60
done
kill -0 $ARM 2>/dev/null || { log "K16's arm already exited -- nothing to do"; log "END"; exit 0; }
/usr/bin/grep -q "CAPTURED own entry" $S2 && { log "K16 captured -- nothing to do"; log "END"; exit 0; }
INF="$(~/.local/share/uv/tools/vcc-cli/bin/python -c "from vcc import api, auth; print(api.get_limits('https://virtualcellchallenge.org', auth.resolve_token('default').token).get('in_flight'))" 2>&1 | tail -1)"
[ "$INF" = "1" ] || { log "in_flight=$INF -- nothing to cancel"; log "END"; exit 0; }
VP="$(ps -axo pid=,args= | awk '/vcc submit submissions\/446bj2zk16/ && !/awk/ {print $1; exit}')"
[ -n "$VP" ] || { log "no K16 vcc process found -- stopping"; log "END"; exit 1; }
kill -TERM "$VP"; log "SIGTERMed K16's vcc wait (pid $VP)"
for i in $(seq 1 30); do EID="$(/usr/bin/grep -oE 'entry: [A-Za-z0-9]{20}' $S2 | tail -1 | awk '{print $2}')"; [ -n "$EID" ] && break; sleep 5; done
[ -n "$EID" ] || { log "no entry id in $S2 after the SIGTERM -- stopping (cancel by hand)"; log "END"; exit 1; }
ST="$(vcc status "$EID" 2>&1 | awk '/status:/ {print $2; exit}')"
log "K16 entry $EID status: ${ST:-?}"
if [ "$ST" = "launching" ]; then OUT="$(vcc cancel "$EID" --yes 2>&1)"; log "vcc cancel rc=$?: $(echo "$OUT" | tr '\n' ' ' | cut -c1-200)"
else log "status is not 'launching' -- not cancelling"; fi
log "END"
