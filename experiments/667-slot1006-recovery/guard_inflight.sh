#!/usr/bin/env bash
# 667 generic in-flight guard (coordinator, 2026-10-06 ~15:35Z), generalising guard1006.sh. Our high-nnz scoring jobs got no machine
# on 10-06 (524dr queued 5.5 h until cancelled; K16 queued since 05:31Z). A queued entry left at the reset blocks the next day's
# 00:01Z submits ("one in flight").
# usage: guard_inflight.sh <AT yyyymmddHHMM> <slot log> [<slot log> ...]
# At AT, ONLY IF limits show in_flight 1 and a `vcc submit` process is running: SIGTERM it (click.echo already flushed
# "entry: <id>" into its arm's capture), read the newest "entry: <id>" from the given logs, `vcc status <id>`, and
# `vcc cancel <id> --yes` only if the status is still "launching" (queued; a running job cannot be cancelled; cancelling a queued
# entry spends no slot). Never submits anything.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
AT="${1:?AT}"; shift; LOGS=("$@"); X=experiments/667-slot1006-recovery; LOG=$X/work/guard_inflight-$AT.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/guard_inflight-$AT.pid
log "armed (pid $$): at $AT cancel our in-flight entry if it is still queued; logs: ${LOGS[*]}"
until [ "$(date -u +%Y%m%d%H%M)" -ge "$AT" ]; do sleep 60; done
INF="$(~/.local/share/uv/tools/vcc-cli/bin/python -c "from vcc import api, auth; print(api.get_limits('https://virtualcellchallenge.org', auth.resolve_token('default').token).get('in_flight'))" 2>&1 | tail -1)"
[ "$INF" = "1" ] || { log "in_flight=$INF -- nothing to cancel"; log "END"; exit 0; }
VP="$(ps -axo pid=,args= | awk '/\/vcc submit / && !/awk/ {print $1; exit}')"
[ -n "$VP" ] || { log "in_flight 1 but no vcc submit process -- stopping (check by hand)"; log "END"; exit 1; }
kill -TERM "$VP"; log "SIGTERMed vcc wait (pid $VP)"
EID=""; for i in $(seq 1 30); do EID="$(cat "${LOGS[@]}" 2>/dev/null | /usr/bin/grep -oE 'entry: [A-Za-z0-9]{20}' | tail -1 | awk '{print $2}')"; [ -n "$EID" ] && break; sleep 5; done
[ -n "$EID" ] || { log "no entry id in the logs after the SIGTERM -- stopping (cancel by hand)"; log "END"; exit 1; }
ST="$(vcc status "$EID" 2>&1 | awk '/status:/ {print $2; exit}')"; log "entry $EID status: ${ST:-?}"
if [ "$ST" = "launching" ]; then OUT="$(vcc cancel "$EID" --yes 2>&1)"; log "vcc cancel rc=$?: $(echo "$OUT" | tr '\n' ' ' | cut -c1-200)"
else log "status is not 'launching' -- not cancelling"; fi
log "END"
