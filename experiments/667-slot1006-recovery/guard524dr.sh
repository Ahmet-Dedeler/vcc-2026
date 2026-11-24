#!/usr/bin/env bash
# 667: 23:30Z guard for 524dr's resumed entry KTzIDR0Occ2jYzklkARE (coordinator, 2026-10-06 ~22:15Z). If its status is still
# "launching" (queued, no machine) at 23:30Z, it stops resume524dr.sh (so it cannot retry), SIGTERMs the vcc wait, and cancels the
# entry (`vcc cancel --yes` spends no slot), so CPH's 10-07 00:01Z submit is not blocked. Anything else: nothing to do.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/667-slot1006-recovery; LOG=$X/work/guard524dr.log; EID=KTzIDR0Occ2jYzklkARE
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/guard524dr.pid
log "armed (pid $$): at 23:30Z cancel $EID if still launching"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610062330 ]; do [ "$(date -u +%Y%m%d%H%M)" -ge 202610070001 ] && { log "woke after 00:01Z"; log "END"; exit 0; }; sleep 30; done
ST="$(vcc status "$EID" 2>&1 | awk '/status:/ {print $2; exit}')"; log "$EID status: ${ST:-?}"
if [ "$ST" = "launching" ] || [ "$ST" = "uploading" ]; then
  RP="$(cat $X/work/resume524dr.pid 2>/dev/null)"; [ -n "$RP" ] && ps -o args= -p "$RP" | /usr/bin/grep -q resume524dr && kill "$RP" && log "stopped resume524dr (pid $RP)"
  VP="$(ps -axo pid=,args= | awk '/vcc submit submissions\/524dr-october-recipe-abc/ && !/awk/ {print $1; exit}')"; [ -n "$VP" ] && kill -TERM "$VP" && log "SIGTERMed vcc (pid $VP)"
  OUT="$(vcc cancel "$EID" --yes 2>&1)"; log "vcc cancel rc=$?: $(echo "$OUT" | tr '\n' ' ' | cut -c1-200)"
else log "not launching/uploading -- nothing to do"; fi
log "END"
