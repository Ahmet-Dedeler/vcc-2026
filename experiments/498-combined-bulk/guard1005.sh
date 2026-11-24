#!/usr/bin/env bash
# One-shot guard (coordinator, 2026-10-04): if 10-04 slot 2's K20 upload (entry kJfe2dP9LEJawDEf0pMm) is STILL pending at
# 10-04 23:30Z, cancel it with `vcc cancel` so 10-05's 00:01Z submit (496a, the hook read that decides mde_sign) is not
# blocked by "one in flight at a time". K20's package stays novel (it was never recorded as submitted), so it can be
# re-submitted on a later day. Does nothing if the entry was captured or is no longer pending.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/498-combined-bulk/work; LOG=$W/submit1004s2.log; EID=kJfe2dP9LEJawDEf0pMm; AT=202610042330
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/guard1005.pid
until [ "$(date -u +%Y%m%d%H%M)" -ge $AT ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610050001 ] && exit 0   # woke too late: 10-05's own waiter is already deciding
  sleep 60
done
captured(){ awk '/ARMED slot 2:/ {b = ""} {b = b $0 "\n"} END {printf "%s", b}' "$LOG" | /usr/bin/grep -q "CAPTURED own entry"; }
pending(){ python3 - "$EID" <<'PY'
import json, os, sys
try:
    d = json.load(open(os.path.expanduser("~/.config/vcc/state.json")))
except Exception:
    sys.exit(1)
sys.exit(0 if sys.argv[1] in d.get("profiles", {}).get("default", {}).get("pending_uploads", {}) else 1)
PY
}
if captured || ! pending; then echo "$(date -u +%FT%TZ) guard1005: K20 captured or not pending -- nothing to do" >> $W/guard1005.out; exit 0; fi
log "guard1005: K20 entry $EID still pending at 23:30Z -- cancelling it to free the in-flight slot for 10-05"
OUT="$(vcc cancel "$EID" 2>&1)"; RC=$?
log "guard1005: vcc cancel rc=$RC: $(echo "$OUT" | tail -c 200 | tr '\n' ' ')"
log "END"
