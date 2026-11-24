#!/usr/bin/env bash
# Patient follow-up for 10-04 slot 2 (446bj2zk20, entry kJfe2dP9LEJawDEf0pMm), coordinator 2026-10-04.
# arm1004.sh (pid 18260) retries the resume until its 12:00Z cap. Name resolution for storage.googleapis.com keeps failing
# on the iPhone hotspot, so the upload may still be pending then. This waits for that waiter to exit. If no capture
# followed, it keeps resuming THE SAME entry until 10-04 23:00Z: the entry is in 10-04's allowance, so finishing it
# later spends no new slot. It must be done before 10-05's 00:01Z submit (one in flight at a time). It never starts a
# fresh submission: no pending upload and no capture means it stops.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/498-combined-bulk/work; LOG=$W/submit1004s2.log
PKG=submissions/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc
RL=experiments/428-zero-gene-audit/work/build_446bj2zk20/run446bj2zk20.log
NAME="You're behind—very behind."
EID=kJfe2dP9LEJawDEf0pMm; WAITER=18260; CAP=202610042300
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
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
echo $$ > $W/resume1004s2.pid
while kill -0 $WAITER 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge $CAP ] && { log "resume1004s2: cap reached while arm1004 still runs -- stopping"; exit 0; }
  sleep 30
done
if captured; then echo "$(date -u +%FT%TZ) arm1004 captured the entry -- nothing to do" >> $W/resume1004s2.out; exit 0; fi
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$SHA" ] && [ "$GOT" = "$SHA" ] || { log "resume1004s2 refused: sha mismatch got=$GOT want=$SHA"; log "END"; exit 1; }
log "resume1004s2 (coordinator): arm1004 exited without a capture; patient resume of $EID until $CAP"
N=0; RC=1; OUT=""
while [ "$(date -u +%Y%m%d%H%M)" -lt $CAP ]; do
  if ! pending; then log "resume1004s2: $EID is no longer a pending upload and was not captured -- stopping (no fresh submit)"; break; fi
  N=$((N + 1))
  OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$EID" 2>&1)"; RC=$?
  [ $RC -eq 0 ] && break
  log "resume1004s2 attempt $N failed rc=$RC: $(echo "$OUT" | tail -c 160 | tr '\n' ' ')"
  sleep 120
done
log "SUBMIT rc=$RC (resume1004s2, $N attempts)"; printf '%s\n' "$OUT" >> "$LOG"
if [ $RC -eq 0 ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$EID" --note "446bj2zk20 slot 2 on 10-04 (resumed): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$EID" --package "$PKG" --sha256 "$SHA" --note "446bj2zk20 slot 2 on 10-04 (resumed)" \
      --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $EID -> data/board-entries/$EID.json" || log "WARN capture of $EID failed; re-run scripts/capture_board_entry.py $EID"
else
  log "ERROR resume1004s2 gave up at the cap; the pending upload stays for a manual --resume"
fi
log "END"
