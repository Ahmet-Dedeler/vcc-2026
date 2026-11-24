#!/usr/bin/env bash
# arm1015f.sh (690f, [17adab] 2026-10-14 for review 720 F1): arm1015.sh + a FINAL_REUPLOAD=1 branch. arm1015.sh itself is
# untouched (a live waiter executes it). Without FINAL_REUPLOAD the behaviour is arm1015.sh's, byte for byte below that line.
#   FINAL_REUPLOAD=1: the final phase's "re-upload the best scored package unchanged" (DEF-RUNBOOK §6, FINAL-PHASE rule 2).
#     Taken only if ALL hold, else NOT SUBMITTING: the clock is past FINALPHASE_FROM (default 202610220700); the name is
#     the reserved final name from docs/MODEL-NAMES.md, or FINAL_REUPLOAD_TEST=1 (the one early harmless test of the
#     server's behaviour on identical bytes, never first on Nov 4); the novelty check REFUSES the package (it IS ledgered)
#     and the ledgered entry's capture (data/board-entries/<id>.json) is a FINAL-phase entry (is_final or partition != val).
#     After the submit it records the re-upload under its own ledger key "<sha>#reupload-<new id>" with reupload_of, so
#     the original record is kept. Every decision is logged.
# 670 (coordinator [400b6e], 2026-10-14): a generic slot waiter for 10-15 onward, built from 660/arm1007.sh (gates, 536b retry +
# resume block, novelty record, capture) with the readiness and prior-entry waits made explicit.
#   bash arm1015.sh <slot-label> <TAG> <ready-log> <run-log> <package.vcc> "<joke name>" "<ledger note>" <prior-log> <FLOOR> <CAP>
#     ready-log : a build chain's log; the waiter needs a line "READY <TAG> <sha>" in it (a "NOT READY" line stops it)
#     run-log   : run_arm_generic.sh's run log for TAG (gates: "[TAG] PASS", "[TAG] CLIFF OK", "novelty rc=0", package sha)
#     prior-log : the log of the previous slot's waiter/arm. We wait until it has a terminal line (" END", "SUBMITTED ok" with
#                 its own capture, "NOT SUBMITTING", "ERROR"), then capture its entry id if it has one (never overwrite raw members).
#     FLOOR/CAP : UTC yyyymmddHHMM; never submit before FLOOR or after CAP.
#     [READYTAG]: the tag in the ready-log's READY line when it differs from the run log's TAG (dn671.sh prints "READY dn671").
# work/NOGO-<slot-label> stops it before the submit. -d stays EMPTY (AGENTS §14).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SL="${1:?slot label}"; TAG="${2:?tag}"; RDY="${3:?ready log}"; RL="${4:?run log}"; PKG="${5:?package}"; NAME="${6:?name}"; NOTE="${7:?note}"
PRIOR="${8:?prior log}"; FLOOR="${9:?floor}"; CAP="${10:?cap}"; RTAG="${11:-$TAG}"
W="$REPO/experiments/670-coord-1015/work"; mkdir -p "$W"; LOG="$W/submit-$SL.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
stop(){ log "NOT SUBMITTING: $*"; log "END"; exit 1; }
now(){ date -u +%Y%m%d%H%M; }
echo $$ > "$W/submit-$SL.pid"
log "ARMED $SL: $TAG pkg=$PKG name='$NAME' floor=$FLOOR cap=$CAP prior=$PRIOR (pid $$)"
# 1. the build must print READY
until /usr/bin/grep -qE "READY $RTAG [0-9a-f]{64}" "$RDY" 2>/dev/null; do
  /usr/bin/grep -q "NOT READY" "$RDY" 2>/dev/null && stop "build log says: $(/usr/bin/grep 'NOT READY' "$RDY" | tail -1)"
  [ -f "$W/NOGO-$SL" ] && stop "NOGO-$SL present"
  [ "$(now)" -ge "$CAP" ] && stop "no READY by CAP"
  sleep 30; done
SHA="$(/usr/bin/grep -oE "READY $RTAG [0-9a-f]{64}" "$RDY" | tail -1 | awk '{print $3}')"
log "build READY: $SHA"
/usr/bin/grep -q "^\[$TAG\] PASS" "$RL" || stop "config diff did not PASS in $RL"
/usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" || stop "cliff check not OK"
/usr/bin/grep -q "novelty rc=0" "$RL" || stop "novelty not rc=0"
# 2. the prior slot must be terminal; capture its entry before ours can replace it on the board API
until /usr/bin/grep -qE " END$|SUBMIT rc=|NOT SUBMITTING|ERROR " "$PRIOR" 2>/dev/null; do
  [ -f "$W/NOGO-$SL" ] && stop "NOGO-$SL present"
  [ "$(now)" -ge "$CAP" ] && stop "prior never terminal by CAP"
  sleep 30; done
PID_PRIOR="$(/usr/bin/grep -oE 'pid=[0-9]+|pid [0-9]+' "$PRIOR" | head -1 | /usr/bin/grep -oE '[0-9]+')"
while [ -n "$PID_PRIOR" ] && kill -0 "$PID_PRIOR" 2>/dev/null; do [ "$(now)" -ge "$CAP" ] && stop "prior pid $PID_PRIOR alive at CAP"; sleep 30; done
PID_ID="$(/usr/bin/grep -oE 'entry[: ]+[A-Za-z0-9]{20}|entry id: [A-Za-z0-9]{20}' "$PRIOR" | tail -1 | /usr/bin/grep -oE '[A-Za-z0-9]{20}$')"
if [ -n "$PID_ID" ]; then
  if [ -s "data/board-entries/$PID_ID.json" ] || uv run python scripts/capture_board_entry.py "$PID_ID" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1; then
    log "prior entry $PID_ID captured"
  else stop "could not capture prior entry $PID_ID (would erase its raw members)"; fi
else log "prior log has no entry id (prior did not submit): going alone"; fi
# 3. floor and cap
until [ "$(now)" -ge "$FLOOR" ]; do [ -f "$W/NOGO-$SL" ] && stop "NOGO-$SL present"; sleep 20; done
[ "$(now)" -ge "$CAP" ] && stop "clock past CAP $CAP"
[ -f "$W/NOGO-$SL" ] && stop "NOGO-$SL present"
# 4. the limits endpoint must show nothing in flight (a stuck entry blocks us; cancel it by hand, see 667)
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; [ "$GOT" = "$SHA" ] || stop "sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | /usr/bin/grep -q '"can_submit": true' || stop "account cannot submit"
BE="${BOARD_ENTRIES691:-data/board-entries}"; REUP_OF=""; LEDGER="${LEDGER691:-submissions/SUBMITTED_LEDGER.json}"   # *691: tests only
if [ "${FINAL_REUPLOAD:-0}" = 1 ]; then
  [ "$(now)" -ge "${FINALPHASE_FROM:-202610220700}" ] || stop "FINAL_REUPLOAD before the final phase (${FINALPHASE_FROM:-202610220700})"
  RES="$(python3 -c 'import re,sys; m=re.search(r"\*\*Final Nov 5 entry:\*\* `([^`]+)`", open(sys.argv[1]).read()); print(m.group(1) if m else "")' docs/MODEL-NAMES.md)"
  if [ "$NAME" != "$RES" ]; then
    [ "${FINAL_REUPLOAD_TEST:-0}" = 1 ] || stop "FINAL_REUPLOAD needs the reserved final name '$RES' (or FINAL_REUPLOAD_TEST=1 for the early test)"
    log "FINAL_REUPLOAD_TEST=1: the early identical-bytes test under a queue name"
  fi
  NOV="$(uv run python scripts/check_submission_novelty.py --ledger "$LEDGER" check "$PKG" 2>&1)"; NRC=$?
  printf '%s\n' "$NOV" >> "$LOG"
  [ "$NRC" = 1 ] && echo "$NOV" | /usr/bin/grep -q "ALREADY BEEN SCORED" || stop "FINAL_REUPLOAD: the package is not a ledgered one (novelty rc $NRC): use the normal path"
  REUP_OF="$(echo "$NOV" | /usr/bin/grep -oE 'entry +[A-Za-z0-9]{20}' | head -1 | awk '{print $2}')"
  [ -n "$REUP_OF" ] || stop "FINAL_REUPLOAD: no entry id in the novelty refusal"
  python3 -c 'import json,sys; e=json.load(open(sys.argv[1]))["entry"]; sys.exit(0 if (e.get("is_final") is True or e.get("partition") not in (None, "", "val")) else 1)' "$BE/$REUP_OF.json" 2>/dev/null \
    || stop "FINAL_REUPLOAD: ledgered entry $REUP_OF has no FINAL-phase capture in $BE (an A/B/C package is never a final re-upload)"
  log "FINAL_REUPLOAD: re-uploading the scored final package of entry $REUP_OF, unchanged, as '$NAME'"
else
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || stop "novelty check refused"
fi
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
_pend536(){ python3 - "$PKG" <<'PY'
import json, os, sys
try:
    d = json.load(open(os.path.expanduser("~/.config/vcc/state.json")))
except Exception:
    sys.exit(0)
for eid, r in d.get("profiles", {}).get("default", {}).get("pending_uploads", {}).items():
    if os.path.realpath(str(r.get("local_path", ""))) == os.path.realpath(sys.argv[1]):
        print(eid)
PY
}
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; N536=1
CAPS536=$(date -j -u -f %Y%m%d%H%M "$CAP" +%s 2>/dev/null); [ -n "$CAPS536" ] || CAPS536=$(( $(date +%s) + 21600 ))
B536=${VCC_RETRY_BASE_S:-90}
MAXTRY=60; [ "${FINAL_REUPLOAD_TEST:-0}" = 1 ] && MAXTRY=2   # 720 R6: the early test reads a deterministic refusal in minutes
while [ $RC -ne 0 ] && [ $N536 -lt $MAXTRY ] && [ "$(now)" -lt "$CAP" ]; do
  S536=$(( B536 << (N536 < 5 ? N536 - 1 : 4) )); L536=$(( CAPS536 - $(date +%s) - 60 ))
  if [ $L536 -le 0 ]; then log "submit attempt $N536 failed rc=$RC; no time left before the cap: stopping"; break; fi
  [ $S536 -gt $L536 ] && S536=$L536
  P536="$(_pend536 | head -1)"
  if [ -n "$P536" ]; then
    log "submit attempt $N536 failed rc=$RC; resuming the pending upload $P536 in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$P536" 2>&1)"; RC=$?
  elif echo "$OUT" | /usr/bin/grep -qi "entry\|scor"; then
    log "submit attempt $N536 failed rc=$RC after an entry existed and no upload is pending: not retrying: $(echo "$OUT" | tail -c 200)"; break
  else
    log "submit attempt $N536 failed rc=$RC before any entry; fresh retry in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
  fi
  N536=$((N536 + 1)); done
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | /usr/bin/grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | /usr/bin/grep -o '\b[A-Za-z0-9]\{20\}\b' | /usr/bin/grep -v -E '^[0-9]+$' | /usr/bin/grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  if [ -n "$REUP_OF" ]; then
    python3 - "$PKG" "$ID" "$REUP_OF" "$NOTE: $NAME" "$LEDGER" <<'PY' >> "$LOG" 2>&1 && log "ledger: re-upload of $REUP_OF recorded as $ID" || log "WARN ledger record of the re-upload failed"
import datetime, hashlib, json, sys
pkg, eid, of, note, L = sys.argv[1:6]
h = hashlib.sha256()
with open(pkg, "rb") as f:
    for b in iter(lambda: f.read(1 << 24), b""):
        h.update(b)
d = h.hexdigest()
led = json.load(open(L))
led[f"{d}#reupload-{eid}"] = {"package": pkg, "entry_id": eid, "reupload_of": of, "note": "FINAL RE-UPLOAD: " + note,
                              "date": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
json.dump(led, open(L, "w"), indent=2, sort_keys=True); open(L, "a").write("\n")
print(f"recorded {d[:16]}#reupload-{eid} (reupload_of {of})")
PY
  else
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "$NOTE: $NAME" >> "$LOG" 2>&1
  fi
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "$NOTE" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
