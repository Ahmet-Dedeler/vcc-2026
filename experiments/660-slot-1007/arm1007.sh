#!/usr/bin/env bash
# 2026-10-07 slot waiter (one script, two roles), launched by 660/chain1007.sh at 10-06 20:00Z (653's v4 rule, fixed 10-02)
# (docs/NEXT-SLOT-DECISION.md, "the 10-07 / 10-08 rule, FIXED before any 10-06 entry exists").
#   slot 1: bash arm1007.sh 1 <tag> <build run log> <package.vcc> "<joke name>" "<ledger note>"
#           waits for the 10-07 00:01Z reset (cap 12:00Z 10-07), re-verifies the gates and the package sha, submits
#   slot 2: bash arm1007.sh 2 <tag> <build run log> <package.vcc> "<joke name>" "<ledger note>"
#           waits for slot 1's "CAPTURED own entry" in work/submit1007s1.log, then the same
# work/NOGO-1007-s<n> stops either before its submit. -d is empty; -m is the joke name (AGENTS §14).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
S="${1:?slot}"; TAG="${2:?tag}"; RL="${3:?run log}"; PKG="${4:?package}"; NAME="${5:?name}"; NOTE="${6:?note}"
W="$REPO/experiments/660-slot-1007/work"; LOG="$W/submit1007s$S.log"; S1="$W/submit1007s1.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
stop(){ log "NOT SUBMITTING: $*"; exit 1; }
log "ARMED slot $S: $TAG pkg=$PKG name='$NAME' (pid $$)"
grep -q "^\[$TAG\] PASS" "$RL" || stop "config diff did not PASS in $RL"
grep -q "^\[$TAG\] CLIFF OK" "$RL" || stop "cliff check not OK"
grep -q "novelty rc=0" "$RL" || stop "novelty not rc=0"
if [ "$S" = "1" ]; then
  until [ "$(date -u +%Y%m%d%H%M)" -ge 202610070001 ]; do
    [ -f "$W/NOGO-1007-s1" ] && stop "NOGO-1007-s1 present"
    [ "$(date -u +%Y%m%d%H%M)" -ge 202610071200 ] && stop "clock passed 12:00Z 10-07 without reaching the slot"
    sleep 20; done
else
  # slot 2 goes after slot 1's entry is captured, or after slot 1 has finished without one (NOT SUBMITTING, or its END line):
  # each arm here is read against 446bj2 on its own, so a failed slot 1 must not waste slot 2 (2026-09-27)
  # read only slot 1's current run: the lines from its LAST "ARMED slot 1:" line on (an earlier, killed run's lines in the
  # same log do not count; a slot 1 that fails its gates within milliseconds of arming still counts)
  new1(){ awk '/ARMED slot 1:/ {b = ""} {b = b $0 "\n"} END {printf "%s", b}' "$S1" 2>/dev/null; }
  until new1 | grep -q "CAPTURED own entry"; do
    if new1 | grep -qE "NOT SUBMITTING| END$"; then
      until [ "$(date -u +%Y%m%d%H%M)" -ge 202610070001 ]; do sleep 20; done
      log "slot 1 finished without a captured entry -- slot 2 goes alone"; break
    fi
    [ -f "$W/NOGO-1007-s2" ] && stop "NOGO-1007-s2 present"
    [ "$(date -u +%Y%m%d%H%M)" -ge 202610071200 ] && stop "slot 1 never captured by 12:00Z 10-07"
    sleep 30; done
fi
# 501 section 1 (2026-09-27): slot 1's cap sits inside a loop that stops running at the reset, so it never fired on a
# late wake. Check the clock again: never submit this day's package after 12:00Z 10-07 (or on a later day).
[ "$(date -u +%Y%m%d%H%M)" -ge 202610071200 ] && stop "clock is past 12:00Z 10-07: not submitting this day's package late"
[ -f "$W/NOGO-1007-s$S" ] && stop "NOGO-1007-s$S present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || stop "no package sha in $RL"
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; [ "$GOT" = "$SHA" ] || stop "sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || stop "novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
# 536 (2026-09-30): retry loop -- resume an interrupted upload of THIS package, retry fresh only before any entry existed,
# stop if an entry exists without a pending upload; up to 120 attempts (network-patient since 10-04: overnight hotspot drops killed both 10-04 uploads after 6 quick tries), never past this day's cap (202610071200).
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
# 536b (2026-10-04): back off 90 s -> 24 min and keep trying until this day's cap; never sleep across the cap
CAPS536=$(date -j -u -f %Y%m%d%H%M 202610071200 +%s 2>/dev/null); [ -n "$CAPS536" ] || CAPS536=$(( $(date +%s) + 21600 ))
B536=${VCC_RETRY_BASE_S:-90}
while [ $RC -ne 0 ] && [ $N536 -lt 60 ] && [ "$(date -u +%Y%m%d%H%M)" -lt 202610071200 ]; do
  S536=$(( B536 << (N536 < 5 ? N536 - 1 : 4) )); L536=$(( CAPS536 - $(date +%s) - 60 ))
  if [ $L536 -le 0 ]; then log "submit attempt $N536 failed rc=$RC; no time left before the cap: stopping"; break; fi
  [ $S536 -gt $L536 ] && S536=$L536
  P536="$(_pend536 | head -1)"
  if [ -n "$P536" ]; then
    log "submit attempt $N536 failed rc=$RC; resuming the pending upload $P536 in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$P536" 2>&1)"; RC=$?
  elif echo "$OUT" | grep -qi "entry\|scor"; then
    log "submit attempt $N536 failed rc=$RC after an entry existed and no upload is pending: not retrying (a fresh submit could double-submit): $(echo "$OUT" | tail -c 200)"; break
  else
    log "submit attempt $N536 failed rc=$RC before any entry; fresh retry in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
  fi
  N536=$((N536 + 1)); done
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "$NOTE: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "$NOTE" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
