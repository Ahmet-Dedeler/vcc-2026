#!/usr/bin/env bash
# 2026-09-28 slot waiter (one script, two roles), launched by 456/chain0928.sh after the 09-27 reads fix the branch
# (docs/NEXT-SLOT-DECISION.md, "the 09-28 rule, FIXED before any 09-27 entry exists").
#   slot 1: bash arm0928.sh 1 <tag> <build run log> <package.vcc> "<joke name>" "<ledger note>"
#           waits for the 09-28 00:01Z reset (cap 12:00Z 09-28), re-verifies the gates and the package sha, submits
#   slot 2: bash arm0928.sh 2 <tag> <build run log> <package.vcc> "<joke name>" "<ledger note>"
#           waits for slot 1's "CAPTURED own entry" in work/submit0928s1.log, then the same
# work/NOGO-0928-s<n> stops either before its submit. -d is empty; -m is the joke name (AGENTS §14).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
S="${1:?slot}"; TAG="${2:?tag}"; RL="${3:?run log}"; PKG="${4:?package}"; NAME="${5:?name}"; NOTE="${6:?note}"
W="$REPO/experiments/456-zs-riders/work"; LOG="$W/submit0928s$S.log"; S1="$W/submit0928s1.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
stop(){ log "NOT SUBMITTING: $*"; exit 1; }
log "ARMED slot $S: $TAG pkg=$PKG name='$NAME' (pid $$)"
grep -q "^\[$TAG\] PASS" "$RL" || stop "config diff did not PASS in $RL"
grep -q "^\[$TAG\] CLIFF OK" "$RL" || stop "cliff check not OK"
grep -q "novelty rc=0" "$RL" || stop "novelty not rc=0"
if [ "$S" = "1" ]; then
  until [ "$(date -u +%Y%m%d%H%M)" -ge 202609280001 ]; do
    [ -f "$W/NOGO-0928-s1" ] && stop "NOGO-0928-s1 present"
    [ "$(date -u +%Y%m%d%H%M)" -ge 202609281200 ] && stop "clock passed 12:00Z 09-28 without reaching the slot"
    sleep 20; done
else
  until grep -q "CAPTURED own entry" "$S1" 2>/dev/null; do
    grep -q "NOT SUBMITTING" "$S1" 2>/dev/null && { log "slot 1 did not submit -- slot 2 is the PI's call; exiting"; exit 0; }
    [ -f "$W/NOGO-0928-s2" ] && stop "NOGO-0928-s2 present"
    [ "$(date -u +%Y%m%d%H%M)" -ge 202609282330 ] && stop "slot 1 never captured by 23:30Z 09-28 (cap moved from 12:00Z after the 09-27 20:23Z reboot)"
    sleep 30; done
fi
[ -f "$W/NOGO-0928-s$S" ] && stop "NOGO-0928-s$S present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || stop "no package sha in $RL"
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; [ "$GOT" = "$SHA" ] || stop "sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || stop "novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
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
