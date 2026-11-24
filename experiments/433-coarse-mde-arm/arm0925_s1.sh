#!/usr/bin/env bash
# 2026-09-25 FIRST slot (rule fixed 01:30Z 09-24 after the 09-24 reads: 417e stack good, 433a MDE transfers):
#   433h = 433a + zero repair -5 + 417h's re-level (paired: repair -> nmae/fid/reach/jac, re-level -> mse/pds).
# Bins: experiments/433-coarse-mde-arm/PREREGISTRATION-0925.md. Waits for the 09-25 00:01Z reset (clock loop capped by a
# timestamp: 12:00Z 09-25), re-verifies the build gates and the package sha, submits, records the ledger, captures.
# work/NOGO-0925-s1 stops it at any point before the submit.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/433-coarse-mde-arm/work"; LOG="$W/submit433h.log"
RL="$REPO/experiments/428-zero-gene-audit/work/build_433h/run433h.log"
PKG="submissions/433h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel.vcc"
NAME="the flow cytometer took up jazz"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fallback(){ log "NOT SUBMITTING: $*"; exit 1; }
log "WAITING for the 09-25 00:01Z reset (pid $$)"
grep -q "^\[433h\] PASS" "$RL" || fallback "config diff did not PASS"
grep -q "^\[433h\] CLIFF OK" "$RL" || fallback "cliff check not OK"
grep -q "novelty rc=0" "$RL" || fallback "novelty not rc=0"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609250001 ]; do
  [ -f "$W/NOGO-0925-s1" ] && fallback "NOGO-0925-s1 present"
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251200 ] && fallback "clock passed 12:00Z 09-25 without reaching the slot"
  sleep 20; done
[ -f "$W/NOGO-0925-s1" ] && fallback "NOGO-0925-s1 present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || fallback "no package sha"
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; [ "$GOT" = "$SHA" ] || fallback "sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || fallback "novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "433h slot 1 on 09-25 (433a + zero repair + re-level): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "433h slot 1 on 09-25" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
