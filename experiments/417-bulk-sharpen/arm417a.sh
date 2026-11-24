#!/usr/bin/env bash
# 417a -- the 2026-09-23 FIRST slot (PI decision 20:30Z, docs/NEXT-SLOT-DECISION.md): 410b + bulk-only sharpening q 1.25.
# Waits for run417a_build.sh (config diff PASS, CLIFF OK, novelty rc=0, package sha), then the 00:01Z reset, then submits,
# records the ledger ("entry <id>"-anchored parser) and captures the record.
# FALLBACK: if the build stops, the config diff or the cliff check fails, or work/NOGO-417a exists, this execs the untouched
# 410c waiter (experiments/410-pertscale-revisit/arm410c.sh), which submits 410c after 00:01Z under its own GO-410c.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/417-bulk-sharpen/work"; RL="$W/build_a/run417a.log"; LOG="$W/submit417a.log"
PKG="submissions/417a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125.vcc"; NAME="the pigeon passed peer review"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fallback(){ log "FALLBACK to 410c: $*"; exec bash "$REPO/experiments/410-pertscale-revisit/arm410c.sh"; }
log "WAITING for 417a build DONE (pid $$)"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RL" 2>/dev/null; then fallback "BUILD STOPPED: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RL" | tail -1)"; fi
  n=$((n+1)); [ $n -ge 720 ] && fallback "build not done after 6 h"; sleep 30; done
grep -q "^\[417a\] PASS" "$RL" || fallback "config diff did not PASS"
grep -q "^\[417a\] CLIFF OK" "$RL" || fallback "cliff check not OK: $(grep '\[417a\] CLIFF' "$RL" | tail -1)"
grep -q "novelty rc=0" "$RL" || fallback "novelty not rc=0"
log "build gates passed: $(grep '\[417a\] CLIFF' "$RL" | tail -1)"
n=0; until [ "$(date -u +%Y%m%d%H%M)" -ge 202609230001 ]; do
  [ -f "$W/NOGO-417a" ] && fallback "NOGO-417a present"; n=$((n+1)); [ $n -ge 1440 ] && { log "clock never reached 00:01Z"; exit 1; }; sleep 20; done
[ -f "$W/NOGO-417a" ] && fallback "NOGO-417a present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || fallback "no package sha"
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || fallback "sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting (no fallback: the account state is the problem)"; exit 1; }
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
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "417a bulk-only sharpening q1.25 on 410b: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "417a bulk-only sharpening q1.25 on 410b" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
