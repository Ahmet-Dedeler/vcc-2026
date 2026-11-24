#!/usr/bin/env bash
# 417e -- the 2026-09-24 FIRST slot (plan fixed 01:00Z 09-23, docs/NEXT-SLOT-DECISION.md): 417a + stack LOAD.
# s1b (06:25Z 09-23): replaces arm0924_s1.sh, whose clock loop gave up after 8 h (1440 x 20 s), i.e. ~09:00Z 09-23,
# long before the 00:01Z reset. Same gates and log file; the clock loop now waits up to 30 h. Old waiter pid 27780 killed.
# Waits for 417e's build (config diff PASS, CLIFF OK, novelty rc=0, package sha), then the 09-24 00:01Z reset, then submits,
# records the ledger ("entry <id>"-anchored parser) and captures the record.
# If a gate fails or work/NOGO-417e exists, it does NOT submit (no fallback: the slot is the PI's call).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/417-bulk-sharpen/work"; RL="$W/build_417e/run417e.log"; LOG="$W/submit417e.log"
PKG="submissions/417e-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload.vcc"; NAME="the cell line got a lawyer"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fallback(){ log "NOT SUBMITTING: $*"; exit 1; }
log "WAITING for 417e build DONE (pid $$)"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RL" 2>/dev/null; then fallback "BUILD STOPPED: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RL" | tail -1)"; fi
  n=$((n+1)); [ $n -ge 720 ] && fallback "build not done after 6 h"; sleep 30; done
grep -q "^\[417e\] PASS" "$RL" || fallback "config diff did not PASS"
grep -q "^\[417e\] CLIFF OK" "$RL" || fallback "cliff check not OK: $(grep '\[417e\] CLIFF' "$RL" | tail -1)"
grep -q "novelty rc=0" "$RL" || fallback "novelty not rc=0"
log "build gates passed: $(grep '\[417e\] CLIFF' "$RL" | tail -1)"
n=0; until [ "$(date -u +%Y%m%d%H%M)" -ge 202609240001 ]; do
  [ -f "$W/NOGO-417e" ] && fallback "NOGO-417e present"; n=$((n+1)); [ $n -ge 5400 ] && fallback "clock never reached 00:01Z in 30 h"; sleep 20; done
[ -f "$W/NOGO-417e" ] && fallback "NOGO-417e present"
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
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "417e stack LOAD (H1 x2 + CD4 x1.5) on 417a: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "417e stack LOAD on 417a" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
