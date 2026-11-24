#!/usr/bin/env bash
# 395 -- wait for run395.sh to pass its gates, then submit immediately (no entry in flight; both 2026-09-16 slots are
# unspent), record the ledger with the entry id parsed from the "entry <id>" line (the 20-char grep in arm_submit.sh
# matches the lowercase scoring-job token first), and capture the full record.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
source .venv/bin/activate 2>/dev/null || true
W="$REPO/experiments/395-protected-signs-arm/work"; RL="$W/build/run395.log"; LOG="$W/submit395.log"
PKG="submissions/395-cap15-partial50-designclean4q80.vcc"; NAME="the thermostat has a union rep"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 395 build DONE (pid $$)"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do   # (b) the earlier over-strict stop line is expected; only a stop AFTER FINISH395 counts
  if sed -n "/FINISH395/,$p" "$RL" | grep -qE "stopping|PREP FAILED"; then log "FINISH STOPPED: $(sed -n "/FINISH395/,$p" "$RL" | grep -E "stopping|PREP FAILED" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[395\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
grep -q "novelty rc=0" "$RL" || { log "novelty not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha"; exit 1; }
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { log "ERROR sha mismatch got=$GOT want=$SHA -- NOT submitting"; exit 1; }
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || { log "ERROR novelty check refused -- NOT submitting"; exit 1; }
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "395 clean4q80 on 392: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "395 q80-protected block signs on 392" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
