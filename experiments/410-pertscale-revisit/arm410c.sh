#!/usr/bin/env bash
# 410c -- wait for run410c_build.sh to pass its gates (config diff PASS, novelty rc=0, package sha), then submit (2026-09-22 has
# both slots unspent and nothing in flight), record the ledger with the "entry <id>"-anchored parser, capture the record.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/410-pertscale-revisit/work"; RL="$W/build_c/run410c.log"; LOG="$W/submit410c.log"
PKG="submissions/410c-cap15-partial50-clean4q80-pertscaleinf-ballast.vcc"; NAME="the kettle applied for tenure"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 410 build DONE (pid $$)"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[410c\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
# gates added for 410b: 410 must be terminal (captured) and the PI must have written work/GO-410b after reading H1's ballast rows
n=0; until grep -q "CAPTURED own entry" "$W/submit410b.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 480 ] && { log "410 never captured -- not submitting"; exit 1; }; sleep 30; done
log "410b is captured; waiting for work/GO-410c (or NOGO-410c) and the 2026-09-23 00:00Z reset"
n=0; until [ -f "$W/GO-410c" ] && [ "$(date -u +%Y%m%d%H%M)" -ge 202609230001 ]; do [ -f "$W/NOGO-410c" ] && { log "NOGO-410c present -- not submitting"; exit 1; }; n=$((n+1)); [ $n -ge 1440 ] && { log "no GO after 8 h"; exit 1; }; sleep 20; done
log "GO-410c: $(cat "$W/GO-410c")"
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
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "410c uncapped per-target bulk key through ballast: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "410c uncapped per-target key through ballast" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
