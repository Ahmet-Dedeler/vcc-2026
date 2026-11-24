#!/usr/bin/env bash
# 402 -- wait for run402.sh to pass its gates, then for 399's entry id in the ledger, then hand off to arm_submit_after.sh
# (captures 399's full record, waits until 399 is TERMINAL, re-verifies the hash, submits, records, captures its own entry).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/402-ctxA-jurkat-knn/work"; RL="$W/build/run402.log"; LOG="$W/submit402.log"
PKG="submissions/402-cap15-partial50-knn25o400jaq80.vcc"; NAME="the paperclip demanded a raise"
SHA399="$(grep 'package sha256:' experiments/399-knn-tilt-arm/work/build/run399.log | tail -1 | awk '{print $NF}')"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 399 build DONE (pid $$); 399 package sha ${SHA399:-<unknown yet>}"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|GATE never opened|NO ARTIFACT" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|GATE never|NO ARTIFACT' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[402\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
grep -q "novelty rc=0" "$RL" || { log "novelty not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha"; exit 1; }
[ -n "$SHA399" ] || SHA399="$(grep 'package sha256:' experiments/399-knn-tilt-arm/work/build/run399.log | tail -1 | awk '{print $NF}')"
log "build DONE; package sha $SHA; waiting for 399's entry id (ledger key $SHA399)"
n=0; PRIOR="${PRIOR399:-}"
until [ -n "$PRIOR" ]; do
  PRIOR="$(python3 -c "import json; d=json.load(open('submissions/SUBMITTED_LEDGER.json')); print((d.get('$SHA399') or {}).get('entry_id') or '')" 2>/dev/null)"
  [ -n "$PRIOR" ] && break
  if grep -q "ERROR" experiments/399-knn-tilt-arm/work/submit399.log 2>/dev/null; then log "399 submit logged an ERROR -- not arming"; exit 1; fi
  n=$((n+1)); [ $n -ge 480 ] && { log "no 399 entry id after 4 h"; exit 1; }; sleep 30
done
FLOOR="$(date -u +%Y%m%d%H%M)"
log "399 entry $PRIOR; handing off to arm_submit_after.sh floor=$FLOOR"
exec scripts/arm_submit_after.sh "$PRIOR" "$FLOOR" "$PKG" "$SHA" "$NAME" "$LOG"
