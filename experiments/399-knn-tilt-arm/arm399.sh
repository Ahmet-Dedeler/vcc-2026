#!/usr/bin/env bash
# 399 -- wait for run399.sh to pass its gates, then for 395's entry id in the ledger, then hand off to arm_submit_after.sh
# (captures 395's full record, waits until 395 is TERMINAL, re-verifies the hash, submits, records, captures its own entry).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/399-knn-tilt-arm/work"; RL="$W/build/run399.log"; LOG="$W/submit399.log"
PKG="submissions/399-cap15-partial50-knn25o400q80.vcc"; NAME="the fax machine joined a book club"
SHA395="$(grep 'package sha256:' experiments/395-protected-signs-arm/work/build/run395.log | tail -1 | awk '{print $NF}')"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 399 build DONE (pid $$); 395 package sha ${SHA395:-<unknown yet>}"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|GATE never opened|NO ARTIFACT" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|GATE never|NO ARTIFACT' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[399\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
grep -q "novelty rc=0" "$RL" || { log "novelty not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha"; exit 1; }
[ -n "$SHA395" ] || SHA395="$(grep 'package sha256:' experiments/395-protected-signs-arm/work/build/run395.log | tail -1 | awk '{print $NF}')"
log "build DONE; package sha $SHA; waiting for 395's entry id (ledger key $SHA395)"
n=0; PRIOR="${PRIOR395:-}"
until [ -n "$PRIOR" ]; do
  PRIOR="$(python3 -c "import json; d=json.load(open('submissions/SUBMITTED_LEDGER.json')); print((d.get('$SHA395') or {}).get('entry_id') or '')" 2>/dev/null)"
  [ -n "$PRIOR" ] && break
  if grep -q "ERROR" experiments/395-protected-signs-arm/work/submit395.log 2>/dev/null; then log "395 submit logged an ERROR -- not arming"; exit 1; fi
  n=$((n+1)); [ $n -ge 480 ] && { log "no 395 entry id after 4 h"; exit 1; }; sleep 30
done
FLOOR="$(date -u +%Y%m%d%H%M)"
log "395 entry $PRIOR; handing off to arm_submit_after.sh floor=$FLOOR"
exec scripts/arm_submit_after.sh "$PRIOR" "$FLOOR" "$PKG" "$SHA" "$NAME" "$LOG"
