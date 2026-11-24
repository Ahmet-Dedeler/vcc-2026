#!/usr/bin/env bash
# 392 -- wait for run392.sh to pass its gates, then for 387's entry id, then hand off to arm_submit_after.sh
# (which waits until 387 is TERMINAL, captures its full record, re-verifies the hash and submits).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/392-de-sign-prior-arm/work"; RL="$W/build/run392.log"; LOG="$W/submit392.log"
PKG="submissions/392-cap15-partial50-designclean4.vcc"; NAME="the stapler went on sabbatical"
SHA387="986f7665488a54e3e542e71f3e779c26901fd5002e69ec0f7e49454904f7061c"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 392 build DONE"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|GATE never opened|NO ARTIFACT" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|GATE never|NO ARTIFACT' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[392\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
grep -q "novelty rc=0" "$RL" || { log "novelty not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha"; exit 1; }
log "build DONE; package sha $SHA; waiting for 387's entry id in SUBMITTED_LEDGER.json"
n=0; PRIOR="${PRIOR387:-}"   # passed explicitly: the ledger's parsed id can be the scoring JOB token
[ -n "$PRIOR" ] && log "prior entry id given explicitly: $PRIOR"
until [ -n "$PRIOR" ]; do
  PRIOR="$(python3 -c "import json; d=json.load(open('submissions/SUBMITTED_LEDGER.json')); print((d.get('$SHA387') or {}).get('entry_id') or '')" 2>/dev/null)"
  [ -n "$PRIOR" ] && break
  if grep -q "ERROR" "$W/../../387-sign-repair-arm/work/submit387.log" 2>/dev/null; then log "387 submit logged an ERROR -- not arming"; exit 1; fi
  n=$((n+1)); [ $n -ge 480 ] && { log "no 387 entry id after 4 h"; exit 1; }; sleep 30
done
FLOOR="$(date -u +%Y%m%d%H%M)"
log "387 entry $PRIOR; handing off to arm_submit_after.sh floor=$FLOOR"
exec scripts/arm_submit_after.sh "$PRIOR" "$FLOOR" "$PKG" "$SHA" "$NAME" "$LOG"
