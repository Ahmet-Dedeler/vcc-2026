#!/usr/bin/env bash
# 387 -- wait for run387.sh to finish cleanly, capture the prior entry's full record, then submit at once (slot open).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/387-sign-repair-arm/work"; RL="$W/build/run387.log"; LOG="$W/submit387.log"
PKG="submissions/387-cap15-partial50-signcons.vcc"; PRIOR="lqWNXX5YaYArdglu6Stn"; NAME="the lint roller has opinions"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for build DONE in $RL"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|GATE never opened|NO ARTIFACT" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|GATE never|NO ARTIFACT' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[387\] PASS" "$RL" || { log "one-key config diff did not PASS -- not submitting"; exit 1; }
grep -q "novelty rc=0" "$RL" || { log "novelty check not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha in build log"; exit 1; }
log "build DONE; package sha $SHA"
uv run python scripts/capture_board_entry.py "$PRIOR" --note "319p2, captured before 387" >> "$LOG" 2>&1 || log "capture of $PRIOR returned rc=$? (may already exist)"
[ -f "data/board-entries/$PRIOR.json" ] || { log "prior capture missing -- not submitting"; exit 1; }
WHEN="$(date -u -v-1M +%Y%m%d%H%M 2>/dev/null || date -u +%Y%m%d%H%M)"
log "ARMING arm_submit.sh when=$WHEN (past => fires now)"
exec scripts/arm_submit.sh "$WHEN" "$PKG" "$SHA" "$NAME" "$LOG"
