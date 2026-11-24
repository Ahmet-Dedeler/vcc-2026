#!/usr/bin/env bash
# After 417f's build restores the builder, build 417g (417a + per-context re-level). Builder must be fc30e8b0.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain417g.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for 417f's builder restore (pid $$)"
n=0; until grep -qE "builder after restore|stopping|GATE never" "$W/build_417f/run417f.log" 2>/dev/null || grep -qE "NOT built|END$" "$W/chain_loadarms417.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "417f never restored"; exit 1; }; sleep 30; done
sleep 30
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "fc30e8b0fd884867aa65299ff63b236ceb2f251ea85990a1eeef6a509b6e0eb4" ] || { log "builder is not fc30e8b0 -- 417g NOT launched"; exit 1; }
nohup bash "$E/run417g_build.sh" > "$W/run417g.nohup" 2>&1 &
log "417g build launched pid $!"
