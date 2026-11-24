#!/usr/bin/env bash
# After 417c's build restores the builder (and made an artifact), pin 417d to the same window-2 builder and launch it.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain417d.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for 417c's builder restore (pid $$)"
n=0; until grep -qE "builder after restore|stopping|GATE never" "$W/build_c/run417c.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "417c never restored"; exit 1; }; sleep 30; done
grep -q "builder after restore" "$W/build_c/run417c.log" || { log "417c stopped before restore -- 417d NOT launched"; exit 1; }
sleep 5; grep -q "NO ARTIFACT" "$W/build_c/run417c.log" && { log "417c made no artifact -- 417d NOT launched"; exit 1; }
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
WANT="$(grep -o 'B" = "[0-9a-f]\{64\}' "$E/run417c_build.sh" | head -1 | cut -d'"' -f3)"
[ -n "$WANT" ] && [ "$B" = "$WANT" ] || { log "builder $B != 417c's pinned $WANT -- 417d NOT launched"; exit 1; }
[ -z "$(git status --short scripts/make_submission_pm_cd4.py)" ] || { log "builder dirty -- 417d NOT launched"; exit 1; }
sed -i '' "s/__BUILDER_SHA__/$B/g" "$E/run417d_build.sh"
git add "$E/run417d_build.sh" && git commit -q -m "417d build pinned to builder ${B:0:8}

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
nohup bash "$E/run417d_build.sh" > "$W/run417d.nohup" 2>&1 &
log "417d build launched pid $! on builder ${B:0:8}"
