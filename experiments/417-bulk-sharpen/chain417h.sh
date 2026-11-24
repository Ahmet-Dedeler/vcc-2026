#!/usr/bin/env bash
# After 417g's build restores the builder, build 417h = 417e (417a + stack LOAD) + per-context re-level from 417e's artifact
# energy (ctxlevel417.py EV-optimal: A 1.0088 / B 1.0500 / C 1.0736, stable over board raw 0.69-0.70).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain417h.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for 417g's builder restore (pid $$)"
n=0; until grep -qE "builder after restore|stopping|GATE never" "$W/build_g/run417g.log" 2>/dev/null || grep -q "NOT launched" "$W/chain417g.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "417g never restored"; exit 1; }; sleep 30; done
sleep 30
B="fc30e8b0fd884867aa65299ff63b236ceb2f251ea85990a1eeef6a509b6e0eb4"
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder is not fc30e8b0 -- 417h NOT launched"; exit 1; }
[ -s "$REPO/submissions/417e-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload.h5ad.meta.json" ] || { log "no 417e meta -- 417h NOT launched"; exit 1; }
bash "$E/run_loadctx_arm417.sh" 417h 1.25 "$W/dump_h1x2_cd4x15.npz" "A=1.0088,B=1.0500,C=1.0736" \
  "$REPO/submissions/417e-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload.h5ad.meta.json" \
  reviewer-three-is-a-houseplant 417h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-ctxlevel "$B" >> "$LOG" 2>&1
log "417h: $(grep -E '^\[417h\] (PASS|FAIL|CLIFF)' "$W/build_417h/run417h.log" | tr '\n' ' ')"
log "END"
