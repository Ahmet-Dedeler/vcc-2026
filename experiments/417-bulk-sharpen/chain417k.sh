#!/usr/bin/env bash
# After the window-3 chain ends, build 417k = 417f (417b + stack LOAD, q 1.5) + its own per-context re-level
# (ctxlevel417.py EV-optimal on 417f's artifact energy: A 1.0464 / B 1.0920 / C 1.1140). Read against 417h it is the
# LEVELLED q step (q 1.25 -> 1.5 with the re-level re-derived for q), 420: about +0.004 at k_sh 0.55.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain417k.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for the window-3 chain to end (pid $$)"
n=0; until grep -qE "END$|STOP" "$W/chain_window3.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "window 3 never ended"; exit 1; }; sleep 30; done
sleep 60
B="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder file != HEAD -- 417k NOT launched"; exit 1; }
[ -s "$REPO/submissions/417f-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload.h5ad.meta.json" ] || { log "no 417f meta -- 417k NOT launched"; exit 1; }
log "building 417k on builder ${B:0:8}"
bash "$E/run_loadctx_arm417.sh" 417k 1.5 "$W/dump_h1x2_cd4x15.npz" "A=1.0464,B=1.0920,C=1.1140" \
  "$REPO/submissions/417f-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload.h5ad.meta.json" \
  the-negative-control-went-positive 417k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel "$B" >> "$LOG" 2>&1
log "417k: $(grep -E '^\[417k\] (PASS|FAIL|CLIFF)' "$W/build_417k/run417k.log" | tr '\n' ' ')"
log "END"
