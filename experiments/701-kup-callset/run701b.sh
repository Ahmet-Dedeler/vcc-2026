#!/usr/bin/env bash
# 701 POST HOC H1 size curve (coordinator [2d6f46], 10-16): the 711 October K rule's condition (i).
#   nohup bash scripts/detach.sh bash experiments/701-kup-callset/run701b.sh > experiments/701-kup-callset/work/run701b.nohup 2>&1 &
# Waits for R7a's wrapper (run by 552/work/r7a_after705.sh) to write a NEW "END rc=" and exit; then, serially, each through
# gate602 (8 GB label; FIFO; it refuses beside a builder): run537.py cpm-n<CPM> base <seed> for CPM 3000 3800 3400 2800
# and seeds 0 1 2 (existing metrics are skipped; the base `cpm` runs exist from 09-30), then reach701.py.
# Release guard: RELEASE-DETECTED before a run -> stop. Cap for starting a run: 10-19 12:00Z.
# Terminal line, always: "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/701-kup-callset; W=$X/work; LOG=$W/run701b.log
T0="$(date -u +%FT%TZ)"; CAP=202610191200
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
released(){ [ -e experiments/691-release-watcher/work/RELEASE-DETECTED ]; }
echo $$ > "$W/run701b.pid"
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0): waiting for R7a's wrapper to end"
RW=experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log
until newline "$T0" "$RW" " END rc=" && ! alive experiments/552-r7-reserve/work/r7a_access_1015b.sh \
      && ! alive experiments/552-r7-reserve/work/r7a_after705.sh; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "R7a never ended by $CAP: stopping"; exit 1; }
  sleep 60
done
log "R7a ended: $(tail -1 "$RW" | cut -c1-160)"
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for cpm in 3000 3800 3400 2800; do
  for s in 0 1 2; do
    M="cpm-n$cpm"
    if [ -s "experiments/537-arm-pipeline/work/h1_$M/metrics-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json" ]; then
      log "skip $M seed $s (metrics exist)"; continue; fi
    released && { log "skipped: release detected"; exit 0; }
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "cap reached before $M seed $s"; exit 1; }
    log "START $M seed $s"
    GATE602_QUIET=1 bash experiments/602-council-gate/gate602.sh "701-$M-s$s" 8 -- uv run python experiments/537-arm-pipeline/run537.py "$M" base "$s" $ARGS \
      > "$W/run_${M}_s$s.log" 2>&1
    rc=$?
    log "END rc=$rc $M seed $s: $(tail -1 "$W/run_${M}_s$s.log" | cut -c1-160)"
    [ "$rc" = 0 ] || { log "run failed: stopping"; exit 1; }
  done
done
log "reading (POST HOC)"
uv run python $X/reach701.py > "$W/reach701.out" 2>&1
log "reach701 rc=$?: $(/usr/bin/grep -E '^f |CONDITION' "$W/reach701.out" | tr '\n' ' ' | cut -c1-900)"
