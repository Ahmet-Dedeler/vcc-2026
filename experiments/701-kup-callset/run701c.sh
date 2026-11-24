#!/usr/bin/env bash
# 701 POST HOC H1 size curve, resumed (science [d80a10], 10-16 ~07:35Z): run701b's remaining runs. run701b was stopped at
# 07:15Z because its back-to-back 8 GB runs kept readiness's waiting R-NEW build under its memory floor.
#   nohup bash scripts/detach.sh bash experiments/701-kup-callset/run701c.sh > experiments/701-kup-callset/work/run701c.nohup 2>&1 &
# 1. Waits for R7a's resume waiter (552/work/r7a_after690.sh) to write a NEW " END" and exit (R7a done, or not started).
# 2. Each missing run537.py cpm-n<CPM> base <seed> (CPM 3000 3800 3400 2800, seeds 0 1 2; existing metrics are skipped) goes
#    through gate602 (8 GB), and starts only when no builder chain is alive (bash run_arm_generic.sh / def_day.sh / run690.sh
#    / rabc694.sh / build707.sh, python r7a551.py, matched as interpreter + path). gate602 itself yields to those since the
#    coordinator's 10-16 fix; this also covers the R7a access and build707's whole life.
# 3. reach701.py (condition (i) of the October K rule, fixed in 711's PREREGISTRATION before these runs).
# Release guard: RELEASE-DETECTED before a run -> stop. Cap for starting a run: 10-19 12:00Z. Terminal line, always: "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/701-kup-callset; W=$X/work; LOG=$W/run701c.log
T0="$(date -u +%FT%TZ)"; CAP=202610191200
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
pyalive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^(python[0-9.]*|uv)$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
released(){ [ -e experiments/691-release-watcher/work/RELEASE-DETECTED ]; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
builders(){ local b=""
  for c in run_arm_generic.sh scripts/def_day.sh run690.sh rabc694.sh experiments/707-strong-tier-direction/build707.sh; do
    alive "$c" && b="$b $c"; done
  pyalive experiments/551-r6-c10b/r7a551.py && b="$b r7a551.py"
  echo "$b"; }
echo $$ > "$W/run701c.pid"
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0): waiting for R7a's resume waiter to end"
RA=experiments/552-r7-reserve/work/r7a_after690.log
until newline "$T0" "$RA" " END$" && ! alive experiments/552-r7-reserve/work/r7a_after690.sh; do
  capped && { log "R7a's waiter never ended by $CAP: stopping"; exit 1; }
  sleep 60
done
log "R7a's waiter ended: $(tail -2 "$RA" | head -1 | cut -c1-200)"
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for cpm in 3000 3800 3400 2800; do
  for s in 0 1 2; do
    M="cpm-n$cpm"
    if [ -s "experiments/537-arm-pipeline/work/h1_$M/metrics-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json" ]; then
      log "skip $M seed $s (metrics exist)"; continue; fi
    LAST=""
    while B="$(builders)"; [ -n "$B" ]; do
      [ "$B" = "$LAST" ] || log "waiting for builders:$B"; LAST="$B"
      capped && { log "cap reached while waiting for builders before $M seed $s"; exit 1; }
      sleep 60
    done
    released && { log "skipped: release detected"; exit 0; }
    capped && { log "cap reached before $M seed $s"; exit 1; }
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
