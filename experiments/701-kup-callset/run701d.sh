#!/usr/bin/env bash
# 701 EXPLORATORY H1 point at CPM 4600 = K 5,000 (science [d80a10], 10-16 ~09:25Z), for the coordinator's 10-19/20 slot ask:
# if 711 (K 4,200) lands at a bin edge, a K 5,000 dose partner needs its held-out price first (11.1c.3). NOT part of the
# October K rule's condition (i), which stays reach701.py's 5 sizes. Reported by reach701d.py, labelled exploratory.
#   nohup bash scripts/detach.sh bash experiments/701-kup-callset/run701d.sh > experiments/701-kup-callset/work/run701d.nohup 2>&1 &
# 1. Waits for run701c.sh to write a NEW " END" and exit.
# 2. run537.py cpm-n4600 base <seed> for seeds 0 1 2 (existing metrics skipped) through gate602 (8 GB), each only when no
#    builder chain is alive (same list as run701c).
# 3. reach701d.py. Release guard; cap for starting a run 10-19 12:00Z. Terminal line, always: "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/701-kup-callset; W=$X/work; LOG=$W/run701d.log
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
echo $$ > "$W/run701d.pid"
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0): waiting for run701c to end"
RC=$W/run701c.log
until newline "$T0" "$RC" " END$" && ! alive experiments/701-kup-callset/run701c.sh; do
  capped && { log "run701c never ended by $CAP: stopping"; exit 1; }
  sleep 60
done
log "run701c ended: $(tail -2 "$RC" | head -1 | cut -c1-200)"
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
M=cpm-n4600
for s in 0 1 2; do
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
log "reading (EXPLORATORY)"
uv run python $X/reach701d.py > "$W/reach701d.out" 2>&1
log "reach701d rc=$?: $(/usr/bin/grep -E '^f |4600' "$W/reach701d.out" | tr '\n' ' ' | cut -c1-700)"
