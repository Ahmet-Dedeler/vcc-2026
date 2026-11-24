#!/bin/bash
# 537: serial, gated driver for run537.py (real production-emitter H1 runs, ~7.5 GB each, ~5 min).
# The queue is work/queue537.txt, re-read before every run (so it can grow without editing this live script):
#   one run per line: "<mode> <seed> [ka]"; blank lines and lines starting with # are skipped.
# A run is skipped when its zs metrics file exists. A line tagged "ka" runs ka537.py after the run; a FAIL writes
# work/KA-FAIL and stops the driver (nothing downstream is read until the known answer reproduces).
# Gate before EVERY run (COORDINATOR AMENDMENT 1, 2026-09-30 11:45Z, replacing the mandate's unreachable 33 GiB):
#   memory_pressure free >= 40 % AND vm_stat free+inactive >= 24 GiB AND swap free >= 400 MiB AND
#   def_pipeline.other_builders() == 0 AND no other run537.py alive (one heavy job at a time), and never between
#   10-05 01:30Z and 04:30Z (mandate).
# Peak RSS <= 8 GB: a watchdog samples the run's python RSS every 5 s and kills it above 8.0e9 bytes (marked FAILED).
# Stop: touch work/STOP537. Cap: 2026-10-05 18:00Z.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/537-arm-pipeline"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive537.log"; Q="$HERE/work/queue537.txt"
echo $$ > "$HERE/work/drive537.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo args= | awk -v c="$1" -v a="$2" '$1 ~ c && $0 ~ a {n++} END {exit n ? 0 : 1}'; }
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
MET="metrics-base_s%s_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"
RSS_CAP_KB=7812500     # 8.0e9 bytes in KiB (ps reports rss in KiB)
log "armed (pid $$), amendment-1 gate: mp>=40% fi>=24GiB swapfree>=400MiB builders==0"
last_why=""; last_t=0
gate(){
  while true; do
    [ -e "$HERE/work/STOP537" ] && { log "STOP537 present"; log "END"; exit 0; }
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202610051800 ] && { log "past the 10-05 18:00Z cap"; log "END"; exit 0; }
    why=""
    if [ "$now" -ge 202610050130 ] && [ "$now" -lt 202610050430 ]; then why="$why window-10-05"; fi
    alive "python[0-9.]*$" "537-arm-pipeline/run537[.]py" && why="$why run537-alive"
    nb=$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)
    [ "${nb:-9}" = "0" ] || why="$why builders${nb}"
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    sw=$(sysctl -n vm.swapusage | awk '{for(i=1;i<=NF;i++) if($i=="free") {v=$(i+2); u=substr(v,length(v)); x=substr(v,1,length(v)-1); if(u=="G") x*=1024; printf "%d", x}}')
    [ "${fp:-0}" -ge 40 ] || why="$why free${fp}%"
    [ "${fi:-0}" -ge 24 ] || why="$why fi${fi}<24GiB"
    [ "${sw:-0}" -ge 400 ] || why="$why swapfree${sw}MiB"
    if [ -z "$why" ]; then log "gate ok: free ${fp}% free+inactive ${fi} GiB swap free ${sw} MiB builders 0"; return 0; fi
    t=$(date +%s)
    if [ "$why" != "$last_why" ] || [ $((t - last_t)) -ge 1800 ]; then log "gate closed:$why"; last_why="$why"; last_t=$t; fi
    sleep 120
  done
}
watchdog(){
  # $1 = the tag of the run; kills the run537 python if its RSS exceeds the cap; writes the peak seen
  local tag=$1 peak=0 r pid
  while true; do
    pid=$(ps -axo pid=,args= | awk '$2 ~ /python[0-9.]*$/ && /537-arm-pipeline\/run537[.]py/ {print $1; exit}')
    if [ -z "$pid" ]; then sleep 5; [ -e "$HERE/work/.run_done" ] && break; continue; fi
    r=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d ' ')
    [ -n "$r" ] && [ "$r" -gt "$peak" ] && peak=$r
    if [ -n "$r" ] && [ "$r" -gt "$RSS_CAP_KB" ]; then
      log "WATCHDOG: $tag rss ${r} KiB > cap ${RSS_CAP_KB} KiB: killing pid $pid"; kill "$pid"; echo "$peak" > "$HERE/work/.peak"; break
    fi
    echo "$peak" > "$HERE/work/.peak"
    [ -e "$HERE/work/.run_done" ] && break
    sleep 5
  done
}
next_line(){
  [ -f "$Q" ] || return 1
  while read -r mode seed tag; do
    case "$mode" in ""|\#*) continue;; esac
    m="$HERE/work/h1_$mode/$(printf "$MET" "$seed")"
    [ -s "$m" ] && continue
    [ -e "$HERE/work/h1_$mode/FAILED_s$seed" ] && continue
    echo "$mode $seed $tag"; return 0
  done < "$Q"
  return 1
}
while true; do
  line=$(next_line) || { log "queue empty; waiting"; sleep 600; [ -e "$HERE/work/STOP537" ] && { log "END"; exit 0; }; continue; }
  set -- $line; mode=$1; seed=$2; tag=$3
  gate
  line2=$(next_line) || continue
  [ "$line2" = "$line" ] || continue
  log "START $mode seed $seed"
  rm -f "$HERE/work/.run_done" "$HERE/work/.peak"
  watchdog "$mode-s$seed" &
  WD=$!
  /usr/bin/time -l uv run python "$HERE/run537.py" "$mode" base "$seed" $ARGS > "$LOGD/$mode-s$seed.log" 2>&1; rc=$?
  touch "$HERE/work/.run_done"; wait "$WD" 2>/dev/null
  pk=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-s$seed.log" | awk '{printf "%.2f GB", $1/1e9}')
  wpk=$(cat "$HERE/work/.peak" 2>/dev/null)
  log "END rc=$rc $mode seed $seed peakRSS(time -l)=$pk watchdog-peak=$(awk -v k="${wpk:-0}" 'BEGIN{printf "%.2f GB", k*1024/1e9}')"
  if [ $rc -ne 0 ]; then mkdir -p "$HERE/work/h1_$mode"; touch "$HERE/work/h1_$mode/FAILED_s$seed"; continue; fi
  if [ "$tag" = "ka" ]; then
    out=$(uv run python "$HERE/ka537.py" "$mode" "$seed" 2>&1); krc=$?
    while IFS= read -r l; do log "$l"; done <<< "$out"
    if [ $krc -ne 0 ]; then log "KNOWN ANSWER FAILED: stopping"; touch "$HERE/work/KA-FAIL"; log "END"; exit 1; fi
  fi
done
