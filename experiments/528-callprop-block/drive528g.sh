#!/bin/bash
# 528g: after drive528f.sh: more seeds for the two interaction pairs that decide the recommendation -- the kept block on top
# of 527's hook (446bj2zmz's DE construction; seed 0 read +0.0083 at f 1) and at 523's K 1,600. Same gates as 528e/f.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/528-callprop-block"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive528g.log"
echo $$ > "$HERE/work/drive528g.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo args= | awk -v c="$1" -v a="$2" '$1 ~ c && $0 ~ a {n++} END {exit n ? 0 : 1}'; }
log "armed (pid $$): waiting for drive528e.sh and drive528f.sh to exit"
while alive "bash$" "drive528[ef][.]sh"; do sleep 60; done
log "drive528e/f exited"
gate(){
  while true; do
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202609300600 ] && { log "past 06:00Z 09-30"; log "END"; exit 0; }
    why=""
    alive "python[0-9.]*$" "517-flex-programme/run517[.]py" && why="$why run517"
    alive "python[0-9.]*$" "528-callprop-block/run528[.]py" && why="$why run528"
    need=23
    alive "python[0-9.]*$" "make_submission_pm_cd4[.]py|def_pipeline[.]py" && need=25
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] || why="$why free${fp}%"
    [ "${fi:-0}" -ge "$need" ] || why="$why fi${fi}<${need}GiB"
    [ -z "$why" ] && { log "gate ok: free ${fp}% free+inactive ${fi} GiB (need ${need})"; return 0; }
    sleep 60
  done
}
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
run(){
  local mode=$1 s=$2
  local m="$HERE/work/h1_$mode/metrics-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"
  gate
  [ -s "$m" ] && { log "SKIP $mode seed $s (metrics exist)"; return 0; }
  log "START $mode seed $s"
  /usr/bin/time -l uv run python "$HERE/run528.py" "$mode" base "$s" $ARGS > "$LOGD/$mode-s$s.log" 2>&1; local rc=$?
  log "END rc=$rc $mode seed $s peakRSS=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-s$s.log" | awk '{printf "%.2f GB", $1/1e9}')"
}
for s in 1 2; do
  run "cpm+mzg" "$s"
  run "lcap-kept+mzg" "$s"
done
run "cpm-n1200" 1
run "lcap-kept-n1200" 1
log "QUEUE DONE"
log "END"
