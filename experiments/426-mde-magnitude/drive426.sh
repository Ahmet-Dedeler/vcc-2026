#!/bin/bash
# 426 serial H1 driver: each queue line = arguments to run426.py; memory gate before each run; peak RSS via time -l.
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 1
Q="$1"; LOGD="$REPO/experiments/426-mde-magnitude/work/logs"; mkdir -p "$LOGD"
gate(){
  while true; do
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ] && return 0
    echo "$(date -u +%FT%TZ) gate hold: free ${fp}% free+inactive ${fi} GiB"; sleep 60
  done
}
while IFS= read -r line; do
  [ -z "$line" ] && continue
  gate
  tag=$(echo "$line" | tr ' -' '__')
  echo "$(date -u +%FT%TZ) START $line"
  /usr/bin/time -l uv run python experiments/426-mde-magnitude/run426.py $line > "$LOGD/$tag.log" 2>&1
  echo "$(date -u +%FT%TZ) END rc=$? $line  peakRSS=$(grep 'maximum resident' "$LOGD/$tag.log" | awk '{printf "%.2f GB", $1/1e9}')"
done < "$Q"
echo "$(date -u +%FT%TZ) QUEUE DONE"
