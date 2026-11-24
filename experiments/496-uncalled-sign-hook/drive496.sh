#!/bin/bash
# 496 serial H1 driver (435s drive435.sh; RUN496_EMITTER=hook; known answer checked by ka496.py)
# aborts the queue if a MUST known answer fails.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/496-uncalled-sign-hook"
cd "$REPO" || exit 1
Q="$1"; LOGD="$HERE/work/logs"; mkdir -p "$LOGD"
gate(){
  while true; do
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ] && { echo "$(date -u +%FT%TZ) gate ok: free ${fp}% free+inactive ${fi} GiB"; return 0; }
    echo "$(date -u +%FT%TZ) gate hold: free ${fp}% free+inactive ${fi} GiB"; sleep 60
  done
}
while IFS= read -r line; do
  [ -z "$line" ] && continue
  case "$line" in \#*) continue;; esac
  gate
  tag=$(echo "$line" | tr ' -' '__')
  echo "$(date -u +%FT%TZ) START $line"
  /usr/bin/time -l uv run python "$HERE/run496.py" $line > "$LOGD/$tag.log" 2>&1
  rc=$?
  echo "$(date -u +%FT%TZ) END rc=$rc $line  peakRSS=$(grep 'maximum resident' "$LOGD/$tag.log" | awk '{printf "%.2f GB", $1/1e9}')"
  [ $rc -ne 0 ] && { echo "RUN FAILED: abort"; tail -5 "$LOGD/$tag.log"; exit 1; }
done < "$Q"
echo "$(date -u +%FT%TZ) QUEUE DONE"
