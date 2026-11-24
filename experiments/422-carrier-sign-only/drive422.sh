#!/bin/bash
# 422 serial driver under the AGENTS.md memory gate. Known answer first: base seed 0 must reproduce 416's base_s0 (|d| <= 1e-12).
cd "$(dirname "$0")" || exit 1
for line in "base 0" "s000 0" "s030 0" "s000 1" "s030 1" "s000 2" "s030 2"; do
  while true; do
    free=$(memory_pressure | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free/{f=$3} /Pages inactive/{i=$3} END{gsub("\\.","",f); gsub("\\.","",i); printf "%d", (f+i)*16384/1073741824}')
    if [ "${free:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]; then break; fi
    echo "$(date -u +%FT%TZ) gate hold free=$free% free+inactive=${fi}GiB"; sleep 60
  done
  name=$(echo "$line" | tr ' ' '_')
  echo "$(date -u +%FT%TZ) START $line"
  /usr/bin/time -l uv run python run422.py $line > "work/logs/$name.log" 2>&1
  echo "$(date -u +%FT%TZ) END rc=$? maxrss=$(awk '/maximum resident set size/{print $1}' "work/logs/$name.log") $line"
  if [ "$line" = "base 0" ]; then
    uv run python -c "
import json,sys
a=json.load(open('work/h1/metrics-base_s0.json'))['raw']; b=json.load(open('../416-de-carrier-menu/work/h1/metrics-base_s0.json'))['raw']
d=max(abs(a[k]-b[k]) for k in a); print('KNOWN ANSWER base_s0 vs 416 base_s0 max|d| =', d); sys.exit(0 if d <= 1e-12 else 1)" || { echo "KNOWN ANSWER FAILED: abort"; exit 1; }
  fi
done
echo "$(date -u +%FT%TZ) queue done"
