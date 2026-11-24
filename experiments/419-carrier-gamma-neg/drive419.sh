#!/bin/bash
# 419 serial driver: waits for 416's H1 driver to finish (one H1 harness at a time), then runs the queue under the
# AGENTS.md memory gate. First job is the known answer: base seed 0 must reproduce 416's metrics-base_s0.json exactly.
cd "$(dirname "$0")" || exit 1
while kill -0 "$(cat ../416-de-carrier-menu/work/drive416.pid 2>/dev/null)" 2>/dev/null; do sleep 60; done
echo "$(date -u +%FT%TZ) 416 driver finished; starting 419"
for line in "base 0" "gm030 0" "gm065 0" "gm030 1" "gm065 1" "gm030 2" "gm065 2"; do
  while true; do
    free=$(memory_pressure | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free/{f=$3} /Pages inactive/{i=$3} END{gsub("\\.","",f); gsub("\\.","",i); printf "%d", (f+i)*16384/1073741824}')
    if [ "${free:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]; then break; fi
    echo "$(date -u +%FT%TZ) gate hold free=$free% free+inactive=${fi}GiB"; sleep 60
  done
  name=$(echo "$line" | tr ' ' '_')
  echo "$(date -u +%FT%TZ) START $line"
  /usr/bin/time -l uv run python run419.py $line > "work/logs/$name.log" 2>&1
  echo "$(date -u +%FT%TZ) END rc=$? maxrss=$(awk '/maximum resident set size/{print $1}' "work/logs/$name.log") $line"
  if [ "$line" = "base 0" ]; then
    uv run python -c "
import json,sys
a=json.load(open('work/h1/metrics-base_s0.json'))['raw']; b=json.load(open('../416-de-carrier-menu/work/h1/metrics-base_s0.json'))['raw']
ok = a == b; print('KNOWN ANSWER base_s0 == 416 base_s0:', ok, a if not ok else ''); sys.exit(0 if ok else 1)" || { echo "KNOWN ANSWER FAILED: abort"; exit 1; }
  fi
done
echo "$(date -u +%FT%TZ) queue done"
