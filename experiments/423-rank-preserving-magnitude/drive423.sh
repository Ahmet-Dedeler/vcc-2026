#!/bin/bash
# 423 serial driver under the AGENTS.md memory gate. One process per (seed, bulk): the base blocks are emitted once and every
# post-pass variant is scored from them. Known answer: base_s<seed> must reproduce 416's metrics-base_s<seed>[_bps].json
# (max|d| <= 1e-12) or the queue aborts.
# usage: bash drive423.sh QUEUEFILE   (lines: "SEED BULK VARIANTS")
cd "$(dirname "$0")" || exit 1
while read -r seed bulk variants; do
  [ -z "$seed" ] && continue
  case "$seed" in \#*) continue;; esac
  while true; do
    free=$(memory_pressure | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free/{f=$3} /Pages inactive/{i=$3} END{gsub("\\.","",f); gsub("\\.","",i); printf "%d", (f+i)*16384/1073741824}')
    if [ "${free:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]; then break; fi
    echo "$(date -u +%FT%TZ) gate hold free=$free% free+inactive=${fi}GiB"; sleep 60
  done
  sfx=""; [ "$bulk" = "ps" ] && sfx="_bps"
  name="s${seed}${sfx}"
  echo "$(date -u +%FT%TZ) START seed=$seed bulk=$bulk variants=$variants (free=$free% fi=${fi}GiB)"
  /usr/bin/time -l uv run python run423.py "$seed" --bulk "$bulk" --variants "$variants" > "work/logs/run_${name}.log" 2>&1
  rc=$?
  echo "$(date -u +%FT%TZ) END rc=$rc maxrss=$(awk '/maximum resident set size/{print $1}' "work/logs/run_${name}.log") $name"
  uv run python -c "
import json,sys
a=json.load(open('work/h1/metrics-base_${name}.json'))['raw']; b=json.load(open('../416-de-carrier-menu/work/h1/metrics-base_${name}.json'))['raw']
d=max(abs(a[k]-b[k]) for k in a); print('KNOWN ANSWER base_${name} vs 416 max|d| =', d); sys.exit(0 if d <= 1e-12 else 1)" || { echo "KNOWN ANSWER FAILED: abort"; exit 1; }
done < "$1"
echo "$(date -u +%FT%TZ) queue done"
