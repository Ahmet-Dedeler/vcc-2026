#!/bin/bash
# 416 serial driver. usage: bash drive416.sh QUEUEFILE
# Reads QUEUEFILE one line at a time (lines may be appended while it runs; '#' lines skipped), holds each job until the
# AGENTS.md memory gate passes (memory_pressure free >= 40 % AND vm_stat free+inactive >= 20 GiB), runs it under
# /usr/bin/time -l, logs peak RSS. A failing ka416.py line (known answer) aborts the queue.
cd "$(dirname "$0")" || exit 1
Q=$1; n=1
while true; do
  line=$(sed -n "${n}p" "$Q")
  if [ -z "$line" ]; then echo "$(date -u +%FT%TZ) queue exhausted at line $n"; break; fi
  n=$((n+1))
  case "$line" in \#*) continue;; esac
  while true; do
    free=$(memory_pressure | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free/{f=$3} /Pages inactive/{i=$3} END{gsub("\\.","",f); gsub("\\.","",i); printf "%d", (f+i)*16384/1073741824}')
    if [ "${free:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]; then break; fi
    echo "$(date -u +%FT%TZ) gate hold free=$free% free+inactive=${fi}GiB"; sleep 60
  done
  name=$(echo "$line" | tr ' ' '_' | tr -d '-' | sed 's/\.py//')
  echo "$(date -u +%FT%TZ) START $line (free=$free% free+inactive=${fi}GiB)"
  /usr/bin/time -l uv run python $line > "work/logs/$name.log" 2>&1
  rc=$?
  rss=$(awk '/maximum resident set size/{print $1}' "work/logs/$name.log")
  echo "$(date -u +%FT%TZ) END rc=$rc maxrss=$rss $line"
  case "$line" in ka416.py*) if [ $rc -ne 0 ]; then echo "KNOWN ANSWER FAILED: abort"; exit 1; fi;; esac
done
