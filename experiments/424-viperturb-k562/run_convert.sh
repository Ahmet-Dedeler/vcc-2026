#!/bin/bash
# 424: convert one bin with a memory watchdog (kills R above LIMIT_GB RSS). Usage: run_convert.sh binC
set -u
cd "$(dirname "$0")"
b=$1; LIMIT_GB=${LIMIT_GB:-22}
mkdir -p work/conv/$b
/usr/bin/time -l Rscript convert_bin.R work/dl/genome_wide_$b.RDS work/conv/$b > work/conv/$b.log 2>&1 &
tp=$!
peak=0
while kill -0 $tp 2>/dev/null; do
  rp=$(pgrep -P $tp -f "R" | head -1)
  if [ -n "$rp" ]; then
    kb=$(ps -o rss= -p $rp 2>/dev/null | tr -d ' ')
    kb=${kb:-0}
    [ "$kb" -gt "$peak" ] && peak=$kb
    if [ "$kb" -gt $((LIMIT_GB*1024*1024)) ]; then
      echo "$(date -u +%FT%TZ) WATCHDOG: R rss ${kb} kB > ${LIMIT_GB} GB, killing" >> work/conv/$b.log
      kill $rp
    fi
  fi
  sleep 3
done
wait $tp; rc=$?
echo "$(date -u +%FT%TZ) rc=$rc watchdog_peak_rss_kB=$peak" >> work/conv/$b.log
