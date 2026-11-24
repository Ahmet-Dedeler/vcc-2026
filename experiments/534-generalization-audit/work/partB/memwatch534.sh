#!/bin/bash
# 534 Part B: sample memory every 30 s (evidence for the gate's reachability). Light. Stops at 2026-09-30T00:30Z.
OUT=~/Code/virtual-cell-challenge-2026/experiments/534-generalization-audit/work/partB/memwatch534.log
while [ "$(date -u +%Y%m%d%H%M)" -lt 202609300030 ]; do
  pct=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
  vm_stat | awk -v t="$(date -u +%FT%TZ)" -v p="$pct" '/Pages free/{f=$3}/Pages inactive/{i=$3}/Pages active/{a=$3}/Pages wired/{w=$4}/occupied by compressor/{c=$5}END{gsub(/\./,"",f);gsub(/\./,"",i);gsub(/\./,"",a);gsub(/\./,"",w);gsub(/\./,"",c);g=16384/2^30;printf "%s pressure_free %s%% free+inactive %.1f free %.1f inactive %.1f active %.1f wired %.1f compressor %.1f\n",t,p,(f+i)*g,f*g,i*g,a*g,w*g,c*g}' >> "$OUT"
  sleep 30
done
