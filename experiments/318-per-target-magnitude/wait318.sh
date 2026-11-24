#!/usr/bin/env bash
# block until monitor318.sh has collected all six agg_results.csv (or 3.5 h passes)
cd ~/Code/virtual-cell-challenge-2026
OUT=experiments/318-per-target-magnitude/cloud-output
for i in $(seq 1 210); do
  N=0
  for A in ptrule ptshuf1 ptshuf2 ptshuf3 ptamp6 ptunif; do
    [ -f "$OUT/s318-$A/$A/run/agg_results.csv" ] && N=$((N+1))
  done
  if [ "$N" -ge 6 ]; then echo "$(date -u +%FT%TZ) ALL 6 AGG PRESENT"; exit 0; fi
  sleep 60
done
echo "$(date -u +%FT%TZ) TIMEOUT with $N/6"
