#!/usr/bin/env bash
# Snapshot the public leaderboard into data/leaderboard/ on a fixed cadence.
#
# Why (2026-09-04): experiment 215 decomposed every team's board history into channels and could
# resolve every question except one -- whether cqawesome's x2.06 cosine doubling was one arm or
# five -- because FIVE of their submissions fell inside a single 2.5-day snapshot gap. 215's own
# recommendation: "pull the board every 6 h instead of every 12-24. Fix is free."
#
# At 2 submissions/team/day, a 3-hour cadence resolves every arm individually for every team.
# Snapshots are ~1 MB; 60 days at 3 h is ~480 files, under 500 MB.
#
#   nohup bash scripts/board_snapshot.sh > /dev/null 2>&1 &
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
INTERVAL="${BOARD_SNAPSHOT_INTERVAL:-10800}"   # 3 h
mkdir -p data/leaderboard
while true; do
  TS=$(date -u +%Y-%m-%dT%H%MZ)
  OUT="data/leaderboard/${TS}.json"
  if [ ! -f "$OUT" ]; then
    if curl -sf --max-time 120 https://virtualcellchallenge.org/api/leaderboard -o "$OUT.tmp"; then
      # only keep it if it parses and is non-trivial, so a truncated fetch never poisons the series
      if python3 -c "import json,sys; d=json.load(open('$OUT.tmp')); sys.exit(0 if d.get('num_entries',0)>100 else 1)" 2>/dev/null; then
        mv "$OUT.tmp" "$OUT"
      else
        rm -f "$OUT.tmp"
      fi
    else
      rm -f "$OUT.tmp"
    fi
  fi
  sleep "$INTERVAL"
done
