#!/bin/bash
# Resume lab workers 803 and 805 after the 23:30-01:45Z window (charter Amendment 1: at most 2 lab agents in it).
# Waits for the target time (capped: gives up after 2026-10-16T04:00Z), then resumes each by session id.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
T=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-10-16T01:47:00Z" +%s); CAP=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-10-16T04:00:00Z" +%s)
while [ "$(date -u +%s)" -lt "$T" ]; do sleep 30; done
[ "$(date -u +%s)" -lt "$CAP" ] || { echo "$(date -u +%FT%TZ) past cap, not resuming"; exit 0; }
for e in 803-de-engine 805-cd4-flex-surface; do
  if [ -f "experiments/$e/STOP" ]; then echo "$(date -u +%FT%TZ) $e has STOP, skipped"; continue; fi
  bash experiments/800-checker/lab_subagent.sh --resume "$e" experiments/800-checker/resume/amend1-later.md
  echo "$(date -u +%FT%TZ) resumed $e"
done
echo "$(date -u +%FT%TZ) END"
