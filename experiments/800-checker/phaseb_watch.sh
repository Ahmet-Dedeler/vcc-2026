#!/bin/bash
# phaseb_watch.sh — the lab lead's event stream: new or changed PROPOSAL / DOSSIER / judge VERDICT files and every lab worker's
# alive/EXITED flip (all experiments/8NN-*/claude-run.pid, found each poll). One line per change, every 20 s. Read-only.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
shopt -s nullglob
snap(){
  cur=""
  for f in experiments/8[0-9][0-9]-*/PROPOSAL*.md experiments/8[0-9][0-9]-*/DOSSIER.md experiments/819-judge/*/VERDICT.md; do
    cur="$cur $f:$(stat -f %m "$f")"
  done
  for p in experiments/8[0-9][0-9]-*/claude-run.pid; do
    e=$(basename "$(dirname "$p")")
    # alive = a live process named claude (pids are recycled; 10-17: 805's stale pid became a `sleep 30`)
    if [ "$(ps -o comm= -p "$(cat "$p")" 2>/dev/null | tr -d ' ')" = "claude" ]; then cur="$cur $e:alive"; else cur="$cur $e:EXITED"; fi
  done
  echo "$cur"
}
prev=$(snap)
while true; do
  sleep 20; now=$(snap)
  if [ "$now" != "$prev" ]; then
    comm -13 <(echo "$prev" | tr " " "\n" | sort) <(echo "$now" | tr " " "\n" | sort) | sed "s/^/CHANGE /"
    prev=$now
  fi
done
