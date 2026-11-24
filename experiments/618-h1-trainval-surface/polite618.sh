#!/bin/bash
# 618: enter gate602 for a BIG job (>= 12 GB ticket) only when it can start at once, so it never heads the FIFO queue and
# blocks the council's smaller jobs (19:30-19:51Z the 14 GB operator ticket did exactly that). Waits outside the queue
# until free+inactive >= max_gb + 13 GiB AND at most one gate slot is held, then execs gate602.sh with the same arguments.
#   bash polite618.sh <label> <max_gb> -- <command ...>
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
LABEL=$1; GB=$2; shift 2
need=$((GB + 13))
while :; do
  fi=$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')
  held=$(ls -d experiments/602-council-gate/work/slot*.lock 2>/dev/null | while read d; do o=$(cat "$d"/owner 2>/dev/null | awk '{print $1}'); [ -n "$o" ] && ps -p "$o" >/dev/null 2>&1 && echo x; done | wc -l | tr -d ' ')
  if [ "$fi" -ge "$need" ] && [ "$held" -le 1 ]; then break; fi
  sleep 60
done
printf '%s POLITE %s entering gate602: fi %s GiB, slots held %s\n' "$(date -u +%FT%TZ)" "$LABEL" "$fi" "$held" >> experiments/618-h1-trainval-surface/work/logs/polite618.log
exec bash experiments/602-council-gate/gate602.sh "$LABEL" "$GB" "$@"
