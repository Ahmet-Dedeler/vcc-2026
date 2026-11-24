#!/bin/bash
# 493 runner: one panel/fold per process (comp493.py, ~6-7 GiB like comp438), strictly serial. Memory agreement
# (2026-09-26): before 04:30Z start only at vm_stat free+inactive >= 33 GiB; after 04:30Z at >= 26 GiB. Waits are capped
# by a timestamp. macOS /bin/bash 3.2. Log: work/run493.log.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/493-pertarget-q-october/work
LOG=$W/run493.log
echo $$ > $W/run493.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
OPEN=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T04:30:00Z" +%s)
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T12:00:00Z" +%s)
echo "$(date -u +%FT%TZ) runner start pid $$" >> $LOG
for pf in abc/orion_hct116 efp/orion_hct116 efr/orion_hct116 real/orion_hct116 efp/k562_gwps efr/k562_gwps real/k562_gwps efp/orion_hek293t efr/orion_hek293t real/orion_hek293t abc/k562_gwps abc/orion_hek293t; do
  p=${pf%%/*}; f=${pf##*/}
  if [ -f "$W/comp493-$p-$f.json" ]; then echo "$(date -u +%TZ) $p/$f exists, skip" >> $LOG; continue; fi
  while :; do
    now=$(date +%s); fi=$(fi_gib)
    if [ "$fi" -ge 33 ]; then break; fi
    if [ "$now" -ge "$OPEN" ] && [ "$fi" -ge 26 ]; then break; fi
    if [ "$now" -ge "$DEADLINE" ]; then echo "$(date -u +%TZ) deadline before $p/$f, stopping" >> $LOG; exit 0; fi
    sleep 60
  done
  echo "$(date -u +%TZ) start $p/$f (free+inactive ${fi} GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/493-pertarget-q-october/comp493.py "$p" "$f" >> "$W/comp493-$p-$f.log" 2>&1
  echo "$(date -u +%TZ) end $p/$f rc $?" >> $LOG
done
echo "$(date -u +%TZ) RUN493 END" >> $LOG
