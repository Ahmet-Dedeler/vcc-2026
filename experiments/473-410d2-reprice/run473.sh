#!/bin/bash
# 473b runner: one PLOLO fold per process (arms473.py, ~6.5 GiB like comp420), strictly serial. Memory agreement with the
# peer/coordinator sessions (2026-09-26): before 04:30Z a process over 4 GB starts only when vm_stat free+inactive >= 33 GiB
# (the coordinator's builds and 474's def_day rehearsals re-check memory and do not wait); after 04:30Z the usual <= 10 GB
# rule, and this runner still keeps >= 26 GiB free+inactive at start. Waits are capped by a timestamp, never a count.
# Written for macOS /bin/bash 3.2 (no associative arrays). Log: work/run473.log; pid: work/run473.pid.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/473-410d2-reprice/work
LOG=$W/run473.log
echo $$ > $W/run473.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
OPEN=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T04:30:00Z" +%s)
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T16:00:00Z" +%s)
echo "$(date -u +%FT%TZ) runner start pid $$" >> $LOG
for f in orion_hct116 k562_gwps orion_hek293t kolf2; do
  if [ -f "$W/arms473-$f.json" ]; then echo "$(date -u +%TZ) $f exists, skip" >> $LOG; continue; fi
  while :; do
    now=$(date +%s); fi=$(fi_gib)
    if [ "$fi" -ge 33 ]; then break; fi
    if [ "$now" -ge "$OPEN" ] && [ "$fi" -ge 26 ]; then break; fi
    if [ "$now" -ge "$DEADLINE" ]; then echo "$(date -u +%TZ) deadline reached before $f, stopping" >> $LOG; exit 0; fi
    sleep 60
  done
  echo "$(date -u +%TZ) start $f (free+inactive ${fi} GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/473-410d2-reprice/arms473.py "$f" >> "$W/arms473-$f.log" 2>&1
  echo "$(date -u +%TZ) end $f rc $?" >> $LOG
done
echo "$(date -u +%TZ) RUN473 FOLDS END" >> $LOG
