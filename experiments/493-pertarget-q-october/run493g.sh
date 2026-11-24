#!/bin/bash
# 493g runner: after 492g's runner (pid in 492's work/run492g.pid) has exited, one panel/fold per process (comp493g.py,
# ~6-7 GiB like comp438), strictly serial, each started at vm_stat free+inactive >= 26 GiB (after 04:30Z; the
# coordinator's 471c build runs ~04:15-04:50Z and the gate simply waits). Waits capped at 12:00Z. macOS /bin/bash 3.2.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/493-pertarget-q-october/work
LOG=$W/run493g.log
echo $$ > $W/run493g.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
OPEN=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T04:30:00Z" +%s)
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T12:00:00Z" +%s)
P492=$(cat experiments/492-per-target-sharpening/work/run492g.pid)
echo "$(date -u +%FT%TZ) runner start pid $$; waiting for 492g runner $P492" >> $LOG
while kill -0 "$P492" 2>/dev/null; do
  [ "$(date +%s)" -ge "$DEADLINE" ] && { echo "$(date -u +%TZ) 492g runner still alive at deadline, stopping" >> $LOG; exit 0; }
  sleep 30
done
for pf in efr/orion_hct116 real/orion_hct116 efp/orion_hct116 efr/k562_gwps real/k562_gwps efp/k562_gwps efr/orion_hek293t real/orion_hek293t efp/orion_hek293t; do
  p=${pf%%/*}; f=${pf##*/}
  if [ -f "$W/comp493g-$p-$f.json" ]; then echo "$(date -u +%TZ) $p/$f exists, skip" >> $LOG; continue; fi
  while :; do
    now=$(date +%s); fi=$(fi_gib)
    if [ "$now" -ge "$OPEN" ] && [ "$fi" -ge 26 ]; then break; fi
    if [ "$now" -ge "$DEADLINE" ]; then echo "$(date -u +%TZ) deadline before $p/$f, stopping" >> $LOG; exit 0; fi
    sleep 60
  done
  echo "$(date -u +%TZ) start $p/$f (free+inactive ${fi} GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/493-pertarget-q-october/comp493g.py "$p" "$f" >> "$W/comp493g-$p-$f.log" 2>&1
  echo "$(date -u +%TZ) end $p/$f rc $?" >> $LOG
done
uv run python experiments/493-pertarget-q-october/analyse493g.py > $W/analyse493g.log 2>&1
echo "$(date -u +%TZ) RUN493G END" >> $LOG
