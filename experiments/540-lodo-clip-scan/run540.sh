#!/bin/bash
# 540 PLOLO runner: after 493g's runner (pid in 493's work/run493g.pid) has exited, arms540.py per fold (~5.5 GiB, ~20 s),
# strictly serial, each started at vm_stat free+inactive >= 26 GiB and not before 04:30Z; then the pricer. Capped at
# 12:00Z. macOS /bin/bash 3.2.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/540-lodo-clip-scan/work
LOG=$W/run540.log
echo $$ > $W/run540.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
OPEN=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T04:30:00Z" +%s)
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T12:00:00Z" +%s)
P493=$(cat experiments/493-pertarget-q-october/work/run493g.pid)
echo "$(date -u +%FT%TZ) runner start pid $$; waiting for 493g runner $P493" >> $LOG
while kill -0 "$P493" 2>/dev/null; do
  [ "$(date +%s)" -ge "$DEADLINE" ] && { echo "$(date -u +%TZ) 493g runner still alive at deadline, stopping" >> $LOG; exit 0; }
  sleep 30
done
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  if [ -f "$W/arms540-$f.json" ]; then echo "$(date -u +%TZ) $f exists, skip" >> $LOG; continue; fi
  while :; do
    now=$(date +%s); fi=$(fi_gib)
    if [ "$now" -ge "$OPEN" ] && [ "$fi" -ge 26 ]; then break; fi
    if [ "$now" -ge "$DEADLINE" ]; then echo "$(date -u +%TZ) deadline before $f, stopping" >> $LOG; exit 0; fi
    sleep 60
  done
  echo "$(date -u +%TZ) start $f (free+inactive ${fi} GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/540-lodo-clip-scan/arms540.py "$f" > "$W/arms540-$f.log" 2>&1
  echo "$(date -u +%TZ) end $f rc $?" >> $LOG
done
uv run python experiments/540-lodo-clip-scan/analyse540.py > $W/analyse540.log 2>&1
echo "$(date -u +%TZ) RUN540 END" >> $LOG
