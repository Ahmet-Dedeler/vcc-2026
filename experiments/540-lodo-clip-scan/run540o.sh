#!/bin/bash
# 540o runner: comp540o.py per panel/fold (~6 GiB, ~30 s), strictly serial, each at vm_stat free+inactive >= 26 GiB;
# then the 540c PLOLO interaction (arms540c.py per fold) and its pricer. Capped at 12:00Z. macOS /bin/bash 3.2.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/540-lodo-clip-scan/work
LOG=$W/run540o.log
echo $$ > $W/run540o.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T12:00:00Z" +%s)
gate() { while :; do fi=$(fi_gib); [ "$fi" -ge 26 ] && return 0; [ "$(date +%s)" -ge "$DEADLINE" ] && return 1; sleep 60; done; }
echo "$(date -u +%FT%TZ) runner start pid $$" >> $LOG
for pf in efr/orion_hct116 real/orion_hct116 efp/orion_hct116 efr/k562_gwps real/k562_gwps efp/k562_gwps efr/orion_hek293t real/orion_hek293t efp/orion_hek293t; do
  p=${pf%%/*}; f=${pf##*/}
  if [ -f "$W/comp540o-$p-$f.json" ]; then echo "$(date -u +%TZ) $p/$f exists, skip" >> $LOG; continue; fi
  gate || { echo "$(date -u +%TZ) deadline before $p/$f, stopping" >> $LOG; exit 0; }
  echo "$(date -u +%TZ) start $p/$f (free+inactive $(fi_gib) GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/540-lodo-clip-scan/comp540o.py "$p" "$f" >> "$W/comp540o-$p-$f.log" 2>&1
  echo "$(date -u +%TZ) end $p/$f rc $?" >> $LOG
done
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  [ -f experiments/540-lodo-clip-scan/arms540c.py ] || break
  if [ -f "$W/arms540c-$f.json" ]; then echo "$(date -u +%TZ) 540c $f exists, skip" >> $LOG; continue; fi
  gate || { echo "$(date -u +%TZ) deadline before 540c $f, stopping" >> $LOG; exit 0; }
  echo "$(date -u +%TZ) start 540c $f (free+inactive $(fi_gib) GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/540-lodo-clip-scan/arms540c.py "$f" > "$W/arms540c-$f.log" 2>&1
  echo "$(date -u +%TZ) end 540c $f rc $?" >> $LOG
done
echo "$(date -u +%TZ) RUN540O END" >> $LOG
