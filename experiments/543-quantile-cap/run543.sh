#!/bin/bash
# 543 runner: arms543.py on the 4 A/B/C folds, then comp543o.py on 9 October-like panel/folds; strictly serial, each at
# vm_stat free+inactive >= 26 GiB. Capped at 12:00Z. macOS /bin/bash 3.2.
R=~/Code/virtual-cell-challenge-2026
cd "$R" || exit 1
W=experiments/543-quantile-cap/work
LOG=$W/run543.log
echo $$ > $W/run543.pid
fi_gib() { vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}'; }
DEADLINE=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-09-27T12:00:00Z" +%s)
gate() { while :; do [ "$(fi_gib)" -ge 26 ] && return 0; [ "$(date +%s)" -ge "$DEADLINE" ] && return 1; sleep 60; done; }
echo "$(date -u +%FT%TZ) runner start pid $$" >> $LOG
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  gate || { echo "$(date -u +%TZ) deadline, stopping" >> $LOG; exit 0; }
  echo "$(date -u +%TZ) start $f ($(fi_gib) GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/543-quantile-cap/arms543.py "$f" > "$W/arms543-$f.log" 2>&1
  echo "$(date -u +%TZ) end $f rc $?" >> $LOG
done
for pf in efr/orion_hct116 real/orion_hct116 efp/orion_hct116 efr/k562_gwps real/k562_gwps efp/k562_gwps efr/orion_hek293t real/orion_hek293t efp/orion_hek293t; do
  p=${pf%%/*}; f=${pf##*/}
  gate || { echo "$(date -u +%TZ) deadline, stopping" >> $LOG; exit 0; }
  echo "$(date -u +%TZ) start $p/$f ($(fi_gib) GiB)" >> $LOG
  /usr/bin/time -l uv run python experiments/543-quantile-cap/comp543o.py "$p" "$f" > "$W/comp543o-$p-$f.log" 2>&1
  echo "$(date -u +%TZ) end $p/$f rc $?" >> $LOG
done
echo "$(date -u +%TZ) RUN543 END" >> $LOG
