#!/bin/bash
# 522 — gated, strictly serial runner for the nine heavy dumps (dump522.py, ~6 GB each) and their light arm screens.
# Gate before EVERY dump (the mandate): vm_stat free+inactive >= 33 GiB until 02:00Z 2026-09-29, then >= 20 GiB
# (AGENTS 9) with <= 10 GB per process. Never two dumps at once. Skips surfaces whose dump already exists.
# Writes only under experiments/522-october-bulk-key/work. Ends with "<ts> END" on every path.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/522-october-bulk-key/work
LOG=$W/run_dumps522.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
freeinact() { vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/1073741824}'; }
need() { if [ "$(date -u +%Y%m%d%H%M)" -ge 202609290200 ]; then echo 20; else echo 33; fi; }
echo "$(ts) start pid $$" >> $LOG
for spec in real:k562_gwps real:orion_hct116 real:orion_hek293t efr:k562_gwps efr:orion_hct116 efr:orion_hek293t efp:k562_gwps efp:orion_hct116 efp:orion_hek293t; do
  p=${spec%%:*}; f=${spec##*:}
  if [ ! -f $W/dump522-$p-$f.npz ]; then
    ok=0
    while true; do                       # the gate must hold on 2 consecutive checks 30 s apart (a dump takes ~15 s)
      fi_=$(freeinact); nd=$(need)
      if [ "$fi_" -ge "$nd" ]; then ok=$((ok+1)); else ok=0; fi
      if [ $ok -ge 2 ]; then break; fi
      sleep 30
    done
    echo "$(ts) dump $p $f (free+inactive ${fi_} GiB >= ${nd})" >> $LOG
    /usr/bin/time -l uv run python experiments/522-october-bulk-key/dump522.py $p $f > $W/dump522-$p-$f.log 2>&1
    rc=$?
    echo "$(ts) dump $p $f rc $rc $(grep 'maximum resident' $W/dump522-$p-$f.log | awk '{printf "%.2f GiB", $1/1073741824}')" >> $LOG
    if [ $rc -ne 0 ]; then echo "$(ts) dump FAILED $p $f; stopping" >> $LOG; echo "$(ts) END" >> $LOG; exit 1; fi
  fi
  if [ ! -f $W/arms522-$p-$f.json ]; then
    /usr/bin/time -l uv run python experiments/522-october-bulk-key/arms522.py $p $f > $W/arms522-$p-$f.log 2>&1
    rc=$?
    echo "$(ts) arms $p $f rc $rc $(grep 'maximum resident' $W/arms522-$p-$f.log | awk '{printf "%.2f GiB", $1/1073741824}')" >> $LOG
    if [ $rc -ne 0 ]; then echo "$(ts) arms FAILED $p $f; stopping" >> $LOG; echo "$(ts) END" >> $LOG; exit 1; fi
  fi
done
echo "$(ts) all nine surfaces dumped and screened" >> $LOG
echo "$(ts) END" >> $LOG
