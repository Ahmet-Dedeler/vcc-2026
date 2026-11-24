#!/bin/bash
# 618 (4d0207's rule, 602 README f75e591a): before ref_tv's quiet start, compare the tv operator's END peakRSS with its
# declared 14 GB. If it is <= 14, do nothing (ref_tv goes as planned at 17). If it is above 14, the extrapolation is running
# hot: PAUSE chain618ab4 (SIGSTOP; it cannot be edited while live) so ref_tv can be re-declared by the same ratio by hand.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/watch618ref.log; G=experiments/602-council-gate/work/gate602.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
CH=$(cat $X/work/chain618a.pid); log "watching for END 618-tv-operator (chain $CH)"
while :; do
  line=$(/usr/bin/grep -E "END 618-tv-operator " $G | tail -1)
  [ -n "$line" ] && break
  ps -p "$CH" >/dev/null 2>&1 || { log "chain exited before the tv operator ended: nothing to do"; exit 0; }
  sleep 3
done
pk=$(printf '%s' "$line" | sed -E 's/.*peakRSS=([0-9]+\.[0-9]+).*/\1/')
log "tv operator END: peakRSS $pk GB (declared 14)"
if awk -v p="$pk" 'BEGIN{exit !(p > 14.0)}'; then
  kill -STOP "$CH" && log "PAUSED chain618ab4 (SIGSTOP): tv operator ran hot ($pk > 14); re-declare ref_tv at 17 x $pk / 14 before resuming"
else
  log "tv operator within its declaration: ref_tv proceeds at 17 as planned"
fi
