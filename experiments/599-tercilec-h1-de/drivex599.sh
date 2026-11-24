#!/bin/bash
# 599 extension (preregistered): t and tc for seeds 3 and 4, one gate602 ticket per run (8 GB). Starts only once
# prepx599.npz (the exact K of those seeds) exists. Ends with a terminal "<ts> END" line on every path. Named drivex599,
# not chain5NN: build guards treat a "chain5NN*.sh" process as a running build.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/599-tercilec-h1-de; L=$X/work/logs/drivex599.log; mkdir -p $X/work/logs
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
ARGS="base %s --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
POST="bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs"
echo $$ > $X/work/drivex599.pid
log "drivex599 start (pid $$)"
[ -f $X/work/prepx599.npz ] || { log "no prepx599.npz (the registered K of seeds 3-4) -- stopping"; log "END"; exit 1; }
for s in 3 4; do
  for mode in t tc; do
    if [ -f $X/work/h1_$mode/metrics-base_s${s}_$POST.json ]; then log "skip $mode s$s (done)"; continue; fi
    log "run $mode s$s"
    bash experiments/602-council-gate/gate602.sh "599-$mode-s$s" 8 -- uv run python $X/run599.py $mode $(printf "$ARGS" "$s") \
      > $X/work/logs/$mode-s$s.log 2>&1
    log "done $mode s$s rc=$? ($(tail -1 $X/work/logs/$mode-s$s.log | cut -c1-100))"
  done
done
log "END"
