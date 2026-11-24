#!/bin/bash
# (named drive599, not chain599: build guards treat any "chain5NN*.sh" process as a running build, AGENTS 6.2)
# 599: the known answer (mode ka, seed 0, vs 435's base run) and then the arms t / tc / s175 x seeds 0 1 2, serial, one
# gate602 ticket per run (the council's FIFO memory gate, 8 GB each). The arms start only once PREREGISTRATION.md is
# committed. Every path ends with a terminal "<ts> END" line. Re-runnable: finished runs are skipped.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/599-tercilec-h1-de; L=$X/work/logs/drive599.log; mkdir -p $X/work/logs
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
ARGS="base %s --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
POST="bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs"
echo $$ > $X/work/drive599.pid
log "drive599 start (pid $$)"
run(){ local mode=$1 s=$2 rc
  log "run $mode s$s"
  bash experiments/602-council-gate/gate602.sh "599-$mode-s$s" 8 -- uv run python $X/run599.py $mode $(printf "$ARGS" "$s") \
    > $X/work/logs/$mode-s$s.log 2>&1
  rc=$?
  log "done $mode s$s rc=$rc ($(tail -1 $X/work/logs/$mode-s$s.log | cut -c1-120))"
  return $rc; }
if [ ! -f $X/work/h1_ka/metrics-base_s0_$POST.json ]; then
  run ka 0 || { log "KA run failed -- stopping"; log "END"; exit 1; }
fi
uv run python $X/ka599.py 0 >> "$L" 2>&1 || { log "KA599 FAIL -- stopping"; log "END"; exit 1; }
log "KA599 PASS"
git log --oneline -1 -- $X/PREREGISTRATION.md | grep -q . || { log "PREREGISTRATION.md not committed -- stopping before the arms"; log "END"; exit 1; }
for s in 0 1 2; do
  for mode in t tc s175; do
    if [ -f $X/work/h1_$mode/metrics-base_s${s}_$POST.json ]; then log "skip $mode s$s (done)"; continue; fi
    run $mode $s || log "WARN $mode s$s rc != 0 (continuing)"
  done
done
bash experiments/602-council-gate/gate602.sh "599-score" 6 -- uv run python $X/score599.py > $X/work/logs/score599.log 2>&1
log "score599 rc=$? ($(grep -c . $X/work/logs/score599.log) lines)"
log "END"
