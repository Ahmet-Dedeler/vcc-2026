#!/usr/bin/env bash
# 450c H1 queue, serial, after drive450b's QUEUE DONE: the per-target same-line prior on top of the row blend (raw rows),
# the same with programme-centred rows (RUN450_CENTRE=1; signs still from raw rows), and the per-target prior with the
# CLEAN4 default (to split the programme's share). Re-armed 2026-09-25 03:1xZ with the centred arms added.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; LOG=$X/work/drive450c.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q "QUEUE DONE\|KA FAIL\|stopping\|never ended" $X/work/drive450b.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /drive450b\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251700 ] && { log "drive450b never ended by 17:00Z"; exit 1; }; sleep 30; done
grep -q "QUEUE DONE" $X/work/drive450b.log || { log "drive450b did not finish cleanly -- stopping"; exit 1; }
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local p=$1 s=$2 c=$3; local lf="$LOGD/c_w1_${p}_s${s}_c${c}.log"; log "START w=1 prior=$p seed=$s centred=$c (per-target same-line prior)"
       RUN450_CENTRE=$c RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450c.py 1 "$p" base "$s" $A > "$lf" 2>&1; local rc=$?
       log "END rc=$rc w=1 prior=$p seed=$s centred=$c peakRSS=$(grep 'maximum resident' "$lf" | awk '{printf "%.2f GB", $1/1e9}')"; return $rc; }
for s in 0 1 2; do run pan $s 0 || exit 1; run pan $s 1 || exit 1; done
run clean4 0 0 || exit 1
for s in 0 1 2; do lf="$LOGD/b_w1_pan_s${s}_centred.log"; log "START (450b) w=1 prior=pan seed=$s centred rows"
  RUN450_CENTRE=1 RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450b.py 1 pan base "$s" $A > "$lf" 2>&1 || { log "450b centred failed"; exit 1; }
  log "END (450b) centred seed=$s"; done
log "QUEUE DONE"
