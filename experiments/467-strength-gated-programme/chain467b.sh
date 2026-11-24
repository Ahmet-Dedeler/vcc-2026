#!/usr/bin/env bash
# 467b: does the incumbent CLEAN4 prior pay reach on weak targets too? c4off33/c4off50 (CLEAN4 on strong targets, the
# weak share keeps design signs via protect_q 0, MDE kept) and arcoff33 (Arc on strong, none on weak); 3 seeds each.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/467-strength-gated-programme; LOGD=$X/work/logs; LOG=$X/work/chain467b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local m=$1 s=$2; log "START $m seed $s"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run467.py "$m" base "$s" $ARGS > "$LOGD/$m-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $m seed $s $(grep '^\[467\]' "$LOGD/$m-s$s.log" | tail -1 | cut -c1-160)"; return $rc; }
for m in c4off33 c4off50 arcoff33; do for s in 0 1 2; do run $m $s || exit 1; done; done
uv run python $X/summary467.py >> "$LOG" 2>&1
log "END"
