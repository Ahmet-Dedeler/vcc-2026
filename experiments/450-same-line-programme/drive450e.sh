#!/usr/bin/env bash
# 450e: the best w for a WEAKER same-line bank. Pan rows + noise (x1: gate signs 0.664; x2: 0.609; clean 0.764), full mode
# (rows + per-target prior, CLEAN4 default), w 0.5 / 1 / 2, seed 0. Serial; the machine is otherwise idle.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; LOG=$X/work/drive450e.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
for n in 1 2; do for w in 0.5 1 2; do lf="$LOGD/e_n${n}_w${w}.log"; log "START noise $n w $w"
  RUN450_NOISE=$n RUN450_NOISE_SEED=0 RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450c.py "$w" clean4 base 0 $A > "$lf" 2>&1 \
    || { log "failed noise $n w $w"; exit 1; }
  log "END noise $n w $w"; done; done
log "QUEUE DONE"
