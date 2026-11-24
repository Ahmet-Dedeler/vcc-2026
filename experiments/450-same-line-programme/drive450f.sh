#!/usr/bin/env bash
# 450f: the best w for a WEAKER same-line bank, weakened by random sign flips (p 0.2 -> gate signs ~0.66, p 0.3 -> ~0.61;
# clean 0.764). Full mode (rows + per-target prior, CLEAN4 default), w 0.5 / 1 / 2, seed 0. Serial.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; LOG=$X/work/drive450f.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
for p in 0.2 0.3; do for w in 0.5 1 2; do lf="$LOGD/f_p${p}_w${w}.log"; log "START flip $p w $w"
  RUN450_FLIP=$p RUN450_NOISE_SEED=0 RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450c.py "$w" clean4 base 0 $A > "$lf" 2>&1 \
    && log "END flip $p w $w" || log "FAILED flip $p w $w (continuing)"; done; done
log "QUEUE DONE"
