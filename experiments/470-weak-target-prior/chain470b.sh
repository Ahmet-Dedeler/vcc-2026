#!/usr/bin/env bash
# 470b: strong targets' programme from Arc's non-essential vs essential perturbations (H1). Serial after chain471.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/470-weak-target-prior; LOGD=$X/work/logs; LOG=$X/work/chain470b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain470b.pid
log "armed (pid $$): waiting for chain471 to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END( builder .*)?|.*stopping.*)$" experiments/471-h1full-union/work/chain471.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain471\.sh|run_arm_generic\.sh|run_ka_471\.sh|run_dump471\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "471 never ended by 12:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local m=$1 s=$2; log "START $m seed $s"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run470b.py "$m" base "$s" $ARGS > "$LOGD/$m-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $m seed $s $(grep '^\[470\] '"$m" "$LOGD/$m-s$s.log" | tail -1 | cut -c1-160)"; return $rc; }
for m in sarcne sarcess; do for s in 0 1 2; do run $m $s || exit 1; done; done
uv run python $X/summary470.py >> "$LOG" 2>&1
log "END"
