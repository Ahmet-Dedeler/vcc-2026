#!/usr/bin/env bash
# 475: 484's q80 finding priced on H1 (production emitter vs the patched copy). Serial after chain474. `orig` seed 0 is the
# known answer (stats == 435 seed 0 bit for bit); then `fix` x seeds 0-2; then the summary with target-bootstrap intervals.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/475-q80-protect; LOGD=$X/work/logs; LOG=$X/work/chain475.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain475.pid
log "armed (pid $$): waiting for chain474 to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END" experiments/474-october-hardening/work/chain474.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain474\.sh|def_day\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271600 ] && { log "474 never ended by 16:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local m=$1 s=$2; log "START $m seed $s"; /usr/bin/time -l uv run python $X/run475.py "$m" base "$s" $ARGS > "$LOGD/$m-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $m seed $s $(grep '^\[475\]' "$LOGD/$m-s$s.log" | tail -1 | cut -c1-120)"; return $rc; }
run orig 0 || { log "orig failed"; exit 1; }
uv run python - <<'PY' >> "$LOG"
import numpy as np, glob
a=np.load(glob.glob("experiments/475-q80-protect/work/h1_orig/stats-base_s0_*coarse.npz")[0])
b=np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA475 (orig == 435 seed 0, stats bit-equal):", "PASS" if ok else "FAIL")
open("experiments/475-q80-protect/work/ka475.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
/usr/bin/grep -q PASS $X/work/ka475.txt || { log "KA FAIL -- stopping"; exit 1; }
for s in 0 1 2; do run fix $s || exit 1; done
uv run python $X/summary475.py >> "$LOG" 2>&1
log "END"
