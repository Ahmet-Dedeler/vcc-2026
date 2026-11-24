#!/usr/bin/env bash
# 450 H1 queue, serial, after 448b's lvl3 build ends (heavy jobs strictly serial). Known answer first; abort if it fails.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; mkdir -p "$LOGD"; LOG=$X/work/drive450.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/448-def-rehearsal-q175/work/run448b.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /run448b.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251500 ] && { log "448b never ended by 15:00Z"; exit 1; }; sleep 30; done
ARGS="base 0 --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local p=$1; shift; log "START $p $*"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450.py "$p" "$@" > "$LOGD/$p-$(echo "$*" | tr ' -' '__').log" 2>&1; local rc=$?
       log "END rc=$rc $p $* peakRSS=$(grep 'maximum resident' "$LOGD/$p-$(echo "$*" | tr ' -' '__').log" | awk '{printf "%.2f GB", $1/1e9}')"; return $rc; }
run clean4 $ARGS || { log "KA run failed"; exit 1; }
uv run python - <<'PY' | tee -a "$LOG"
import numpy as np
a=np.load("experiments/450-same-line-programme/work/h1_clean4/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
b=np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA450 (clean4 through the wrapper == 435 seed 0, stats bit-equal):", "PASS" if ok else "FAIL")
open("experiments/450-same-line-programme/work/ka450.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
grep -q PASS $X/work/ka450.txt || { log "KA FAIL -- stopping"; exit 1; }
for s in 0 1 2; do run pan base $s --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse || exit 1; done
run pan_c4fill $ARGS || exit 1
log "QUEUE DONE"
