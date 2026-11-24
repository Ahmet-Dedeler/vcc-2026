#!/usr/bin/env bash
# 450b H1 queue (same-line per-target rows in the effect), serial, after drive450's QUEUE DONE. KA first.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; mkdir -p "$LOGD"; LOG=$X/work/drive450b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q "QUEUE DONE\|KA FAIL\|failed\|stopping\|never ended" $X/work/drive450.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /drive450.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251600 ] && { log "drive450 never ended by 16:00Z"; exit 1; }; sleep 30; done
grep -q "QUEUE DONE" $X/work/drive450.log || { log "drive450 did not finish cleanly -- stopping"; exit 1; }
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local w=$1 p=$2 s=$3; local lf="$LOGD/b_w${w}_${p}_s${s}.log"; log "START w=$w prior=$p seed=$s"
       RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450b.py "$w" "$p" base "$s" $A > "$lf" 2>&1; local rc=$?
       log "END rc=$rc w=$w prior=$p seed=$s peakRSS=$(grep 'maximum resident' "$lf" | awk '{printf "%.2f GB", $1/1e9}')"; return $rc; }
run 0 clean4 0 || exit 1
uv run python - <<'PY' | tee -a "$LOG"
import numpy as np
a=np.load("experiments/450-same-line-programme/work/h1b_w0_clean4/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
b=np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA450b (w 0 through the blend hook == 435 seed 0):", "PASS" if ok else "FAIL")
open("experiments/450-same-line-programme/work/ka450b.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
grep -q PASS $X/work/ka450b.txt || { log "KA FAIL -- stopping"; exit 1; }
for w in 0.5 1 2; do run $w clean4 0 || exit 1; done
run 1 pan 0 || exit 1
for s in 1 2; do run 1 clean4 $s || exit 1; run 1 pan $s || exit 1; done
log "QUEUE DONE"
