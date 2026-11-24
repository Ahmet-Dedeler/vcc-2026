#!/usr/bin/env bash
# 462l: the H1 six-member profile of Arc's own programme as the block prior (upper reference for 458d), serial after 458e's
# build (heavy jobs strictly serial; each run ~3 min, ~8 GB). Known answer first: clean4 seed 0 through the wrapper must
# equal 435's stats bit for bit; abort otherwise. Then arcfull seeds 0-2, arcbank seed 0, and the paired summary.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/462-broad16-lab-signature; LOGD=$X/work/logs; mkdir -p "$LOGD"; LOG=$X/work/drive462.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/drive462.pid
log "armed (pid $$): waiting for 458e's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|461 never.*|builder .* stopping|artifact sha .* stopping)$" experiments/458-arc-lab-programme/work/chain458e.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain458e\.sh|chain461\.sh|chain458d\.sh|run_arm_generic\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609272000 ] && { log "458e never ended by 20:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local p=$1 s=$2; log "START $p seed $s"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run462.py "$p" base "$s" $ARGS > "$LOGD/$p-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $p seed $s peakRSS=$(grep 'maximum resident' "$LOGD/$p-s$s.log" | awk '{printf "%.2f GB", $1/1e9}')"; return $rc; }
run clean4 0 || { log "KA run failed"; exit 1; }
uv run python - <<'PY' >> "$LOG"
import numpy as np
a=np.load("experiments/462-broad16-lab-signature/work/h1_clean4/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
b=np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA462l (clean4 through the wrapper == 435 seed 0, stats bit-equal):", "PASS" if ok else "FAIL")
open("experiments/462-broad16-lab-signature/work/ka462l.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
/usr/bin/grep -q PASS $X/work/ka462l.txt || { log "KA FAIL -- stopping"; exit 1; }
for s in 0 1 2; do run arcfull $s || exit 1; done
run arcbank 0 || exit 1
uv run python $X/summary462l.py >> "$LOG" 2>&1
log "END"
