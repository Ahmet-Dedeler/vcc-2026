#!/usr/bin/env bash
# 467: H1 price of Arc's programme gated to design-strong targets. Serial after 466. Record (clean4, seed 0) doubles as the
# known answer (must equal 435 seed 0 bit for bit); then gated33/50/67 x seeds 0-2; then the summary.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/467-strength-gated-programme; LOGD=$X/work/logs; mkdir -p "$LOGD"; LOG=$X/work/chain467.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain467.pid
log "armed (pid $$): waiting for 466's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|part .* failed -- stopping|.*stopping)$" experiments/466-lab-signature-cells/work/chain466.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain466\.sh|chain465\.sh|def_day\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609272000 ] && { log "466 never ended by 20:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local m=$1 s=$2; log "START $m seed $s"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run467.py "$m" base "$s" $ARGS > "$LOGD/$m-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $m seed $s peakRSS=$(grep 'maximum resident' "$LOGD/$m-s$s.log" | awk '{printf "%.2f GB", $1/1e9}') $(grep '^\[467\]' "$LOGD/$m-s$s.log" | tail -1 | cut -c1-200)"; return $rc; }
run record 0 || { log "record failed"; exit 1; }
uv run python - <<'PY' >> "$LOG"
import numpy as np, glob
a=np.load(glob.glob("experiments/467-strength-gated-programme/work/h1_record/stats-base_s0_*coarse.npz")[0])
b=np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA467 (record = clean4 through the gated wrapper == 435 seed 0, stats bit-equal):", "PASS" if ok else "FAIL")
open("experiments/467-strength-gated-programme/work/ka467.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
/usr/bin/grep -q PASS $X/work/ka467.txt || { log "KA FAIL -- stopping"; exit 1; }
for q in 33 50 67; do for s in 0 1 2; do run gated$q $s || exit 1; done; done
uv run python $X/summary467.py >> "$LOG" 2>&1
log "END"
