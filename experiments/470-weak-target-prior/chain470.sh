#!/usr/bin/env bash
# 470: which block prior do design-weak targets want (H1)? Serial after chain469 (a def_day build). wc4 seed 0 is the
# known answer (== 467's gated33 seed 0, stats bit-equal); then wdown / warcw / warcs / sarcs x seeds 0-2; the summary.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/470-weak-target-prior; LOGD=$X/work/logs; LOG=$X/work/chain470.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain470.pid
log "armed (pid $$): waiting for chain469 to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|def_day did not finish.*|relabel failed|prep failed|prep dry run failed|.*stopping)$" \
        experiments/469-gated-october-path/work/chain469.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain469\.sh|chain468\.sh|def_day\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "469 never ended by 12:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local m=$1 s=$2; log "START $m seed $s"; RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run470.py "$m" base "$s" $ARGS > "$LOGD/$m-s$s.log" 2>&1; local rc=$?
       log "END rc=$rc $m seed $s peakRSS=$(grep 'maximum resident' "$LOGD/$m-s$s.log" | awk '{printf "%.2f GB", $1/1e9}') $(grep '^\[470\] '"$m" "$LOGD/$m-s$s.log" | tail -1 | cut -c1-160)"; return $rc; }
run wc4 0 || { log "wc4 failed"; exit 1; }
uv run python - <<'PY' >> "$LOG"
import numpy as np, glob
a=np.load(glob.glob("experiments/470-weak-target-prior/work/h1_wc4/stats-base_s0_*coarse.npz")[0])
b=np.load(glob.glob("experiments/467-strength-gated-programme/work/h1_gated33/stats-base_s0_*coarse.npz")[0])
ok=all(np.array_equal(a[k],b[k],equal_nan=True) for k in ("sums","jk","p","q","lfc"))
print("KA470 (wc4 == 467 gated33 seed 0, stats bit-equal):", "PASS" if ok else "FAIL")
open("experiments/470-weak-target-prior/work/ka470.txt","w").write("PASS\n" if ok else "FAIL\n")
PY
/usr/bin/grep -q PASS $X/work/ka470.txt || { log "KA FAIL -- stopping"; exit 1; }
for m in wdown warcw warcs sarcs; do for s in 0 1 2; do run $m $s || exit 1; done; done
uv run python $X/summary470.py >> "$LOG" 2>&1
log "END"
