#!/bin/bash
# 528d: replaces drive528c.sh (killed 13:0xZ between runs; its running child finished on its own). Same gates, plus: no other
# run528.py of mine alive (one heavy step of mine at a time), and the skip check runs AFTER the gate, so an orphaned run's
# outputs are seen. Queue: the main arms (skipping finished runs) with the decoupled arm `lcap-dec` (block_steer="cpm": only
# the call set follows the key) added per seed, then extras, then the interaction checks (517's hook as 527 builds it;
# 523's K 1,600) for both lcap forms. Stops launching new runs after 2026-09-30T06:00Z.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/528-callprop-block"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive528d.log"
echo $$ > "$HERE/work/drive528d.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){  # $1 = pattern on the executable (first word of args), $2 = pattern on the args
  ps -axo args= | awk -v c="$1" -v a="$2" '$1 ~ c && $0 ~ a {n++} END {exit n ? 0 : 1}'
}
gate(){
  local waited=0
  while true; do
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202609300600 ] && { log "past 06:00Z 09-30: not starting new runs"; log "END"; exit 0; }
    why=""
    alive "bash$" "527-arc-hook-arm/build527[.]sh" && why="$why build527"
    alive "python[0-9.]*$" "517-flex-programme/run517[.]py" && why="$why run517"
    alive "python[0-9.]*$" "528-callprop-block/run528[.]py" && why="$why run528"
    need=23
    alive "python[0-9.]*$" "make_submission_pm_cd4[.]py|def_pipeline[.]py" && need=25
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] || why="$why free${fp}%"
    [ "${fi:-0}" -ge "$need" ] || why="$why fi${fi}<${need}GiB"
    if [ -z "$why" ]; then log "gate ok: free ${fp}% free+inactive ${fi} GiB (need ${need})"; return 0; fi
    if [ $((waited % 1800)) -eq 0 ]; then log "waiting:$why"; fi
    sleep 60; waited=$((waited + 60))
  done
}
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
run(){  # mode seed
  local mode=$1 s=$2
  local m="$HERE/work/h1_$mode/metrics-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"
  gate
  [ -s "$m" ] && { log "SKIP $mode seed $s (metrics exist)"; return 0; }
  log "START $mode seed $s"
  /usr/bin/time -l uv run python "$HERE/run528.py" "$mode" base "$s" $ARGS > "$LOGD/$mode-s$s.log" 2>&1; local rc=$?
  log "END rc=$rc $mode seed $s peakRSS=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-s$s.log" | awk '{printf "%.2f GB", $1/1e9}')"
  return $rc
}
log "armed (pid $$): drive528's queue after its KA PASS (+ the decoupled arm), then the interactions"
/usr/bin/grep -q "KA PASS base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs" "$HERE/work/drive528.log" \
  || { log "drive528's KA did not pass: stopping"; log "END"; exit 1; }
for s in 0 1 2; do
  run lcap "$s"
  run shuf "$s"
  run lcap-dec "$s"
  run lcap-k1980 "$s"
  run lcap-k1650 "$s"
  run ltv "$s"
done
run lcap-k1650-dec 0
run lcap-k1100 0
run p2f 0
log "MAIN QUEUE DONE"
HOOK_OK=0
if run "cpm+mzg" 0; then
  PYTHONDONTWRITEBYTECODE=1 uv run python - "$HERE/work/h1_cpm+mzg" "$REPO/experiments/517-flex-programme/work" <<'PY' >> "$LOG" 2>&1
import sys, numpy as np
a, b = sys.argv[1], sys.argv[2]
ok = True
for mine, theirs in (("base_s0_bps_prodem_mde0.1cap0.1min0.3coarse", "h1_c4-mzg0.05/stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse"),
                     ("base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs",
                      "h1_c4-mzg0.05_zs/stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse_dnk250zs")):
    x, y = np.load(f"{a}/stats-{mine}.npz"), np.load(f"{b}/{theirs}.npz")
    eq = all(np.array_equal(x[k], y[k], equal_nan=True) for k in ("sums", "jk", "p", "q", "lfc"))
    print(("KA-hook PASS " if eq else "KA-hook FAIL ") + mine)
    ok &= eq
sys.exit(0 if ok else 1)
PY
  [ $? -eq 0 ] && HOOK_OK=1 || log "KA-hook FAILED: the hook rows are skipped"
fi
[ "$HOOK_OK" = "1" ] && { run "lcap+mzg" 0; run "lcap-dec+mzg" 0; }
run "cpm-n1200" 0
run "lcap-n1200" 0
run "lcap-dec-n1200" 0
log "QUEUE DONE"
log "END"
