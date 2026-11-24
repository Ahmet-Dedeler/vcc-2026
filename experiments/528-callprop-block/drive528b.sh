#!/bin/bash
# 528b: interaction checks, after drive528.sh's queue: the propensity block (lcap) on top of (1) 517's uncalled-gene hook as
# 527 builds it (446bj2zmz), and (2) 523's smaller call set (CallSet(400, 1200): K 1,600). Same gates as drive528.sh.
# KA first: `cpm+mzg` seed 0 == 517's real c4-mzg0.05 seed 0 (pre-zs and zs stats bit-equal); the hook rows abort if it fails.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/528-callprop-block"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive528b.log"
echo $$ > "$HERE/work/drive528b.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){  # $1 = pattern on the executable (first word of args), $2 = pattern on the args
  ps -axo args= | awk -v c="$1" -v a="$2" '$1 ~ c && $0 ~ a {n++} END {exit n ? 0 : 1}'
}
log "armed (pid $$): waiting for drive528.sh to exit"
while alive "bash$" "528-callprop-block/drive528[.]sh"; do sleep 60; done
log "drive528.sh exited: $(tail -1 "$HERE/work/drive528.log")"
gate(){
  local waited=0
  while true; do
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202609300600 ] && { log "past 06:00Z 09-30: not starting new runs"; log "END"; exit 0; }
    why=""
    alive "bash$" "527-arc-hook-arm/build527[.]sh" && why="$why build527"
    alive "python[0-9.]*$" "517-flex-programme/run517[.]py" && why="$why run517"
    alive "python[0-9.]*$" "make_submission_pm_cd4[.]py|def_pipeline[.]py" && why="$why builder"
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] || why="$why free${fp}%"
    [ "${fi:-0}" -ge 23 ] || why="$why fi${fi}GiB"
    if [ -z "$why" ]; then log "gate ok: free ${fp}% free+inactive ${fi} GiB"; return 0; fi
    if [ $((waited % 1800)) -eq 0 ]; then log "waiting:$why"; fi
    sleep 60; waited=$((waited + 60))
  done
}
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
run(){  # mode seed
  local mode=$1 s=$2; gate
  log "START $mode seed $s"
  /usr/bin/time -l uv run python "$HERE/run528.py" "$mode" base "$s" $ARGS > "$LOGD/$mode-s$s.log" 2>&1; local rc=$?
  log "END rc=$rc $mode seed $s peakRSS=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-s$s.log" | awk '{printf "%.2f GB", $1/1e9}')"
  return $rc
}
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
[ "$HOOK_OK" = "1" ] && run "lcap+mzg" 0
run "cpm-n1200" 0
run "lcap-n1200" 0
[ "$HOOK_OK" = "1" ] && { run "cpm+mzg" 1; run "lcap+mzg" 1; }
run "cpm-n1200" 1
run "lcap-n1200" 1
log "QUEUE DONE"
log "END"
