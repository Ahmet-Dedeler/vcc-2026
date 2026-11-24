#!/bin/bash
# 528: serial REAL production-emitter runs on H1 (run528.py; ~5.5 min, ~7.5 GB each, with the zs pass).
# Before EVERY run: 527's build chain has exited, no 517 H1 job and no production builder is alive, and AGENTS' memory gate
# (memory_pressure free >= 40 %, vm_stat free+inactive >= 23 GiB) holds. Stops launching new runs after 2026-09-30T06:00Z.
# Known answer first (cpm seed 0 == 435's production run, pre-zs and zs stats bit-equal); the queue aborts if it fails.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/528-callprop-block"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive528.log"
echo $$ > "$HERE/work/drive528.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){  # $1 = pattern on the executable (first word of args; `comm` is truncated on macOS), $2 = pattern on the args
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
log "armed (pid $$)"
run cpm 0 || { log "KA run failed"; log "END"; exit 1; }
PYTHONDONTWRITEBYTECODE=1 uv run python - "$HERE/work/h1_cpm" "$REPO/experiments/435-zs-on-coarse/work/h1" <<'PY' >> "$LOG" 2>&1
import sys, numpy as np
a, b = sys.argv[1], sys.argv[2]
ok = True
for tag in ("base_s0_bps_prodem_mde0.1cap0.1min0.3coarse", "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs"):
    x, y = np.load(f"{a}/stats-{tag}.npz"), np.load(f"{b}/stats-{tag}.npz")
    eq = all(np.array_equal(x[k], y[k], equal_nan=True) for k in ("sums", "jk", "p", "q", "lfc"))
    print(("KA PASS " if eq else "KA FAIL ") + tag)
    ok &= eq
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || { log "KA FAILED -- stopping"; log "END"; exit 1; }
for s in 0 1 2; do
  run lcap "$s"
  run shuf "$s"
  run lcap-k1980 "$s"
  run lcap-k1650 "$s"
  run ltv "$s"
done
run lcap-k1100 0
run p2f 0
log "QUEUE DONE"
log "END"
