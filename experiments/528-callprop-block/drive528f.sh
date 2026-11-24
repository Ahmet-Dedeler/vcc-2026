#!/bin/bash
# 528f: after drive528e.sh's queue: the FINAL emitter (k_emit528v3.py = patch_emitter528.py's output, the one the builder
# patch pins) re-runs `lcap` (steer key; the queue ran it on v1 = k_emit528.py) and `lcap-dec` (steer cpm; ran on v2) at
# seed 0 into h1_<mode>-v3check; both must equal the v1 / v2 runs bit for bit (pre-zs and zs stats). Same gates as 528e.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/528-callprop-block"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive528f.log"
echo $$ > "$HERE/work/drive528f.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo args= | awk -v c="$1" -v a="$2" '$1 ~ c && $0 ~ a {n++} END {exit n ? 0 : 1}'; }
log "armed (pid $$): waiting for drive528e.sh to exit"
while alive "bash$" "drive528e[.]sh"; do sleep 60; done
log "drive528e.sh exited: $(tail -1 "$HERE/work/drive528e.log")"
gate(){
  while true; do
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202609300600 ] && { log "past 06:00Z 09-30"; log "END"; exit 0; }
    why=""
    alive "python[0-9.]*$" "517-flex-programme/run517[.]py" && why="$why run517"
    alive "python[0-9.]*$" "528-callprop-block/run528[.]py" && why="$why run528"
    need=23
    alive "python[0-9.]*$" "make_submission_pm_cd4[.]py|def_pipeline[.]py" && need=25
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    [ "${fp:-0}" -ge 40 ] || why="$why free${fp}%"
    [ "${fi:-0}" -ge "$need" ] || why="$why fi${fi}<${need}GiB"
    [ -z "$why" ] && { log "gate ok: free ${fp}% free+inactive ${fi} GiB (need ${need})"; return 0; }
    sleep 60
  done
}
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for mode in lcap lcap-dec; do
  gate
  log "START v3check $mode seed 0"
  RUN528_EMITTER=v3 /usr/bin/time -l uv run python "$HERE/run528.py" "$mode" base 0 $ARGS > "$LOGD/$mode-v3check-s0.log" 2>&1
  log "END rc=$? v3check $mode seed 0 peakRSS=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-v3check-s0.log" | awk '{printf "%.2f GB", $1/1e9}')"
  PYTHONDONTWRITEBYTECODE=1 uv run python - "$HERE/work/h1_$mode-v3check" "$HERE/work/h1_$mode" <<'PY' >> "$LOG" 2>&1
import sys, numpy as np
a, b = sys.argv[1], sys.argv[2]
ok = True
for tag in ("base_s0_bps_prodem_mde0.1cap0.1min0.3coarse", "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs"):
    x, y = np.load(f"{a}/stats-{tag}.npz"), np.load(f"{b}/stats-{tag}.npz")
    eq = all(np.array_equal(x[k], y[k], equal_nan=True) for k in ("sums", "jk", "p", "q", "lfc"))
    print(("V3CHECK PASS " if eq else "V3CHECK FAIL ") + b.rsplit("/", 1)[-1] + " " + tag)
    ok &= eq
sys.exit(0 if ok else 1)
PY
done
log "END"
