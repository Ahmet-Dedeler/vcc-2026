#!/bin/bash
# 517: serial REAL emitter runs on H1 (run517.py; ~3 min, ~7.5 GB each), started only after the production build queue:
# waits for >= 2026-09-29T02:00Z AND no builder/rehearsal process AND AGENTS' memory gate, before EVERY run.
# Known answers first; the queue aborts if one fails. Stops launching new runs after 2026-09-29T22:00Z.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/517-flex-programme"
cd "$REPO" || exit 1
LOGD="$HERE/work/logs"; mkdir -p "$LOGD"; LOG="$HERE/work/drive517.log"
echo $$ > "$HERE/work/drive517.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
busy(){ ps -axo comm,args | awk '$1 ~ /^(bash|python[0-9.]*|Python)$/ && /run_arm_generic\.sh|def_day\.sh|run_ka_4[0-9]+\.sh|make_submission_pm_cd4\.py|def_pipeline\.py/' | /usr/bin/grep -q .; }
gate(){
  while true; do
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202609292200 ] && { log "past 22:00Z 09-29: not starting new runs"; log "END"; exit 0; }
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    if [ "$now" -ge 202609290200 ] && ! busy && [ "${fp:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]; then
      log "gate ok: free ${fp}% free+inactive ${fi} GiB"; return 0; fi
    sleep 120
  done
}
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){  # mode seed [extra args]
  local mode=$1 s=$2; shift 2; gate
  log "START $mode seed $s $*"
  /usr/bin/time -l uv run python "$HERE/run517.py" "$mode" base "$s" $ARGS "$@" > "$LOGD/$mode-s$s$(echo "$*" | tr -d ' -').log" 2>&1; local rc=$?
  log "END rc=$rc $mode seed $s $* peakRSS=$(/usr/bin/grep 'maximum resident' "$LOGD/$mode-s$s$(echo "$*" | tr -d ' -').log" | awk '{printf "%.2f GB", $1/1e9}')"
  return $rc
}
ka(){  # mode seed reference-stats
  run "$1" "$2" || { log "KA run failed: $1"; log "END"; exit 1; }
  uv run python - "$HERE/work/h1_$1/stats-base_s$2_bps_hookem_mde0.1cap0.1min0.3coarse.npz" "$3" <<'PY' >> "$LOG" 2>&1
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
ok = all(np.array_equal(a[k], b[k], equal_nan=True) for k in ("sums", "jk", "p", "q", "lfc"))
print(("KA PASS " if ok else "KA FAIL ") + sys.argv[1].split("/")[-2] + " == " + sys.argv[2])
sys.exit(0 if ok else 1)
PY
  [ $? -eq 0 ] || { log "KA FAILED: $1 -- stopping"; log "END"; exit 1; }
}
log "armed (pid $$)"
ka c4-none 0 "$REPO/experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
ka g33-none 0 "$REPO/experiments/467-strength-gated-programme/work/h1_gated33/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
ka c4-c4t0.05 0 "$REPO/experiments/496-uncalled-sign-hook/work/h1/stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.npz"
for s in 0 1 2; do run c4-mzg0.05 $s; done          # NULL/PARTIAL arm, better estimator: file:<meanz_ne artifact>:0.05:0.05:le1
for s in 0 1 2; do run c4-arcg0.05 $s; done          # the same with the production P_down programme (458d's artifact)
run g33-mzg0.05 0                                     # POS arm on 458g with the meanz_ne hook source
for s in 0 1 2; do run g33-arcg0.05 $s; done         # POS arm with 458d's programme as the hook source (programme:0.05)
run mz-none 0                                         # meanz_ne as the UNGATED block prior (analytic: no reach cost)
run c4-arcg0.05-r0.25 0                               # the hook at a null-like programme quality
run c4-arct0.05 0                                     # the ungated file form
run arc-none-q0.7 0                                   # 458d's reach contingency (protect_q 0.7)
run mzg-none 0                                        # meanz_ne as the gated block prior
run c4-arcg0.1 0                                      # tau 0.1
for s in 1 2; do run c4-arcg0.05-r0.25 $s; done       # more seeds at the null-like quality
run g33-arcg0.1 0                                     # tau 0.1 on the 458g base
run c4-arcg0.05 0 --zs dn_k250_zs                     # additivity with the incumbent's zs rider (446bj2z; 496z did this for CLEAN4)
run c4-mzg0.05 0 --zs dn_k250_zs
uv run python "$HERE/summ_real517.py" >> "$LOG" 2>&1
uv run python "$HERE/append_real517.py" > /dev/null 2>> "$LOG"
log "END"
