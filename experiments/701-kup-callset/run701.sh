#!/usr/bin/env bash
# 701 H1 pricing chain (science session [bf003b], 2026-10-14): K-UP sizes on 537's real production-emitter H1 harness.
#   nohup bash experiments/701-kup-callset/run701.sh > experiments/701-kup-callset/work/run701.nohup 2>&1 &
# v2 (10-15 04:5xZ): waits until R7a's access wrapper (552/work/r7a_access_1015b.sh) has EXITED and its log's last line is
# its terminal "END rc=" (coordinator: R7a strictly first). gate602 then refuses while any builder runs, so the runs fill
# builder-free gaps only. (v1 waited for 671's DN chain; stopped 10-14 23:18Z before any run.)
# Then 9 runs, serially, each through gate602 (8 GB label; measured peak 7.0-7.3 GB in 537):
#   run537.py cpm-n<cpm> base <seed>, cpm in 3000 2800 3400 (the arm first), seed 0 1 2; outputs land in
#   experiments/537-arm-pipeline/work/h1_cpm-n<cpm>/ (run537's fixed output dir; new mode dirs only).
# Then price701.py (KA vs 537's recorded K pair first). Terminal line: "<ts> END" on every path. Cap for starting a run:
# 10-15 12:30Z (the 10-16 K-UP build must not wait on this price).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/701-kup-callset; W=$X/work; mkdir -p "$W"; LOG=$W/run701.log
CAP=202610151230
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/run701.pid"
log "armed v2 (pid $$): waiting for R7a's wrapper to exit with its END rc= line; cap $CAP"
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
until ! alive experiments/552-r7-reserve/work/r7a_access_1015b.sh && tail -1 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log 2>/dev/null | /usr/bin/grep -qE "END rc="; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "R7a not finished by the cap: stopping"; log "END"; exit 1; }
  sleep 60
done
log "R7a exited: $(tail -2 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log | tr '\n' ' ' | cut -c1-300)"
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for cpm in 3000 2800 3400; do
  for s in 0 1 2; do
    M="cpm-n$cpm"
    if [ -s "experiments/537-arm-pipeline/work/h1_$M/metrics-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json" ]; then
      log "skip $M seed $s (metrics exist)"; continue; fi
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "cap reached before $M seed $s"; log "END"; exit 1; }
    log "START $M seed $s"
    bash experiments/602-council-gate/gate602.sh "701-$M-s$s" 8 -- uv run python experiments/537-arm-pipeline/run537.py "$M" base "$s" $ARGS \
      > "$W/run_${M}_s$s.log" 2>&1
    rc=$?
    log "END rc=$rc $M seed $s: $(tail -1 "$W/run_${M}_s$s.log" | cut -c1-200)"
    [ "$rc" = 0 ] || { log "run failed: stopping"; log "END"; exit 1; }
  done
done
log "pricing"
uv run python $X/price701.py > "$W/price701.out" 2>&1
rc=$?
log "price701 rc=$rc: $(/usr/bin/grep -E '^CAL|KA vs' "$W/price701.out" | tr '\n' ' ' | cut -c1-900)"
log "END"
