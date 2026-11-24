#!/usr/bin/env bash
# 466: the lab-signature test with cells, full run, serial after 465's rehearsal (each part ~18-20 GB peak; ~35 min).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/466-lab-signature-cells; W=$X/work; LOG=$W/chain466.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain466.pid"
log "armed (pid $$): waiting for 465's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|.*stopping)$" experiments/465-both-switches-efr/work/chain465.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain465\.sh|def_day\.sh|drive462\.sh|chain458e\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609272000 ] && { log "465 never ended by 20:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
for pt in HCT116 HEK293T K562 summary; do
  log "START $pt"
  /usr/bin/time -l uv run python $X/lab466.py --part $pt > "$W/part-$pt.log" 2>&1; rc=$?
  log "END $pt rc=$rc peakRSS=$(grep 'maximum resident' "$W/part-$pt.log" | awk '{printf "%.1f GB", $1/2^30}')"
  [ $rc -eq 0 ] || { log "part $pt failed -- stopping"; exit 1; }
done
log "END"
