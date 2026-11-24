#!/usr/bin/env bash
# 247 item C — the three staging rows, each ALONE, gated on memory (mandate: >= 15 GiB waits for
# memory_pressure >= 40% free AND vm_stat free+inactive >= 30 GiB, and no other builder).
#   off     VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34                     -> must reproduce 87a4b6d1... (222 stage-off)
#   on028   + VCC_H1_UNION_G=0.28                                      -> must reproduce abeb0749... (222 stage-on028)
#   carry   + VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641 -> cd4.fraction == 0.32207201067578145 to <= 1e-9
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
LOGDIR="$HERE/work/stage"; mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage247.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }

gate_ok(){
  local pct fi others
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB  other builders ${others}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && [ "$others" = "0" ]
}

wait_gate(){
  local n=0
  until gate_ok; do
    n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened after $n polls (30 s each); giving up"; exit 3; }
    sleep 30
  done
}

run_row(){
  local tag="$1"; shift
  wait_gate
  log "START $tag  env: $*"
  env "$@" /usr/bin/time -l uv run python "$HERE/stage247.py" --tag "$tag" > "$LOGDIR/$tag.log" 2>&1
  local rc=$?
  local rss; rss="$(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$LOGDIR/$tag.log")"
  log "END   $tag  rc=$rc  peak RSS ${rss} GiB  $(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1)  $(grep -o 'cd4.fraction .*' "$LOGDIR/$tag.log" | head -1)"
  sleep 5
}

log "ARMED pid=$$"
run_row off   VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34
run_row on028 VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
run_row carry VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641
log "ALL DONE"
