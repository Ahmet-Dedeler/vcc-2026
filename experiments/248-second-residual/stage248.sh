#!/usr/bin/env bash
# 248 known-answer rows 1 and 2 at staging level, plus the two arms' effect matrices.  Each row
# ALONE and memory-gated (mandate: memory_pressure >= 40% free AND vm_stat free+inactive >= 30 GiB
# AND no other builder AND no vcc prep), the working-tree builder pinned to this experiment's bytes.
#   off     working tree, union 2 unset                        -> abeb0749...  (KA row 1, 222 stage-on028.json)
#   pre_u2  PRE-248 byte copy, VCC_UNION2_* SET                -> abeb0749...  (the flag does not exist there)
#   swap    working tree, union 1 UNSET, union 2 = the H1 pair -> abeb0749...  (KA row 2, order invariance)
#   cum1    working tree, union 1 = H1 pair @0.28, union 2 = dld1_lowmoi @0.35          (the rider)
#   cum4    working tree, union 1 = H1 pair @0.28, union 2 = 225's four @0.35           (CUM4)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
LOGDIR="$HERE/work/stage"; mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage248.log"
WANT="abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others preps want have bstate
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  want="$(cat "$HERE/work/builder.intended.sha256" 2>/dev/null)"
  have="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  [ "$have" = "$want" ] && bstate="intended" || bstate="NOT-INTENDED(${have:0:8})"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  builder ${bstate}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && [ "$others" = "0" ] && [ "$preps" = "0" ] && [ "$bstate" = "intended" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened; giving up"; exit 3; }; sleep 30; done; }
run_row(){
  local tag="$1" builder="$2"; shift 2
  wait_gate
  log "START $tag  builder=$builder  env: $*"
  env -u VCC_UNION2_SOURCES -u VCC_UNION2_G -u VCC_H1_UNION_G -u VCC_UNION_SOURCES -u VCC_CD4_SIZING -u VCC_CD4_TARGET_G -u VCC_CD4_DEFLATE_K \
    "$@" /usr/bin/time -l uv run python "$HERE/stage248.py" --builder "$builder" --tag "$tag" > "$LOGDIR/$tag.log" 2>&1
  local rc=$?
  local rss got; rss="$(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$LOGDIR/$tag.log")"
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  log "END   $tag  rc=$rc  peak RSS ${rss} GiB  array_sha256 ${got}  $(grep -o 'cd4.fraction [^ ]*' "$LOGDIR/$tag.log" | head -1)  $(grep -o 'union2.solved_fraction .*' "$LOGDIR/$tag.log" | head -1)"
  sleep 5
}
log "ARMED pid=$$"
run_row off    scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
run_row pre_u2 "$HERE/work/make_submission_pm_cd4.PRE.py" VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35
run_row swap   scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_UNION2_SOURCES=gse295214_h1,h1_2025 VCC_UNION2_G=0.28
run_row cum1   scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35
run_row cum4   scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 VCC_UNION2_G=0.35
for tag in off pre_u2 swap; do
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  [ "$got" = "$WANT" ] && log "KA(staging) $tag MATCH $WANT" || log "KA(staging) $tag MISMATCH got ${got:-<none>} want $WANT"
done
log "ALL DONE"
