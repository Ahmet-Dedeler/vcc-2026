#!/usr/bin/env bash
# 255 -- the arm's effect matrix (k = 10) and known-answer row 2 (k = 5 re-staged), each row ALONE,
# memory-gated (mandate: memory_pressure >= 40% free AND vm_stat free+inactive >= 30 GiB AND no other
# builder AND no vcc prep AND the working-tree builder at git HEAD's bytes).
#   svd10  working tree, VCC_CD4_DEFLATE_K=10  -> the arm's staged effect matrix
#   svd5   working tree, VCC_CD4_DEFLATE_K=5   -> a28d5be2...  (249's staged arm = cd4svd5's effects.final_sha256)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
LOGDIR="$HERE/work/stage"; mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage255.log"
WANT_K5="a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others
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
  env -u VCC_CD4_DEFLATE_K -u VCC_CD4_PARTIAL_G -u VCC_CD4_DEFLATE_BASIS -u VCC_UNION2_SOURCES -u VCC_UNION2_G -u VCC_CD4_SIZING -u VCC_CD4_TARGET_G -u VCC_UNION_SOURCES "$@" /usr/bin/time -l uv run python "$HERE/stage255.py" --builder "$builder" --tag "$tag" > "$LOGDIR/$tag.log" 2>&1
  local rc=$?
  local rss got; rss="$(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$LOGDIR/$tag.log")"
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  log "END   $tag  rc=$rc  peak RSS ${rss} GiB  array_sha256 ${got}  $(grep -o 'cd4.fraction [^ ]*' "$LOGDIR/$tag.log" | head -1)"
  sleep 5
}
log "ARMED pid=$$"
run_row svd10 scripts/make_submission_pm_cd4.py VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=10
run_row svd5  scripts/make_submission_pm_cd4.py VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5
got="$(grep -o 'array_sha256 .*' "$LOGDIR/svd5.log" | head -1 | awk '{print $2}')"
[ "$got" = "$WANT_K5" ] && log "KA2(staging) svd5 MATCH $WANT_K5" || log "KA2(staging) svd5 MISMATCH got ${got:-<none>} want $WANT_K5"
log "ALL DONE"
