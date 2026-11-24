#!/usr/bin/env bash
# 251 known-answer rows 1 and 2 at staging level, each row ALONE, memory-gated
# (mandate: memory_pressure >= 40% free AND vm_stat free+inactive >= 30 GiB AND no other builder
#  AND no vcc prep AND the working-tree builder at this experiment's intended bytes).
#   pre     PRE-251 byte copy, VCC_CD4_DEFLATE_BASIS=V5_abc   -> abeb0749...  (the flag does not exist there)
#   off     working tree, flags unset                         -> abeb0749...  (222 stage-on028.json, 247/248/249)
#   k5      working tree, VCC_CD4_DEFLATE_K=5                 -> a28d5be2...  (249's staged arm; the re-fit path is untouched)
#   pinned  working tree, VCC_CD4_DEFLATE_BASIS=V5_abc        -> compared with a28d5be2... by compare251.py
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
LOGDIR="$HERE/work/stage"; mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage251.log"
WANT_OFF="abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315"
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
  env "$@" /usr/bin/time -l uv run python "$HERE/stage251.py" --builder "$builder" --tag "$tag" > "$LOGDIR/$tag.log" 2>&1
  local rc=$?
  local rss got; rss="$(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$LOGDIR/$tag.log")"
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  log "END   $tag  rc=$rc  peak RSS ${rss} GiB  array_sha256 ${got}  $(grep -o 'cd4.fraction [^ ]*' "$LOGDIR/$tag.log" | head -1)"
  sleep 5
}
log "ARMED pid=$$"
run_row pre    "$HERE/work/make_submission_pm_cd4.PRE.py" VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_BASIS=V5_abc
run_row off    scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
run_row k5     scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5
run_row pinned scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_BASIS=V5_abc
for tag in pre off; do
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  [ "$got" = "$WANT_OFF" ] && log "KA1(staging) $tag MATCH $WANT_OFF" || log "KA1(staging) $tag MISMATCH got ${got:-<none>} want $WANT_OFF"
done
got="$(grep -o 'array_sha256 .*' "$LOGDIR/k5.log" | head -1 | awk '{print $2}')"
[ "$got" = "$WANT_K5" ] && log "KA2(re-fit untouched) k5 MATCH $WANT_K5" || log "KA2(re-fit untouched) k5 MISMATCH got ${got:-<none>} want $WANT_K5"
uv run python "$HERE/compare251.py" 2>&1 | tee -a "$LOG"
log "ALL DONE"
