#!/usr/bin/env bash
# 249 known-answer row 1 (staging level) + the arm's effect matrix, each row ALONE, memory-gated
# (mandate: memory_pressure >= 40% free AND vm_stat free+inactive >= 30 GiB AND no other builder).
#   pre5   PRE-249 byte copy, VCC_CD4_DEFLATE_K=5     -> abeb0749...  (the flag does not exist there)
#   off    working tree, flag unset                    -> abeb0749...  (222 stage-on028.json, re-verified by 247)
#   svd5   working tree, VCC_CD4_DEFLATE_K=5           -> the arm
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
LOGDIR="$HERE/work/stage"; mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage249.log"
WANT="abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  # `other_builders()` sees builds and operator fits; a `vcc prep` (22-25 GiB) is not a builder, so
  # count those separately -- the mandate forbids running alongside a build OR a prep.
  # (an agent's own command line quotes "vcc prep" in its mandate text, so `claude` lines are
  # excluded and the token must start a word)
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  # and the working-tree builder must be THIS experiment's intended bytes: 244's rehearsal pins
  # the registry into the same file mid-run (def_pin_panel_artifacts --apply) and restores it after
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
  env "$@" /usr/bin/time -l uv run python "$HERE/stage249.py" --builder "$builder" --tag "$tag" > "$LOGDIR/$tag.log" 2>&1
  local rc=$?
  local rss got; rss="$(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$LOGDIR/$tag.log")"
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  log "END   $tag  rc=$rc  peak RSS ${rss} GiB  array_sha256 ${got}  $(grep -o 'cd4.fraction [^ ]*' "$LOGDIR/$tag.log" | head -1)"
  sleep 5
}
log "ARMED pid=$$"
run_row pre5 "$HERE/work/make_submission_pm_cd4.PRE.py" VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5
run_row off  scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
run_row svd5 scripts/make_submission_pm_cd4.py           VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5
for tag in pre5 off; do
  got="$(grep -o 'array_sha256 .*' "$LOGDIR/$tag.log" | head -1 | awk '{print $2}')"
  [ "$got" = "$WANT" ] && log "KA1(staging) $tag MATCH $WANT" || log "KA1(staging) $tag MISMATCH got ${got:-<none>} want $WANT"
done
log "ALL DONE"
