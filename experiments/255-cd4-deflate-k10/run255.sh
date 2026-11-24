#!/usr/bin/env bash
# 255 -- the serial heavy sequence, each step memory-gated and alone:
#   1. ARM build   VCC_CD4_DEFLATE_K=10 build255.sh submissions/255-cd4svd10-f034.h5ad   (~22 min, 20 GiB)
#   2. postbuild   def_postbuild.py on the arm's meta (--allow-config-diff cd4_deflate_k; the reference
#                  stays the shipped 222 meta, as 249's and 254's did)
#   3. vcc prep    --dry-run, then the package submissions/255-cd4svd10-f034.vcc          (~25 GiB)
#   NO known-answer rebuild of the incumbent: 254 section KA 1 rebuilt 8df4dadd... through these builder
#   bytes (3615d20f...) at 20:36-21:02Z on 2026-09-04; the mandate says cite it, do not repeat it.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run255.log"
ARM="submissions/255-cd4svd10-f034"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  # `other_builders()` sees builds and operator fits; a `vcc prep` (22-25 GiB) is not a builder, so
  # count those separately -- the mandate forbids running alongside a build OR a prep.
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  # and the working-tree builder must be git HEAD's bytes (this experiment edits no production file)
  want="$(cat "$HERE/work/builder.intended.sha256" 2>/dev/null)"
  have="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  [ "$have" = "$want" ] && bstate="intended" || bstate="NOT-INTENDED(${have:0:8})"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  builder ${bstate}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && [ "$others" = "0" ] && [ "$preps" = "0" ] && [ "$bstate" = "intended" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened; giving up"; exit 3; }; sleep 30; done; }
rss_of(){ awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$1"; }
wall_of(){ awk '/ real /{print $1}' "$1" | tail -1; }

log "ARMED pid=$$"
# ---------------------------------------------------------------- 1. the arm
wait_gate
log "ARM BUILD start (VCC_CD4_DEFLATE_K=10)"
VCC_CD4_DEFLATE_K=10 bash "$HERE/build255.sh" "$ARM.h5ad" > "$W/build-arm.log" 2>&1
RC=$?; log "ARM BUILD rc=$RC wall $(wall_of "$W/build-arm.log") s peak RSS $(rss_of "$W/build-arm.log") GiB"
[ $RC -eq 0 ] || { log "ARM BUILD FAILED; see $W/build-arm.log"; exit 1; }
ARM_SHA="$(shasum -a 256 "$ARM.h5ad" | cut -d' ' -f1)"; ARM_B="$(stat -f%z "$ARM.h5ad")"
log "ARM h5ad sha256 $ARM_SHA bytes $ARM_B"
# ---------------------------------------------------------------- 2. postbuild
uv run python scripts/def_postbuild.py --meta "$ARM.h5ad.meta.json" --allow-config-diff cd4_deflate_k --out "$W/postbuild-arm.json" > "$W/postbuild-arm.log" 2>&1
log "ARM postbuild rc=$? $(grep -cE '^ *PASS|PASS ' "$W/postbuild-arm.log") PASS lines, $(grep -c FAIL "$W/postbuild-arm.log") FAIL lines"
# ---------------------------------------------------------------- 3. prep
wait_gate
log "PREP dry-run start"
/usr/bin/time -l vcc prep "$ARM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$ARM.vcc" --dry-run > "$W/prep-dryrun.log" 2>&1
log "PREP dry-run rc=$? peak RSS $(rss_of "$W/prep-dryrun.log") GiB"
wait_gate
log "PREP start"
/usr/bin/time -l vcc prep "$ARM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$ARM.vcc" > "$W/prep.log" 2>&1
RC=$?; log "PREP rc=$RC wall $(wall_of "$W/prep.log") s peak RSS $(rss_of "$W/prep.log") GiB"
if [ -f "$ARM.vcc" ]; then
  VCC_SHA="$(shasum -a 256 "$ARM.vcc" | cut -d' ' -f1)"; VCC_B="$(stat -f%z "$ARM.vcc")"
  log "ARM vcc sha256 $VCC_SHA bytes $VCC_B"
fi
printf '%s  %s\n%s  %s\n' "$ARM_SHA" "$ARM.h5ad" "${VCC_SHA:-<no vcc>}" "$ARM.vcc" > "$W/SHA256SUMS.txt"
log "DONE"
