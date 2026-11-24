#!/usr/bin/env bash
# 249 -- the serial heavy sequence, each step memory-gated and alone:
#   1. ARM build   VCC_CD4_DEFLATE_K=5 build249.sh submissions/249-cd4svd5-f034.h5ad   (~22 min, 20 GiB)
#   2. postbuild   def_postbuild.py on the arm's meta (--allow-config-diff cd4_deflate_k)
#   3. vcc prep    --dry-run, then the package submissions/249-cd4svd5-f034.vcc          (~25 GiB)
#   4. KA build    flag unset -> submissions/249-ka-h1union-g028.h5ad must be 8df4dadd... (row 1);
#                  postbuild; the h5ad is deleted if it matches (a byte duplicate of 222's), kept if not.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run249.log"
ARM="submissions/249-cd4svd5-f034"; KA="submissions/249-ka-h1union-g028"
WANT_H5AD="8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6"; WANT_BYTES=12648807800
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
rss_of(){ awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$1"; }
wall_of(){ awk '/ real /{print $1}' "$1" | tail -1; }

log "ARMED pid=$$"
# ---------------------------------------------------------------- 1. the arm
wait_gate
log "ARM BUILD start"
VCC_CD4_DEFLATE_K=5 bash "$HERE/build249.sh" "$ARM.h5ad" > "$W/build-arm.log" 2>&1
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
# ---------------------------------------------------------------- 4. the known-answer build
wait_gate
log "KA BUILD start (flag unset)"
env -u VCC_CD4_DEFLATE_K bash "$HERE/build249.sh" "$KA.h5ad" > "$W/build-ka.log" 2>&1
RC=$?; log "KA BUILD rc=$RC wall $(wall_of "$W/build-ka.log") s peak RSS $(rss_of "$W/build-ka.log") GiB"
[ $RC -eq 0 ] || { log "KA BUILD FAILED; see $W/build-ka.log"; exit 1; }
KA_SHA="$(shasum -a 256 "$KA.h5ad" | cut -d' ' -f1)"; KA_B="$(stat -f%z "$KA.h5ad")"
if [ "$KA_SHA" = "$WANT_H5AD" ] && [ "$KA_B" = "$WANT_BYTES" ]; then MATCH=MATCH; else MATCH=MISMATCH; fi
log "KA h5ad sha256 $KA_SHA bytes $KA_B  (want $WANT_H5AD $WANT_BYTES)  $MATCH"
uv run python scripts/def_postbuild.py --meta "$KA.h5ad.meta.json" --out "$W/postbuild-ka.json" > "$W/postbuild-ka.log" 2>&1
log "KA postbuild rc=$? $(grep -c FAIL "$W/postbuild-ka.log") FAIL lines"
cp "$KA.h5ad.meta.json" "$W/"
if [ "$MATCH" = "MATCH" ]; then
  rm -f "$KA.h5ad"; log "DELETED $KA.h5ad (byte duplicate of submissions/222-h1union-g028.h5ad; meta kept in $W and submissions/)"
else
  log "KEPT $KA.h5ad for inspection (a MISMATCH is a finding)"
fi
printf '%s  %s\n%s  %s\n' "$ARM_SHA" "$ARM.h5ad" "${VCC_SHA:-<no vcc>}" "$ARM.vcc" > "$W/SHA256SUMS.txt"
log "DONE"
