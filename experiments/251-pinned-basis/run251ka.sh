#!/usr/bin/env bash
# 251 -- known-answer row 1 at PACKAGE level, gated and alone (248's `run248ka.sh` pattern):
#   both flags unset through the changed builder -> work/build/251-ka-h1union-g028.h5ad must be
#   8df4dadd... on 12,648,807,800 bytes (the shipped h1union); def_postbuild against the 222 meta;
#   the h5ad is deleted on MATCH (a byte duplicate of submissions/222-h1union-g028.h5ad), kept if not.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run251ka.log"
KA="$W/251-ka-h1union-g028"
WANT_H5AD="8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6"; WANT_BYTES=12648807800
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
rss_of(){ awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$1"; }
wall_of(){ awk '/ real /{print $1}' "$1" | tail -1; }

log "ARMED pid=$$"
wait_gate
log "KA BUILD start (both flags unset)"
env -u VCC_CD4_DEFLATE_K -u VCC_CD4_DEFLATE_BASIS bash "$HERE/build251.sh" "$KA.h5ad" > "$W/build-ka.log" 2>&1
RC=$?; log "KA BUILD rc=$RC wall $(wall_of "$W/build-ka.log") s peak RSS $(rss_of "$W/build-ka.log") GiB"
[ $RC -eq 0 ] || { log "KA BUILD FAILED; see $W/build-ka.log"; exit 1; }
KA_SHA="$(shasum -a 256 "$KA.h5ad" | cut -d' ' -f1)"; KA_B="$(stat -f%z "$KA.h5ad")"
if [ "$KA_SHA" = "$WANT_H5AD" ] && [ "$KA_B" = "$WANT_BYTES" ]; then MATCH=MATCH; else MATCH=MISMATCH; fi
log "KA h5ad sha256 $KA_SHA bytes $KA_B  (want $WANT_H5AD $WANT_BYTES)  $MATCH"
uv run python scripts/def_postbuild.py --meta "$KA.h5ad.meta.json" --out "$W/postbuild-ka.json" > "$W/postbuild-ka.log" 2>&1
log "KA postbuild rc=$? $(grep -E 'checks pass' "$W/postbuild-ka.log") $(grep -c FAIL "$W/postbuild-ka.log") FAIL lines"
if [ "$MATCH" = "MATCH" ]; then
  rm -f "$KA.h5ad"; log "DELETED $KA.h5ad (byte duplicate of submissions/222-h1union-g028.h5ad; meta kept in $W)"
else
  log "KEPT $KA.h5ad for inspection (a MISMATCH is a finding)"
fi
log "DONE"
