#!/usr/bin/env bash
# 256 -- KNOWN-ANSWER ROW 1 (240's row 1, re-run against the NEW incumbent), gated and alone:
#   the driver's DEFAULTS on data/challenge (both deflation variables unset in the environment) must
#   rebuild the shipped `cd4svd5` byte-for-byte:
#     work/ka1/256-ka-cd4svd5-f034.h5ad  -> aba27758... on 12,648,807,800 B
#     vcc prep of it                     -> 6e8b41c6... (the scored package, entry Xv0ayvuoHpsKYWv2iu42)
#   def_postbuild runs inside the driver (--stages build,postbuild) against the 249 meta (the new default).
#   Both artifacts are deleted on MATCH (byte duplicates of submissions/249-cd4svd5-f034.*), kept if not.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
source .venv/bin/activate 2>/dev/null || true
W="$HERE/work/ka1"; mkdir -p "$W"; LOG="$W/run256ka1.log"
NAME="256-ka-cd4svd5-f034"; KA="$W/$NAME"
WANT_H5AD="aba2775836df3bcf9619a520809fc7a71f325c0a81a881e6f9eb16c970f65faa"; WANT_BYTES=12648807800
WANT_VCC="6e8b41c6244b32dbe2a9765c7d2f2293b2032c49e59b455982b95bc6bb074533"
MIN_FI=23
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others preps
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB (need >= ${MIN_FI})  other builders ${others}  vcc preps ${preps}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" -v m="$MIN_FI" 'BEGIN{exit !(v>=m)}' && [ "$others" = "0" ] && [ "$preps" = "0" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened in 2 h; giving up"; exit 3; }; sleep 30; done; }
rss_of(){ awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$1" | tail -1; }
wall_of(){ awk '/ real /{print $1}' "$1" | tail -1; }

log "ARMED pid=$$  driver $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1 | cut -c1-16)  postbuild $(shasum -a 256 scripts/def_postbuild.py | cut -d' ' -f1 | cut -c1-16)  builder $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1 | cut -c1-16)  git $(git rev-parse --short HEAD)"
log "env: VCC_CD4_DEFLATE_K/BASIS and every VCC_* arm UNSET -- the driver's defaults supply everything"
wait_gate
log "DRIVER start: defaults on data/challenge, --stages build,postbuild"
env -u VCC_CD4_DEFLATE_K -u VCC_CD4_DEFLATE_BASIS -u VCC_UNION_SOURCES -u VCC_UNION2_SOURCES -u VCC_UNION2_G \
    -u VCC_CD4_SIZING -u VCC_CD4_TARGET_G -u VCC_ALPHA -u VCC_CD4_TARGET_F -u VCC_H1_UNION_G -u VCC_MASS_FIX -u VCC_CD4_PARTIAL_G \
    /usr/bin/time -l uv run python scripts/def_pipeline.py --bundle data/challenge \
    --work "$W" --name "$NAME" --stages build,postbuild \
    --whitener-out data/donors/whitener_union_pruned.npz --whitener-id union_pruned \
    --cd4-out data/external/cd4t/cloud-compact/cd4_compact --min-free-gib "$MIN_FI" \
    > "$W/driver.log" 2>&1
RC=$?; log "DRIVER rc=$RC wall $(wall_of "$W/driver.log") s peak RSS $(rss_of "$W/driver.log") GiB"
grep -E '^\[arms\]|^  (OK|NO ) config\.' "$W/driver.log" | tee -a "$LOG" >/dev/null
[ $RC -eq 0 ] || { log "DRIVER FAILED (rc=$RC); see $W/driver.log -- KA ROW 1 NOT ESTABLISHED"; exit 1; }
[ -f "$KA.h5ad" ] || { log "no artifact at $KA.h5ad -- KA ROW 1 NOT ESTABLISHED"; exit 1; }
KA_SHA="$(shasum -a 256 "$KA.h5ad" | cut -d' ' -f1)"; KA_B="$(stat -f%z "$KA.h5ad")"
if [ "$KA_SHA" = "$WANT_H5AD" ] && [ "$KA_B" = "$WANT_BYTES" ]; then MATCH=MATCH; else MATCH=MISMATCH; fi
log "KA h5ad sha256 $KA_SHA bytes $KA_B  (want $WANT_H5AD $WANT_BYTES)  $MATCH"
PB="$(grep -E 'checks pass' "$W/logs/postbuild.log" 2>/dev/null | tail -1)"; [ -n "$PB" ] || PB="$(grep -rhE 'checks pass' "$W" 2>/dev/null | tail -1)"
log "postbuild inside the driver: ${PB:-<no summary line found>}"

wait_gate
log "PREP start"
/usr/bin/time -l vcc prep "$KA.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$KA.vcc" > "$W/prep.log" 2>&1
PRC=$?; log "PREP rc=$PRC wall $(wall_of "$W/prep.log") s peak RSS $(rss_of "$W/prep.log") GiB"
if [ $PRC -eq 0 ] && [ -f "$KA.vcc" ]; then
  V_SHA="$(shasum -a 256 "$KA.vcc" | cut -d' ' -f1)"; V_B="$(stat -f%z "$KA.vcc")"
  [ "$V_SHA" = "$WANT_VCC" ] && VMATCH=MATCH || VMATCH=MISMATCH
  log "KA vcc sha256 $V_SHA bytes $V_B  (want $WANT_VCC)  $VMATCH"
else
  VMATCH=NOPREP; log "PREP FAILED; the .vcc row is not established"
fi
if [ "$MATCH" = "MATCH" ]; then rm -f "$KA.h5ad"; log "DELETED $KA.h5ad (byte duplicate of submissions/249-cd4svd5-f034.h5ad; meta kept)"; fi
if [ "$VMATCH" = "MATCH" ]; then rm -f "$KA.vcc"; log "DELETED $KA.vcc (byte duplicate of submissions/249-cd4svd5-f034.vcc)"; fi
log "KA ROW 1 RESULT: h5ad $MATCH, vcc $VMATCH"
[ "$MATCH" = "MATCH" ] && [ "$VMATCH" = "MATCH" ]
