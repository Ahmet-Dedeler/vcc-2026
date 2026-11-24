#!/usr/bin/env bash
# 495, levelling round 3. Round 2 left context C at rho 0.9627 (-1.7 %, outside the 0.5 % band chain1003 requires), while A
# and B converged. relevel_iter.py's secant (local exponent 2.373) gives C 1.1235 with A/B unchanged: A=1.0680,B=1.1215,
# C=1.1235. Runs after chain479b (process exit). The file is named chain495.sh on purpose: chain476's gate regex waits for
# `chain495\.sh`, so 476 cannot start beside this build. The round-2 build dir and package are moved aside (.r2) and the arm
# is rebuilt under the same tag, name and stem, so chain1003 (10-03), chain498 (495's meta) and read1003.py need no change.
# Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/495-lodo-cap6; LOG=$X/work/chain495r3.log; E=experiments/428-zero-gene-audit; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain495r3.pid
log "armed (pid $$): waiting for chain479b's process to exit (no other build), builder == HEAD"
P479B="$(cat experiments/479-pertarget-sharpen/work/chain479b.pid 2>/dev/null)"
until ! ps -p "${P479B:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain479b \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain479b?\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609291200 ] && { log "479b never ended by 12:00Z 09-29 -- stopping"; log "END"; exit 1; }; sleep 60; done
B="$(shasum -a 256 $BU | cut -d' ' -f1)"; EM="$(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)"
[ "${EM:0:16}" = edb5b8dcbee51968 ] || { log "emitter is ${EM:0:16}, not edb5b8dc (496p ran first?) -- stopping"; log "END"; exit 1; }
STEM=495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5; NM=the-water-bath-filed-for-divorce
N2="$REPO/$X/work/dump_h1x2_cd4x15_lodo05hi6.npz"; DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1)
[ "${DS2:0:16}" = d58750b83915d40d ] || { log "cap-6 dump sha ${DS2:0:16} unexpected -- stopping"; log "END"; exit 1; }
R="$REPO/$(ls submissions/446bj2-*.h5ad.meta.json | head -1)"; F3="A=1.0680,B=1.1215,C=1.1235"
log "released: builder ${B:0:16}; round-2 dir and package moved aside (.r2); building tag 495 at $F3"
[ -d $E/work/build_495 ] && mv $E/work/build_495 $E/work/build_495.r2
[ -f submissions/$STEM.vcc ] && mv submissions/$STEM.vcc submissions/$STEM.r2.vcc
[ -f submissions/$STEM.h5ad.meta.json ] && cp submissions/$STEM.h5ad.meta.json submissions/$STEM.r2.h5ad.meta.json
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN
  export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F3"
  bash $E/run_arm_generic.sh 495 1.75 "$R" "$NM" "$STEM" "$B" "$EM" bulk_ctx_scale,bulk_load,bulk_load_sha256 ) >> $X/work/build495r3.nohup 2>&1
log "495 (round 3): $(/usr/bin/grep -E '^\[495\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping' $E/work/build_495/run495.log | tr '\n' ' ' | cut -c1-500)"
H="$E/work/build_495/$NM.h5ad"
if [ -s "$H" ]; then
  uv run python experiments/417-bulk-sharpen/ntclevel417.py "$H" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1
  log "495 rho (round 3): $(uv run python scripts/def_day_helpers.py rho "${H%.h5ad}.ntclevel.json" 2>&1 | tail -1) (must be within 0.5 % of 0.9793)"
fi
rm -f "$H"
log "END"
