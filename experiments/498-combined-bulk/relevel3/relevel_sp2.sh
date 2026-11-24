#!/usr/bin/env bash
# 498sp2, levelling round 3 (2026-09-28). Round 2 left context C at rho 0.9714 (-0.81 %, outside the 0.5 % band chain1004
# requires) while A (0.9802) and B (0.9779) converged. relevel_iter.py's secant on the two measured rounds (C: 1.1235 ->
# 1.0005, 1.1116 -> 0.9714; local exponent 2.780) gives C 1.1149 with A/B unchanged: A=1.0758,B=1.1163,C=1.1149.
# Runs after chain502 (the end of the 09-28 queue; process exit, no other build, builder and emitter == HEAD). The file
# must NOT be named chain49[5-8].sh: every earlier chain's "no other build" gate matches that pattern, and a WAITING
# process with that name deadlocked the queue (2026-09-28 15:50Z, renamed from chain498.sh). Its build is still visible
# to the gates through run_arm_generic.sh; chain523 waits for this script's own pid (chain498r3.pid). The round-2 build dir and package are
# moved aside (.r2) and the arm is rebuilt under the same tag, name and stem, so chain1004 (10-04), chain1005 and
# rule1004.py need no change. If 496p's emitter pin differs from 495's reference meta, de_sign_emitter_sha256 joins the
# expected config diff (499 F2: the pin moves that key on every build; 496p's KA-1 shows no output change). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/498-combined-bulk; LOG=$X/work/chain498r3.log; E=experiments/428-zero-gene-audit; BU=scripts/make_submission_pm_cd4.py
EMP=src/vcc2026/rc_emit_k.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
echo $$ > $X/work/chain498r3.pid
log "armed (pid $$): waiting for chain502's process to exit (no other build), builder and emitter == HEAD"
P502="$(cat experiments/502-october-de-riders/work/chain502.pid 2>/dev/null)"
until ! ps -p "${P502:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain502 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-7])\.sh|chain496p\.sh|chain497e\.sh|chain51[01]\.sh|chain502\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && eq_head $BU && eq_head $EMP; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610021200 ] && { log "prerequisites not met by 12:00Z 10-02 -- stopping"; log "END"; exit 1; }; sleep 60; done
B="$(shasum -a 256 $BU | cut -d' ' -f1)"; EM="$(shasum -a 256 $EMP | cut -d' ' -f1)"
STEM=498sp2-cap15-partial50-clean4q80-lodo6-ballast-sharpensp2-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5; NM=the-cuvette-learned-to-juggle
N2="$REPO/experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz"; DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1)
[ "${DS2:0:16}" = d58750b83915d40d ] || { log "cap-6 dump sha ${DS2:0:16} unexpected -- stopping"; log "END"; exit 1; }
R="$REPO/submissions/495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json"
[ -s "$R" ] || { log "495's meta missing -- stopping"; log "END"; exit 1; }
F3="A=1.0758,B=1.1163,C=1.1149"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$R")"
SRC_EM="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")' $BU)"
EXPECT="bulk_ctx_scale,bulk_sharpen_pt,bulk_sharpen_q"
[ -n "$SRC_EM" ] && [ "$SRC_EM" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "released: builder ${B:0:16}, emitter ${EM:0:16}; expected diff vs 495: $EXPECT; round-2 dir and package moved aside (.r2); building tag 498sp2 at $F3"
[ -d $E/work/build_498sp2 ] && mv $E/work/build_498sp2 $E/work/build_498sp2.r2
[ -f submissions/$STEM.vcc ] && mv submissions/$STEM.vcc submissions/$STEM.r2.vcc
[ -f submissions/$STEM.h5ad.meta.json ] && cp submissions/$STEM.h5ad.meta.json submissions/$STEM.r2.h5ad.meta.json
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN
  export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F3"
  bash $E/run_arm_generic.sh 498sp2 "stercile:2.5,2.0,1.75" "$R" "$NM" "$STEM" "$B" "$EM" "$EXPECT" ) >> $X/work/build498sp2r3.nohup 2>&1
log "498sp2 (round 3): $(/usr/bin/grep -E '^\[498sp2\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping' $E/work/build_498sp2/run498sp2.log | tr '\n' ' ' | cut -c1-500)"
H="$E/work/build_498sp2/$NM.h5ad"
if [ -s "$H" ]; then
  uv run python experiments/417-bulk-sharpen/ntclevel417.py "$H" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1
  log "498sp2 rho (round 3): $(uv run python scripts/def_day_helpers.py rho "${H%.h5ad}.ntclevel.json" 2>&1 | tail -1) (must be within 0.5 % of 0.9793)"
fi
rm -f "$H"
log "END"
