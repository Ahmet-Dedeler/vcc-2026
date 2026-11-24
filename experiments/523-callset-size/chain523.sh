#!/usr/bin/env bash
# 523: the DE call-set size knob, then the preregistered dose pair (PREREGISTRATION.md, written before any build).
#  1. patch_builder523.py (VCC_RC_CALLSET_CPM; unset = 2200); builder suite 56/0;
#  2. KA0: 446bj2z's env with the knob UNSET rebuilds 446bj2z's output ec3a5ae6 byte for byte (run_ka_523.sh);
#  3. the arms 446bj2zk16 (cpm 1200, K 1,600) and 446bj2zk20 (cpm 1600, K 2,000): 446bj2z's env + the knob, config diff
#     vs 446bj2z exactly rc_callset_cpm (+ de_sign_emitter_sha256 if 496p's pin differs from 446bj2z's recorded one).
# Serial after the 09-28 queue's last build (chain498r3, process exit; no other build; builder and emitter == HEAD).
# Any failure in 1-2 restores the pre-523 builder. commit_after_chain.sh commits the builder (spec commit523.spec).
# Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/523-callset-size; W=$X/work; LOG=$W/chain523.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
E=experiments/428-zero-gene-audit; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
echo $$ > $W/chain523.pid
log "armed (pid $$): waiting for chain498r3's process to exit, no other build, builder and emitter == HEAD"
P3="$(cat experiments/498-combined-bulk/work/chain498r3.pid 2>/dev/null)"
until ! ps -p "${P3:-1}" -o args= 2>/dev/null | /usr/bin/grep -q relevel_sp2 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|chain497e\.sh|chain51[01]\.sh|chain502\.sh|run_ka_4|run_ka_5|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && eq_head $BU && eq_head $EMP; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610030000 ] && { log "prerequisites not met by 00:00Z 10-03 -- stopping"; log "END"; exit 1; }; sleep 60; done
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $W/builder-pre523.py
back(){ cp $W/builder-pre523.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
log "released: builder ${OLD:0:16}, emitter $(shasum -a 256 $EMP | cut -c1-16)"
uv run python $X/patch_builder523.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $W/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $W/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $W/builder-tests-after.log || back "builder suite not 56/0"
REF="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"][:16])' "$REF")" = ec3a5ae6167ab2d6 ] || back "446bj2z's meta is not ec3a5ae6"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || back "LODO dump sha ${DLS:0:16} unexpected"
envz(){ export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_DE_MAG=dn_k250_zs \
               VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681"; unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN; }
# --- 2. KA0 (knob unset)
( envz; unset VCC_RC_CALLSET_CPM; KA_Q=1.75 bash $X/run_ka_523.sh KA0 "$REF" "$NEW" ) >> $W/ka0.nohup 2>&1
KAL="$X/work/build_KA0/runKA0.log"
/usr/bin/grep -q "KA PASS (identical)" "$KAL" 2>/dev/null || back "KA0 did not reproduce 446bj2z: $(/usr/bin/grep -E 'KA FAIL|driver rc|stopping' "$KAL" | tail -2 | tr '\n' ' ')"
eq_head $BU && back "builder equals HEAD after the patch (the patch did not land)"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEW" ] || back "builder moved during KA0"
log "KA0 PASS: 446bj2z's env with VCC_RC_CALLSET_CPM unset rebuilds ec3a5ae6 byte for byte on builder ${NEW:0:16}"
# --- 3. the arms
EM="$(shasum -a 256 $EMP | cut -d' ' -f1)"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REF")"
SRC_EM="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")' $BU)"
EXPECT="rc_callset_cpm"; [ -n "$SRC_EM" ] && [ "$SRC_EM" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "arms: expected diff vs 446bj2z: $EXPECT"
for arm in "k16:1200:the-stir-bar-counted-sheep" "k20:1600:the-vortexer-hummed-a-tune"; do
  a="${arm%%:*}"; rest="${arm#*:}"; cpm="${rest%%:*}"; nm="${rest#*:}"; tag="446bj2z$a"
  stem="$tag-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm$cpm"
  ( envz; export VCC_RC_CALLSET_CPM=$cpm
    bash $E/run_arm_generic.sh "$tag" 1.75 "$REF" "$nm" "$stem" "$NEW" "$EM" "$EXPECT" ) >> $W/build_$tag.nohup 2>&1
  log "$tag: $(/usr/bin/grep -E "^\[$tag\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping" $E/work/build_$tag/run$tag.log | tr '\n' ' ' | cut -c1-600)"
  rm -f "$E/work/build_$tag/$nm.h5ad"
done
log "END builder ${NEW:0:16} (committed by commit_after_chain.sh)"
