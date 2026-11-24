#!/usr/bin/env bash
# 535: 528's call-propensity block (emitter + builder patches) and the arm 446bj2zcp (PREREGISTRATION.md).
#  1. emitter: patch_emitter528.py on src/vcc2026/rc_emit_k.py (1658ae2c) -> must equal 528's tested k_emit528.py (0f56209d);
#     builder: patch_builder528.py on scripts/make_submission_pm_cd4.py (39c1dc66; pins the new emitter);
#     528's emitter tests (26/26, 23/23), test_patch528 on a builder copy, builder suite 56/0;
#  2. KA-1: 446bj2z's env, knob unset -> ec3a5ae6 byte for byte;
#     KA-2: 446bj2z's env + the INERT artifact (e62ee765, CPM-order key, steer kept) -> ec3a5ae6 byte for byte;
#     KA-mz: 446bj2zmz's env (the Arc hook), knob unset -> ea5aaa92 byte for byte;
#  3. the arm 446bj2zcp: 446bj2z's env + the A/B/C artifact (04fb24ec), diff vs 446bj2z exactly rc_callset_block +
#     de_sign_emitter_sha256.
# Starts >= 09-30 06:00Z (after the 09-30 read's possible builds; clear of 534's fits). Serial; any failure in 1-2 restores
# both files. commit_after_chain.sh commits them (commit535.spec). Never submits. Named build535.sh (no chain5NN trap).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/535-callprop-arm; W=$X/work; LOG=$W/build535.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
P528=experiments/528-callprop-block; E=experiments/428-zero-gene-audit; S=submissions; KA=experiments/527-arc-hook-arm/run_ka_527.sh
ART=$P528/work/callblock528_abc.npz; ARTS=04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd
INERT=$P528/work/callblock528_inert_abc.npz; INERTS=e62ee76539f6932728f9c58bd37d72a8e4c8cd74d4657ae7251cecfe62d23b69
EM_NEW=0f56209d57e8385a
mkdir -p $W
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4|build527\.sh|fix526\.sh/' | /usr/bin/grep -q . && return 0
        [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" != 0 ]; }
echo $$ > $W/build535.pid
log "armed (pid $$): waiting for >= 09-30 06:00Z, no other build, builder and emitter == HEAD, free+inactive >= 23 GiB"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609300600 ] && ! busy && eq_head $BU && eq_head $EMP && [ "$(fi_gib)" -ge 23 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610031200 ] && { log "prerequisites not met by 12:00Z 10-03 -- stopping"; log "END"; exit 1; }
  sleep 60; done
for f in $ART:$ARTS $INERT:$INERTS; do [ "$(shasum -a 256 "${f%%:*}" | cut -d' ' -f1)" = "${f##*:}" ] || { log "artifact ${f%%:*} sha differs -- stopping"; log "END"; exit 1; }; done
OLDB="$(shasum -a 256 $BU | cut -d' ' -f1)"; OLDE="$(shasum -a 256 $EMP | cut -d' ' -f1)"; cp $BU $W/builder-pre535.py; cp $EMP $W/emitter-pre535.py
back(){ cp $W/builder-pre535.py $BU; cp $W/emitter-pre535.py $EMP; log "RESTORED builder ${OLDB:0:16} emitter ${OLDE:0:16} ($1)"; log "END"; exit 1; }
log "released: builder ${OLDB:0:16}, emitter ${OLDE:0:16}"
# --- 1. patches + tests. 528's tests run BEFORE either file is swapped in: test_patch528's P5 re-patches the production
# emitter and must see the unpatched one (09-30 06:00Z: the first run swapped first, P5 found no anchor, the chain restored).
uv run python $P528/patch_emitter528.py --src $EMP --out $W/rc_emit_k.new.py >> $LOG 2>&1 || back "emitter patch failed"
[ "$(shasum -a 256 $W/rc_emit_k.new.py | cut -c1-16)" = "$EM_NEW" ] || back "patched emitter is not 528's tested k_emit528.py ($EM_NEW)"
cmp -s $W/rc_emit_k.new.py $P528/k_emit528.py || back "patched emitter differs from k_emit528.py"
cp $BU $W/builder.new.py
uv run python $P528/patch_builder528.py $W/builder.new.py >> $LOG 2>&1 || back "builder patch failed"
[ "$(shasum -a 256 $W/builder.new.py | cut -c1-16)" = 4daae237c5fe4fd2 ] || back "patched builder is not 4daae237"
PYTHONDONTWRITEBYTECODE=1 uv run python $P528/test_emitter528v2.py > $W/test_emitter528v2.log 2>&1; log "test_emitter528v2: $(tail -1 $W/test_emitter528v2.log)"
/usr/bin/grep -q "^26/26 PASS" $W/test_emitter528v2.log || back "test_emitter528v2 failed"
PYTHONDONTWRITEBYTECODE=1 uv run python $P528/test_emitter528.py > $W/test_emitter528.log 2>&1; log "test_emitter528: $(tail -1 $W/test_emitter528.log)"
/usr/bin/grep -q "^23/23 PASS" $W/test_emitter528.log || back "test_emitter528 failed"
T=$W/patchtest; mkdir -p $T; cp $W/builder-pre535.py $T/b.py
PYTHONDONTWRITEBYTECODE=1 uv run python $P528/test_patch528.py $T/b.py > $W/test_patch528.log 2>&1; log "test_patch528: $(tail -1 $W/test_patch528.log)"
/usr/bin/grep -q "^22/22 PASS" $W/test_patch528.log || back "test_patch528 failed"
eq_head $EMP && eq_head $BU || back "production files moved during the tests"
cp $W/rc_emit_k.new.py $EMP.tmp535 && mv $EMP.tmp535 $EMP
cp $W/builder.new.py $BU.tmp535 && mv $BU.tmp535 $BU
NEWB="$(shasum -a 256 $BU | cut -d' ' -f1)"; NEWE="$(shasum -a 256 $EMP | cut -d' ' -f1)"
log "installed: builder ${OLDB:0:16} -> ${NEWB:0:16}, emitter ${OLDE:0:16} -> ${NEWE:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $W/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $W/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $W/builder-tests-after.log || back "builder suite not 56/0"
# --- 2. known answers
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
REFM="$(ls $REPO/$S/446bj2zmz-*.h5ad.meta.json | head -1)"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || back "LODO dump sha ${DLS:0:16} unexpected"
MZ=experiments/517-flex-programme/work/arch1meanz_ne_hook_abc.npz; MZS=165c6b78a16de2990ed87ca0e3436ab8c4023f7dddcdb30cf00b8f5eebda9f6b
base(){ export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" VCC_DE_MAG=dn_k250_zs
        unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK; }
ka(){ local tag="$1" ref="$2" want="$3"; shift 3
  ( base; for kv in "$@"; do export "$kv"; done; KA_Q=1.75 bash $KA "$tag" "$ref" "$NEWB" ) >> $W/ka_$tag.nohup 2>&1
  local l=experiments/527-arc-hook-arm/work/build_$tag/run$tag.log
  /usr/bin/grep -q "KA PASS (identical)" "$l" 2>/dev/null || back "$tag did not reproduce ${want}: $(/usr/bin/grep -E 'KA FAIL|driver rc|stopping' "$l" | tail -2 | tr '\n' ' ')"
  [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEWB" ] || back "builder moved during $tag"
  log "$tag PASS: rebuilds ${want} byte for byte on builder ${NEWB:0:16} / emitter ${NEWE:0:16}"; }
ka KA535a "$REFZ" ec3a5ae6
ka KA535b "$REFZ" "ec3a5ae6 (inert artifact)" "VCC_RC_CALLSET_BLOCK=prop:$INERT:$INERTS"
ka KA535c "$REFM" "ea5aaa92 (446bj2zmz)" "VCC_DE_MDE_SIGN=file:$MZ:$MZS:0.05:0.05:le1"
# --- 3. the arm
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
EXPECT="rc_callset_block"; [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
TAG=446bj2zcp; STEM="$TAG-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cp"
log "arm $TAG: expected diff vs 446bj2z: $EXPECT"
( base; export VCC_RC_CALLSET_BLOCK="prop:$ART:$ARTS"
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" callprop535 "$STEM" "$NEWB" "$NEWE" "$EXPECT" ) >> $W/build_$TAG.nohup 2>&1
log "$TAG: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping" $E/work/build_$TAG/run$TAG.log | tr '\n' ' ' | cut -c1-700)"
rm -f "$E/work/build_$TAG/callprop535.h5ad"
log "END builder ${NEWB:0:16} emitter ${NEWE:0:16} (committed by commit_after_chain.sh)"
