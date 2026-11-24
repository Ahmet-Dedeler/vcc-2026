#!/usr/bin/env bash
# 527: 517's programme-sourced uncalled-gene hook on the board (PREREGISTRATION.md, written before any 09-30 entry).
#  1. patch_builder517.py on the production builder (VCC_DE_MDE_SIGN gains `programme:` and `file:` sources; `clean4:`
#     parses exactly as 496's); builder suite 56/0; 517's test_patch517.py 26/26 on copies of the patched and unpatched
#     builders (ROOT pinned, as 517 ran it);
#  2. KA-1: 446bj2z's env with VCC_DE_MDE_SIGN unset rebuilds 446bj2z's output ec3a5ae6 byte for byte;
#  3. KA-496: 496a's env (446bj2's + clean4:0.05) rebuilds 496a's output e21ff752 byte for byte;
#  4. the arm 446bj2zmz = 446bj2z + VCC_DE_MDE_SIGN=file:<517's meanz_ne hook artifact>:<sha>:0.05:0.05:le1 (DE-only, no
#     re-level), config diff vs 446bj2z exactly de_mde_sign + de_sign_emitter_sha256.
# Serial: no other builder, no build chain, builder and emitter == HEAD. Any failure in 1-3 restores the pre-527 builder.
# commit_after_chain.sh commits the builder (spec commit527.spec). Never submits. Named build527.sh (not chain5NN.sh): the
# 524 rehearsal's gate waits on `chain5..` names, and this script must not wait on itself (AGENTS §6.2).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/527-arc-hook-arm; W=$X/work; LOG=$W/build527.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
E=experiments/428-zero-gene-audit; S=submissions; P517=experiments/517-flex-programme
MZ=experiments/517-flex-programme/work/arch1meanz_ne_hook_abc.npz; MZS=165c6b78a16de2990ed87ca0e3436ab8c4023f7dddcdb30cf00b8f5eebda9f6b
mkdir -p $W
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
echo $$ > $W/build527.pid
log "armed (pid $$): waiting for no other build, builder and emitter == HEAD"
until ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" = 0 ] \
      && eq_head $BU && eq_head $EMP; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610031200 ] && { log "prerequisites not met by 12:00Z 10-03 -- stopping"; log "END"; exit 1; }; sleep 60; done
[ "$(shasum -a 256 $MZ | cut -d' ' -f1)" = "$MZS" ] || { log "hook artifact sha differs from $MZS -- stopping"; log "END"; exit 1; }
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $W/builder-pre527.py
back(){ cp $W/builder-pre527.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
log "released: builder ${OLD:0:16}, emitter $(shasum -a 256 $EMP | cut -c1-16)"
# --- 1. patch + tests
uv run python $P517/patch_builder517.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $W/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $W/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $W/builder-tests-after.log || back "builder suite not 56/0"
T=$W/patchtest; mkdir -p $T; cp $BU $T/b517.py; cp $W/builder-pre527.py $T/head.py
for f in $T/b517.py $T/head.py; do
  uv run python - "$f" <<'PY' || back "ROOT pin failed on $f"
import sys
p = sys.argv[1]; t = open(p).read(); old = 'ROOT = Path(__file__).resolve().parents[1]\n'
assert t.count(old) == 1
open(p, 'w').write(t.replace(old, 'ROOT = Path("~/Code/virtual-cell-challenge-2026")  # 527 test copy only\n'))
PY
done
uv run python $P517/test_patch517.py $T/b517.py $T/head.py > $W/test_patch517.log 2>&1
log "test_patch517: $(tail -1 $W/test_patch517.log)"
/usr/bin/grep -q "ALL 26 PASS" $W/test_patch517.log || back "test_patch517 not 26/26"
rm -f $T/b517.py $T/head.py
# --- 2. KA-1 (knob unset == 446bj2z)
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"][:16])' "$REFZ")" = ec3a5ae6167ab2d6 ] || back "446bj2z's meta is not ec3a5ae6"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || back "LODO dump sha ${DLS:0:16} unexpected"
F0="A=1.0962,B=1.1461,C=1.1681"
base(){ export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0"
        unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_DE_MAG; }
( base; export VCC_DE_MAG=dn_k250_zs; KA_Q=1.75 bash $X/run_ka_527.sh KA1 "$REFZ" "$NEW" ) >> $W/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS (identical)" $X/work/build_KA1/runKA1.log 2>/dev/null || back "KA-1 did not reproduce 446bj2z: $(/usr/bin/grep -E 'KA FAIL|driver rc|stopping' $X/work/build_KA1/runKA1.log | tail -2 | tr '\n' ' ')"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEW" ] || back "builder moved during KA-1"
log "KA-1 PASS: 446bj2z's env with VCC_DE_MDE_SIGN unset rebuilds ec3a5ae6 byte for byte on builder ${NEW:0:16}"
# --- 3. KA-496 (clean4:0.05 on 446bj2's env == 496a)
REF4="$REPO/$S/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05.h5ad.meta.json"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"][:16])' "$REF4")" = e21ff752ce0349d7 ] || back "496a's meta is not e21ff752"
( base; export VCC_DE_MDE_SIGN=clean4:0.05; KA_Q=1.75 bash $X/run_ka_527.sh KA496 "$REF4" "$NEW" ) >> $W/ka496.nohup 2>&1
/usr/bin/grep -q "KA PASS (identical)" $X/work/build_KA496/runKA496.log 2>/dev/null || back "KA-496 did not reproduce 496a: $(/usr/bin/grep -E 'KA FAIL|driver rc|stopping' $X/work/build_KA496/runKA496.log | tail -2 | tr '\n' ' ')"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEW" ] || back "builder moved during KA-496"
log "KA-496 PASS: 496a's env (clean4:0.05) rebuilds e21ff752 byte for byte on builder ${NEW:0:16}"
# --- 4. the arm (the builder stays patched from here on: the committer commits it)
EM="$(shasum -a 256 $EMP | cut -d' ' -f1)"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
EXPECT="de_mde_sign"; [ "$EM" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "arm 446bj2zmz: expected diff vs 446bj2z: $EXPECT"
TAG=446bj2zmz; STEM="$TAG-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-hkmz05le1"
( base; export VCC_DE_MAG=dn_k250_zs VCC_DE_MDE_SIGN="file:$MZ:$MZS:0.05:0.05:le1"
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" hookmz527 "$STEM" "$NEW" "$EM" "$EXPECT" ) >> $W/build_$TAG.nohup 2>&1
log "$TAG: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping" $E/work/build_$TAG/run$TAG.log | tr '\n' ' ' | cut -c1-700)"
log "$TAG de_mde_sign: $(/usr/bin/grep -h '\[de_mde_sign\]' $E/work/build_$TAG/logs/build.log 2>/dev/null | tr '\n' ' ' | cut -c1-400)"
rm -f "$E/work/build_$TAG/hookmz527.h5ad"
log "END builder ${NEW:0:16} (committed by commit_after_chain.sh)"
