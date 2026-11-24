#!/usr/bin/env bash
# build_630arm.sh -- build the 630 arm = 446bj2z's recipe + VCC_DE_UPRN=b_k100_q20r3p10 on the POST-10-06 builder.
# UNARMED: the coordinator launches it (after 10-06 APPLY DAY step 8). It never submits.
#   bash experiments/630-rank-neutral-up/build_630arm.sh          # build
#   DRY=1 bash experiments/630-rank-neutral-up/build_630arm.sh    # print the plan, gates and env.json; build nothing
#
# Gates BEFORE the build (any miss -> NOT READY):
#   * builder (scripts/make_submission_pm_cd4.py) and emitter (src/vcc2026/rc_emit_k.py) == HEAD (the post-apply commit),
#     the builder carries 630's mark "630: rank-neutral UP", and src/vcc2026/rank_up.py is the pinned 7e6b3e52 module;
#   * no other build / vcc prep running (the same patterns build549.sh uses); disk >= 35 GiB; free+inactive >= 23 GiB
#     (run_arm_generic.sh then applies its own memory gate); the time cap CAP (UTC yyyymmddHHMM).
# The build: run_arm_generic.sh with 446bj2z's knobs (LODO bulk load, MDE coarse, zero repair -5, re-level scales, zs) plus
#   VCC_DE_UPRN and VCC_RC_AUDIT_DUMP. run_arm_generic.sh refuses a config diff vs 446bj2z that is not exactly the EXPECT
#   keys (de_uprn, plus de_sign_emitter_sha256 when the emitter pin moved), and runs the zero-gene audit, the cliff check,
#   vcc prep and the novelty check.
# Gates AFTER the build: the run log carries "[630arm] PASS: exactly", "CLIFF OK", the zero-gene audit's
#   "ZERO block total (target excluded): 0 in 0 blocks",
#   novelty rc 0, and a package sha; ka630_audit.py on the rc_audit dump exits 0 (FAILSAFE_TOTAL 0 + capacity floor in
#   A, B and C). The terminal line of the log is exactly one of:  READY 630arm <pkg sha256>  |  NOT READY <reason>
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/630-rank-neutral-up; W=$X/work/build630arm; mkdir -p "$W"; LOG=$W/build_630arm.log
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py; RUP=src/vcc2026/rank_up.py
E=experiments/428-zero-gene-audit; S=submissions
RUP_SHA=7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f
VARIANT=b_k100_q20r3p10
TAG=630arm; STEM="446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100"
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"
AUDIT="$W/$TAG.rc_audit.jsonl"
CAP="${CAP:-202610071800}"
DRY="${DRY:-0}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
notready(){ log "NOT READY $1"; exit "${2:-1}"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4|build527\.sh|fix526\.sh|build549\.sh/' | /usr/bin/grep -q . && return 0
        ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)' && return 0
        [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" != 0 ]; }
gates(){  # prints one reason per failing gate; returns 0 when all pass
  local bad=0
  eq_head $BU || { echo "builder $(sha $BU | cut -c1-16) != HEAD"; bad=1; }
  eq_head $EMP || { echo "emitter $(sha $EMP | cut -c1-16) != HEAD"; bad=1; }
  /usr/bin/grep -q "630: rank-neutral UP" $BU || { echo "builder lacks 630's mark (APPLY DAY not done)"; bad=1; }
  [ -f $RUP ] && [ "$(sha $RUP)" = "$RUP_SHA" ] || { echo "rank_up.py missing or not 7e6b3e52"; bad=1; }
  [ -f "$REFZ" ] || { echo "446bj2z meta missing"; bad=1; }
  [ "$(disk_gib)" -ge 35 ] || { echo "disk $(disk_gib) GiB < 35"; bad=1; }
  [ "$(fi_gib)" -ge 23 ] || { echo "free+inactive $(fi_gib) GiB < 23"; bad=1; }
  busy && { echo "another build or vcc prep is running"; bad=1; }
  return $bad
}
log "build_630arm (pid $$, DRY=$DRY, CAP=$CAP): 446bj2z + VCC_DE_UPRN=$VARIANT; stem $STEM"
[ -e "$S/$STEM.vcc" ] && notready "a package already exists at $S/$STEM.vcc (not overwriting)"
DLS="$(sha "$DL" 2>/dev/null)"; [ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || notready "LODO dump sha ${DLS:0:16} unexpected"
if [ "$DRY" = 1 ]; then
  R="$(gates)"; GRC=$?
  log "DRY: gates $( [ $GRC = 0 ] && echo PASS || echo "would wait/stop: $(echo "$R" | tr '\n' ';')" )"
else
  until R="$(gates)"; do
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && notready "gates not met by CAP $CAP: $(echo "$R" | tr '\n' ';')"
    sleep 60
  done
fi
NEWB="$(sha $BU)"; NEWE="$(sha $EMP)"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
EXPECT="de_uprn"; [ -n "$INC_EM" ] && [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
python3 - "$W/env.json" "$NEWB" "$NEWE" "$EXPECT" "$DL:$DLS" "$AUDIT" "$VARIANT" "$STEM" "$(git rev-parse --short HEAD)" <<'PY'
import json, sys
p, b, e, exp, dl, audit, var, stem, head = sys.argv[1:10]
json.dump({"arm": "630arm", "stem": stem, "git_head": head, "builder_sha256": b, "emitter_sha256": e,
           "rank_up_sha256": "7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f",
           "expect_config_diff_vs_446bj2z": exp.split(","),
           "env": {"VCC_BULK_LOAD": dl, "VCC_DE_MDE": "0.1,0.1,0.3,coarse", "VCC_DE_ZERO_REPAIR": "-5",
                   "VCC_BULK_CTX_SCALE": "A=1.0962,B=1.1461,C=1.1681", "VCC_DE_MAG": "dn_k250_zs", "VCC_DE_UPRN": var,
                   "VCC_RC_AUDIT_DUMP": audit, "q": "1.75"},
           "unset": ["VCC_SAMELINE", "VCC_DE_PRIOR_GATE", "VCC_UNION_SOURCES", "VCC_DE_MDE_SIGN", "VCC_RC_CALLSET_CPM",
                     "VCC_RC_CALLSET_BLOCK", "VCC_DE_DNMAG"]}, open(p, "w"), indent=1)
PY
log "env.json written ($W/env.json); builder ${NEWB:0:16} emitter ${NEWE:0:16}; expected diff vs 446bj2z: $EXPECT"
if [ "$DRY" = 1 ]; then log "DRY: would run run_arm_generic.sh $TAG 1.75 <446bj2z meta> uprn630 $STEM ${NEWB:0:16} ${NEWE:0:16} $EXPECT"; log "NOT READY dry run (nothing built)"; exit 0; fi
rm -f "$AUDIT"          # the builder appends: never let a stale dump stand
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" VCC_DE_MAG=dn_k250_zs
  unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK VCC_DE_DNMAG
  export VCC_DE_UPRN="$VARIANT" VCC_RC_AUDIT_DUMP="$REPO/$AUDIT"
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" uprn630 "$STEM" "$NEWB" "$NEWE" "$EXPECT" ) >> "$W/build_$TAG.nohup" 2>&1
RL=$E/work/build_$TAG/run$TAG.log
log "run log: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|zero-gene audit|package sha|novelty rc|DONE$|stopping|GATE" "$RL" 2>/dev/null | tr '\n' ' ' | cut -c1-800)"
[ "$(sha $BU)" = "$NEWB" ] && [ "$(sha $EMP)" = "$NEWE" ] || notready "builder or emitter moved during the build"
/usr/bin/grep -q "^\[$TAG\] PASS: exactly" "$RL" 2>/dev/null || notready "config diff is not exactly $EXPECT (see $RL)"
/usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" 2>/dev/null || notready "cliff check not OK (see $RL)"
ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || notready "zero-gene audit not 0: ${ZL:-missing}"
NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"
[ "$NOV" = "rc=0" ] || notready "novelty check ${NOV:-missing}"
PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
[ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] || notready "no package (see $RL)"
[ -s "$AUDIT" ] || notready "rc_audit dump $AUDIT missing or empty"
uv run python $X/ka630_audit.py "$AUDIT" --contexts A,B,C > "$W/ka630_audit.out" 2>&1
KRC=$?; tail -8 "$W/ka630_audit.out" | tee -a "$LOG"
[ $KRC -eq 0 ] || notready "ka630_audit rc $KRC (fail-safe trigger or capacity floor; see $W/ka630_audit.out)"
log "package $REPO/$S/$STEM.vcc sha256 $(sha "$S/$STEM.vcc")"
[ "$(sha "$S/$STEM.vcc")" = "$PK" ] || notready "package sha on disk differs from the run log's"
log "READY 630arm $PK"
