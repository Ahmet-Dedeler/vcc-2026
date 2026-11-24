#!/usr/bin/env bash
# 549: build the arm 446bj2zcph = 446bj2z + 528's call-propensity block (the same A/B/C artifact as 446bj2zcp) + 496's CLEAN4
# uncalled-gene hook (VCC_DE_MDE_SIGN=clean4:0.05, as 496a). Worker 537 (v3, §3.2) priced block + hook on H1 at +0.0059 …
# +0.0067 vs 446bj2z (block alone +0.0045 … +0.0053): super-additive, because the hook signs the block's swapped-in genes.
# Launched by hand only if 446bj2zcp reads >= PARTIAL on 10-01 (or NULL/NEG with d jac >= +0.002 and d fid < 0: 537 §0 C).
# Config diff vs 446bj2z must be exactly rc_callset_block + de_mde_sign (+ de_sign_emitter_sha256 when the pin moved).
# Pauses 537's driver (SIGSTOP between jobs) for the build and always resumes it (SIGCONT). Never submits.
# Named build549.sh: matches no other chain's builder pattern (its child run_arm_generic.sh does, as it should).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/549-cph-arm; W=$X/work; LOG=$W/build549.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
E=experiments/428-zero-gene-audit; S=submissions; P528=experiments/528-callprop-block
ART=$P528/work/callblock528_abc.npz; ARTS=04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd
TAG=446bj2zcph; STEM="$TAG-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph"
D537=experiments/537-arm-pipeline/work; DP=""
mkdir -p $W; echo $$ > $W/build549.pid
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
resume537(){ if [ -n "$DP" ] && ps -p "$DP" >/dev/null 2>&1; then kill -CONT "$DP" 2>/dev/null; log "537 driver $DP resumed (SIGCONT)";
             printf '%s resumed by the coordinator after build549\n' "$(date -u +%FT%TZ)" >> $D537/PAUSED-BY-COORDINATOR.txt; fi; }
trap resume537 EXIT
fin(){ log "$1"; log "END"; exit "${2:-1}"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4|build527\.sh|fix526\.sh/' | /usr/bin/grep -q . && return 0
        [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" != 0 ]; }
log "armed (pid $$): builder/emitter == HEAD (4daae237 / 0f56209d), no other build, disk >= 35 GiB (cap 10-06 12:00Z)"
[ -e "$S/$STEM.vcc" ] && fin "a package already exists at $S/$STEM.vcc -- not overwriting it"
until ! busy && eq_head $BU && eq_head $EMP && [ "$(disk_gib)" -ge 35 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610061200 ] && fin "prerequisites not met by 10-06 12:00Z -- stopping"
  sleep 60; done
NEWB="$(shasum -a 256 $BU | cut -d' ' -f1)"; NEWE="$(shasum -a 256 $EMP | cut -d' ' -f1)"
[ "${NEWB:0:16}" = 4daae237c5fe4fd2 ] && [ "${NEWE:0:16}" = 0f56209d57e8385a ] || fin "builder ${NEWB:0:16} / emitter ${NEWE:0:16} are not 4daae237 / 0f56209d -- stopping"
[ "$(shasum -a 256 "$ART" | cut -d' ' -f1)" = "$ARTS" ] || fin "artifact $ART sha differs -- stopping"
# pause 537's driver between jobs: stop the bash driver (its running job finishes), then wait for that job to exit
DP="$(cat $D537/drive537.pid 2>/dev/null)"
if [ -n "$DP" ] && ps -p "$DP" -o args= 2>/dev/null | /usr/bin/grep -q drive537; then
  kill -STOP "$DP" && log "537 driver $DP paused (SIGSTOP); waiting for its running job to exit"
  printf '%s paused by the coordinator for build549 (446bj2zcph); resumed automatically when the build ends\n' "$(date -u +%FT%TZ)" >> $D537/PAUSED-BY-COORDINATOR.txt
  # the job is the driver's direct child /usr/bin/time (-> uv -> python); its per-job watchdog subshell (bash) stays while the
  # driver is stopped (it waits for .run_done), so only the job's own processes are waited for
  # (a finished job stays a zombie until the stopped driver reaps it: state Z is excluded)
  until [ -z "$(pgrep -P "$DP" | xargs -n1 ps -o stat=,comm= -p 2>/dev/null | awk '$1 !~ /^Z/ {print $2}' | /usr/bin/grep -E '(^|/)(time|uv|python[0-9.]*)$')" ]; do sleep 20; done
  log "537's running job has exited"
else DP=""; log "537's driver is not running"; fi
until [ "$(fi_gib)" -ge 23 ]; do [ "$(date -u +%Y%m%d%H%M)" -ge 202610061200 ] && fin "memory never freed by 10-06 12:00Z"; sleep 30; done
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || fin "LODO dump sha ${DLS:0:16} unexpected -- stopping"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
EXPECT="rc_callset_block,de_mde_sign"; [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "released: builder ${NEWB:0:16} emitter ${NEWE:0:16}, disk $(disk_gib) GiB, free+inactive $(fi_gib) GiB; expected diff vs 446bj2z: $EXPECT"
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" VCC_DE_MAG=dn_k250_zs
  unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK
  export VCC_RC_CALLSET_BLOCK="prop:$ART:$ARTS" VCC_DE_MDE_SIGN=clean4:0.05
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" cph549 "$STEM" "$NEWB" "$NEWE" "$EXPECT" ) >> $W/build_$TAG.nohup 2>&1
RL=$E/work/build_$TAG/run$TAG.log
log "$TAG: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping|GATE" $RL | tr '\n' ' ' | cut -c1-700)"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEWB" ] && [ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$NEWE" ] || fin "builder or emitter moved during the build"
PK="$(/usr/bin/grep 'package sha256:' $RL | tail -1 | awk '{print $NF}')"
[ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] || fin "no package (see $RL)"
fin "DONE: package $S/$STEM.vcc sha ${PK:0:16} (the h5ad stays until the arm is scored)" 0
