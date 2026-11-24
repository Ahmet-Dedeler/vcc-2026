#!/usr/bin/env bash
# 535 rebuild of the arm 446bj2zcp (2026-09-30 10:55Z). build535's arm passed every gate at 08:36Z (config diff exactly
# rc_callset_block + de_sign_emitter_sha256, zero-gene audit 0, CLIFF OK; artifact sha 0b3733d5...), but run_arm_generic's
# prep gate never opened (disk 17.7 GiB < 20 GiB) and build535's last line then removed the only copy of the prediction
# (rm -f callprop535.h5ad, which was meant for after the package existed). 47 GiB were freed (417's 09-23 build trees,
# git clean -xdf).
# This rebuilds the arm with build535's env on the same committed builder / emitter (4daae237 / 0f56209d, == HEAD) and
# requires the new artifact to equal the first build's sha byte for byte; on any mismatch it writes
# 464/work/NOGO-1001-s2 so the slot-2 waiter stops before its submit. The first build dir is moved aside
# (build_446bj2zcp.nopkg-0930) so this run's log holds only this run's gate lines. The prediction h5ad stays until the arm
# is scored. Never submits; arm1001cp.sh (re-armed on work/rebuild535.pid) arms 10-01 slot 2 if the package is ready.
# Named rebuild535.sh: matches no other chain's builder pattern (its child run_arm_generic.sh does, as it should).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/535-callprop-arm; W=$X/work; LOG=$W/rebuild535.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
P528=experiments/528-callprop-block; E=experiments/428-zero-gene-audit; S=submissions
ART=$P528/work/callblock528_abc.npz; ARTS=04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd
FIRST=0b3733d5e6e0b2f0c8c7138cbc767f4c42a29a2387a56516f3ef01f907bb2b2c   # build535's artifact (first build dir's run log)
TAG=446bj2zcp; STEM="$TAG-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cp"
BD=$E/work/build_$TAG; RL=$BD/run$TAG.log; NOGO=experiments/464-slot-1001/work/NOGO-1001-s2
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fin(){ log "$1"; log "END"; exit "${2:-1}"; }
nogo(){ printf '%s rebuild535: %s\n' "$(date -u +%FT%TZ)" "$1" > "$NOGO"; fin "NOGO written ($NOGO): $1" 1; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4|build527\.sh|fix526\.sh/' | /usr/bin/grep -q . && return 0
        [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" != 0 ]; }
mkdir -p $W; echo $$ > $W/rebuild535.pid
log "armed (pid $$): waiting for no other build, builder and emitter == HEAD, free+inactive >= 23 GiB, disk >= 35 GiB (cap 10-01 00:00Z)"
until ! busy && eq_head $BU && eq_head $EMP && [ "$(fi_gib)" -ge 23 ] && [ "$(disk_gib)" -ge 35 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610010000 ] && fin "prerequisites not met by 10-01 00:00Z -- stopping (no slot 2)"
  sleep 60; done
NEWB="$(shasum -a 256 $BU | cut -d' ' -f1)"; NEWE="$(shasum -a 256 $EMP | cut -d' ' -f1)"
[ "${NEWB:0:16}" = 4daae237c5fe4fd2 ] && [ "${NEWE:0:16}" = 0f56209d57e8385a ] || fin "builder ${NEWB:0:16} / emitter ${NEWE:0:16} are not 535's committed 4daae237 / 0f56209d -- stopping"
[ "$(shasum -a 256 "$ART" | cut -d' ' -f1)" = "$ARTS" ] || fin "artifact $ART sha differs -- stopping"
[ -e "$S/$STEM.vcc" ] && fin "a package already exists at $S/$STEM.vcc -- not overwriting it; stopping"
if [ -d "$BD" ]; then
  [ -e "$BD/callprop535.h5ad" ] && fin "$BD still holds a prediction -- this script only rebuilds a removed one; stopping"
  mv "$BD" "$BD.nopkg-0930" || fin "could not move $BD aside"; log "moved the first build dir aside: $BD.nopkg-0930"
fi
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "${DLS:0:16}" = 00ccfcba1b67fa24 ] || fin "LODO dump sha ${DLS:0:16} unexpected -- stopping"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
EXPECT="rc_callset_block"; [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "released: builder ${NEWB:0:16} emitter ${NEWE:0:16}, disk $(disk_gib) GiB; arm $TAG, expected diff vs 446bj2z: $EXPECT"
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" VCC_DE_MAG=dn_k250_zs
  unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK
  export VCC_RC_CALLSET_BLOCK="prop:$ART:$ARTS"
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" callprop535 "$STEM" "$NEWB" "$NEWE" "$EXPECT" ) >> $W/rebuild_$TAG.nohup 2>&1
log "$TAG: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping|GATE" $RL | tr '\n' ' ' | cut -c1-700)"
GOT="$(/usr/bin/grep -E ' artifact: [0-9]+ ' $RL | tail -1 | awk '{print $NF}')"
[ "$GOT" = "$FIRST" ] || nogo "the rebuilt artifact ${GOT:0:16} is not the first build's ${FIRST:0:16}"
log "REPRODUCED: the rebuilt artifact equals the first build's ${FIRST:0:16} byte for byte"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEWB" ] && [ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$NEWE" ] || nogo "builder or emitter moved during the rebuild"
PK="$(/usr/bin/grep 'package sha256:' $RL | tail -1 | awk '{print $NF}')"
[ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] || fin "no package after the rebuild (see $RL) -- no slot 2"
log "package $S/$STEM.vcc sha ${PK:0:16}"
fin "DONE (builder ${NEWB:0:16} emitter ${NEWE:0:16})" 0
