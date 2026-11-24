#!/usr/bin/env bash
# 701 (science session [bf003b], 2026-10-14): build the K-UP arm of the FROZEN preregistration
# (experiments/701-kup-callset/PREREGISTRATION.md, 9b641a0a + its H1 amendment), modelled on 523's chain step 3 and
# 671's step 1.
#   nohup bash experiments/701-kup-callset/build701.sh > experiments/701-kup-callset/work/build701.nohup 2>&1 &
# ARM  446bj2zku30 = 446bj2z's exact env + VCC_RC_CALLSET_CPM=3000, built by run_arm_generic.sh against 446bj2z's meta
#      (ec3a5ae6), EXPECT rc_callset_cpm (+ de_sign_emitter_sha256 if the emitter pin differs from 446bj2z's record).
# Gates: builder AND emitter == HEAD (671's builder committed after its KAu), the LODO dump sha, disk >= 35 GiB,
#        free+inactive >= 23 GiB, no other build or vcc prep.
# After: "[ku30] PASS: exactly", CLIFF OK, zero-gene audit 0, novelty rc 0, package sha; pds/mse must be untouched (checked
#        on the board read, as for K16/K20). Terminal lines: "READY ku30 <pkg sha>" | "NOT READY <reason>"; then "END".
# Never submits; the coordinator [400b6e] arms the slot.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/701-kup-callset; W=$X/work; mkdir -p "$W"; LOG=$W/build701.log
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
E=experiments/428-zero-gene-audit; S=submissions
TAG=ku30
STEM="446bj2zku30-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm3000"
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
AUDIT="$REPO/$W/$TAG.rc_audit.jsonl"
CAP=202610161200   # 10-15 04:5xZ: moved from 10-15 14:00Z (K-UP is now a 10-16 s2 candidate, coordinator checklist step 5)
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null \
        || ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > $W/build701.pid
log "build701 armed (pid $$): ku30 = 446bj2z + VCC_RC_CALLSET_CPM=3000; cap $CAP"
gates(){ local bad=0
  eq_head $BU || { echo "builder != HEAD"; bad=1; }
  eq_head $EMP || { echo "emitter != HEAD"; bad=1; }
  [ "$(sha "$DL")" = "$DLS" ] || { echo "LODO dump sha"; bad=1; }
  [ -s "$REFZ" ] || { echo "446bj2z meta missing"; bad=1; }
  [ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"][:16])' "$REFZ")" = ec3a5ae6167ab2d6 ] || { echo "446bj2z meta is not ec3a5ae6"; bad=1; }
  [ "$(disk_gib)" -ge 35 ] || { echo "disk $(disk_gib) GiB < 35"; bad=1; }
  [ "$(fi_gib)" -ge 23 ] || { echo "free+inactive $(fi_gib) GiB < 23"; bad=1; }
  busy && { echo "another build or vcc prep is running"; bad=1; }
  return $bad; }
run(){
  [ -e "$S/$STEM.vcc" ] && { log "NOT READY a package already exists at $S/$STEM.vcc"; return 1; }
  until R="$(gates)"; do
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "NOT READY gates not met by CAP: $(echo "$R" | tr '\n' ';')"; return 1; }
    sleep 60; done
  BUS="$(sha $BU)"; EMS="$(sha $EMP)"
  INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFZ")"
  EXPECT="rc_callset_cpm"; [ -n "$INC_EM" ] && [ "$EMS" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
  log "gates PASS; builder ${BUS:0:16} emitter ${EMS:0:16}; EXPECT vs 446bj2z: $EXPECT"
  rm -f "$AUDIT"
  ( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs
    unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_BLOCK \
          VCC_DE_UPRN VCC_DE_DNMAG VCC_BULK_LOAD_CTX
    export VCC_RC_CALLSET_CPM=3000 VCC_RC_AUDIT_DUMP="$AUDIT"
    bash $E/run_arm_generic.sh "$TAG" 1.75 "$REFZ" kup701 "$STEM" "$BUS" "$EMS" "$EXPECT" ) >> "$W/build_$TAG.nohup" 2>&1
  RL=$E/work/build_$TAG/run$TAG.log
  log "run log: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|zero-gene audit|package sha|novelty rc|DONE$|stopping|GATE" "$RL" 2>/dev/null | tr '\n' ' ' | cut -c1-900)"
  [ "$(sha $BU)" = "$BUS" ] && [ "$(sha $EMP)" = "$EMS" ] || { log "NOT READY builder or emitter moved during the build"; return 1; }
  /usr/bin/grep -q "^\[$TAG\] PASS: exactly" "$RL" 2>/dev/null || { log "NOT READY config diff is not exactly $EXPECT"; return 1; }
  /usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" 2>/dev/null || { log "NOT READY cliff check"; return 1; }
  ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
  echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || { log "NOT READY zero-gene audit: ${ZL:-missing}"; return 1; }
  NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"
  [ "$NOV" = "rc=0" ] || { log "NOT READY novelty ${NOV:-missing}"; return 1; }
  PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
  [ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] && [ "$(sha "$S/$STEM.vcc")" = "$PK" ] || { log "NOT READY no package or sha differs"; return 1; }
  printf '%s READY ku30 %s %s\n' "$(date -u +%FT%TZ)" "$PK" "$S/$STEM.vcc" >> "$X/work/READY-science.log"
  log "READY ku30 $PK"; return 0; }
run || true
log "END"
