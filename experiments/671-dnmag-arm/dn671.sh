#!/usr/bin/env bash
# 671 (coordinator [400b6e], 2026-10-14 ~22:45Z): the DNMAG arm of 631's FROZEN preregistration
# (experiments/631-nmae-headroom/PREREGISTRATION-DNMAG.md), then its KA0 (all-unset rebuild on the patched builder).
#   nohup bash experiments/671-dnmag-arm/dn671.sh > experiments/671-dnmag-arm/work/dn671.nohup 2>&1 &
# The builder carries patch_builder631 (applied by hand 22:40Z, uncommitted until KAu passes): POSTB below.
# Step 1  DN = 630arm's exact env + VCC_DE_DNMAG=d_m200, built by run_arm_generic.sh against 630arm's meta, EXPECT de_dnmag.
#         Gates after: "[dn631] PASS: exactly", CLIFF OK, zero-gene audit 0, novelty rc 0, package sha; then the audit's
#         dn fail-safe counts (dn_invariant_failed / dn_guard_failed must be 0 wherever they appear).
# Step 2  KAu (KA1015u): all-unset on the patched builder == 446bj2z (ec3a5ae6), the byte-identity proof for the unset knob.
# Terminal lines (always one of each step): READY dn671 <pkg sha> | NOT READY <reason>;  KAU PASS | KAU FAIL <why>; then END.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/671-dnmag-arm; W=$X/work; mkdir -p "$W"; LOG=$W/dn671.log
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py; DNM=src/vcc2026/dn_mag.py
POSTB=8f8f6d5c038c53d08f9bc0189a5084ad9df71b2a425cb6f8401dc461193509a8
DNM_SHA=0a8cc560231a16f0fff7fbbf316e3a5e190cbdf20490239ef354b49b3ff60dac
E=experiments/428-zero-gene-audit; S=submissions; K=experiments/527-arc-hook-arm
TAG=dn631
STEM="446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100-dnm200"
REF630="$REPO/$S/446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100.h5ad.meta.json"
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
AUDIT="$REPO/$W/$TAG.rc_audit.jsonl"
CAP=202610151200
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null \
        || ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > $W/dn671.pid
log "dn671 armed (pid $$): DN = 630arm + VCC_DE_DNMAG=d_m200, then KAu; cap $CAP"
gates(){ local bad=0
  [ "$(sha $BU)" = "$POSTB" ] || { echo "builder $(sha $BU | cut -c1-16) != POSTB"; bad=1; }
  eq_head $EMP || { echo "emitter != HEAD"; bad=1; }
  /usr/bin/grep -q "631: DOWN magnitude" $BU || { echo "builder lacks 631's mark"; bad=1; }
  [ -f $DNM ] && [ "$(sha $DNM)" = "$DNM_SHA" ] || { echo "dn_mag.py missing or not 0a8cc560"; bad=1; }
  [ "$(sha "$DL")" = "$DLS" ] || { echo "LODO dump sha"; bad=1; }
  [ -s "$REF630" ] && [ -s "$REFZ" ] || { echo "ref meta missing"; bad=1; }
  [ "$(disk_gib)" -ge 35 ] || { echo "disk $(disk_gib) GiB < 35"; bad=1; }
  [ "$(fi_gib)" -ge 23 ] || { echo "free+inactive $(fi_gib) GiB < 23"; bad=1; }
  busy && { echo "another build or vcc prep is running"; bad=1; }
  return $bad; }
# ---- step 1: DN ----
step1(){
  [ -e "$S/$STEM.vcc" ] && { log "NOT READY a package already exists at $S/$STEM.vcc"; return 1; }
  until R="$(gates)"; do
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "NOT READY gates not met by CAP: $(echo "$R" | tr '\n' ';')"; return 1; }
    sleep 60; done
  NEWE="$(sha $EMP)"
  INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REF630")"
  EXPECT="de_dnmag"; [ -n "$INC_EM" ] && [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
  log "step 1 gates PASS; builder ${POSTB:0:16} emitter ${NEWE:0:16}; EXPECT vs 630arm: $EXPECT"
  rm -f "$AUDIT"
  ( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs
    unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK
    export VCC_DE_UPRN=b_k100_q20r3p10 VCC_DE_DNMAG=d_m200 VCC_RC_AUDIT_DUMP="$AUDIT"
    bash $E/run_arm_generic.sh "$TAG" 1.75 "$REF630" dnmag631 "$STEM" "$POSTB" "$NEWE" "$EXPECT" ) >> "$W/build_$TAG.nohup" 2>&1
  RL=$E/work/build_$TAG/run$TAG.log
  log "run log: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|zero-gene audit|package sha|novelty rc|DONE$|stopping|GATE" "$RL" 2>/dev/null | tr '\n' ' ' | cut -c1-900)"
  [ "$(sha $BU)" = "$POSTB" ] && [ "$(sha $EMP)" = "$NEWE" ] || { log "NOT READY builder or emitter moved during the build"; return 1; }
  /usr/bin/grep -q "^\[$TAG\] PASS: exactly" "$RL" 2>/dev/null || { log "NOT READY config diff is not exactly $EXPECT"; return 1; }
  /usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" 2>/dev/null || { log "NOT READY cliff check"; return 1; }
  ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
  echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || { log "NOT READY zero-gene audit: ${ZL:-missing}"; return 1; }
  NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"
  [ "$NOV" = "rc=0" ] || { log "NOT READY novelty ${NOV:-missing}"; return 1; }
  PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
  [ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] && [ "$(sha "$S/$STEM.vcc")" = "$PK" ] || { log "NOT READY no package or sha differs"; return 1; }
  # dn fail-safes: count any nonzero dn_invariant_failed / dn_guard_failed in the audit dump and the meta
  uv run python - "$AUDIT" "$REPO/$S/$STEM.h5ad.meta.json" > "$W/dn_failsafe.out" 2>&1 <<'PY'
import json, sys, re
bad = 0; seen = 0
for p in sys.argv[1:]:
    try: txt = open(p).read()
    except OSError: continue
    for k in ("dn_invariant_failed", "dn_guard_failed"):
        for m in re.finditer(r'"%s"\s*:\s*([0-9.]+)' % k, txt):
            seen += 1; bad += float(m.group(1)) != 0
print("dn fail-safe fields seen", seen, "nonzero", bad)
sys.exit(1 if bad else 0)
PY
  FRC=$?; log "$(cat $W/dn_failsafe.out | tail -1)"
  [ $FRC -eq 0 ] || { log "NOT READY dn fail-safe triggered (see $W/dn_failsafe.out)"; return 1; }
  log "READY dn671 $PK"; return 0; }
step1 || true
# ---- step 2: KAu on the patched builder ----
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM \
        VCC_RC_CALLSET_BLOCK VCC_DE_UPRN VCC_DE_DNMAG
  until R="$(gates)"; do [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && exit 3; sleep 60; done
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs KA_Q=1.75
  [ -e "$K/work/build_KA1015u" ] && mv "$K/work/build_KA1015u" "$K/work/build_KA1015u.aborted-$(date -u +%Y%m%dT%H%M%SZ)"
  bash $K/run_ka_527.sh KA1015u "$REFZ" "$POSTB" ) >> "$W/KA1015u.nohup" 2>&1
KL=$K/work/build_KA1015u/runKA1015u.log
log "KAu: $(/usr/bin/grep -E 'driver rc=|KA PASS|KA FAIL|builder after restore' $KL 2>/dev/null | tr '\n' ' ' | cut -c1-300)"
if /usr/bin/grep -q "driver rc=0" $KL 2>/dev/null && /usr/bin/grep -q "KA PASS (identical)" $KL; then log "KAU PASS (== ec3a5ae6 on builder ${POSTB:0:16})"
  rm -f $K/work/build_KA1015u/ka527-KA1015u.h5ad
else log "KAU FAIL (see $KL)"; fi
[ "$(sha $BU)" = "$POSTB" ] || log "WARNING builder after KAu $(sha $BU | cut -c1-16) != POSTB (pin not restored?)"
log "END"
