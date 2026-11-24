#!/usr/bin/env bash
# 693b — the two build-level known answers after 693p (dn_mag v2 pinned in the builder). The coordinator chains this right
# after applying + committing 693p. Builds strictly serial, behind the same gates as dn671 (builder == HEAD and carrying
# dn_mag2's pin; nothing else building or prepping; free+inactive >= 23 GiB; disk >= 35 GiB).
#   KAu2   all knobs unset on the bumped builder == 446bj2z's ec3a5ae6 (527's run_ka_527.sh, dn671 step 2's env, verbatim)
#   DNv2   dn671's DN build (630arm + VCC_DE_DNMAG=d_m200, the exact step-1 env) on the v2 tree, new tag dn631v2. PASS iff the
#          novelty checker recognises its package as the already-scored DN entry wzaVLFTmlAMW3Mfw0dD8 (ledger 3c6a51b8:
#          content c5e37d4f…, h5ad 62d52937…): the board ledger is the judge. The meta's de_dnmag.module_sha256 moves to v2
#          by design; the prediction must not.
#   bash experiments/693-dnmag-speed/ka693b.sh            (nohup it; terminal line "<ts> END"; verdicts "KAU2 PASS|FAIL",
#                                                           "DNV2 PASS|FAIL")
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/693-dnmag-speed; W=$X/work; mkdir -p "$W"; LOG=$W/ka693b.log
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py; DNM2=src/vcc2026/dn_mag2.py
E=experiments/428-zero-gene-audit; S=submissions; K=experiments/527-arc-hook-arm
STEM0="446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100"
STEM="$STEM0-dnm200v2"; TAG=dn631v2; DN_ENTRY=wzaVLFTmlAMW3Mfw0dD8
REF630="$REPO/$S/$STEM0.h5ad.meta.json"
REFZ="$REPO/$S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
CAP="${CAP693B:-202610161200}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
trap 'log "END"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null \
        || ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)'; }
V2SHA="$(sha $X/dn_mag2.py)"
gates(){ local bad=0
  [ "$(sha $BU)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || { echo "builder != HEAD"; bad=1; }
  /usr/bin/grep -q "from vcc2026 import dn_mag2 as _dn_mag" $BU && /usr/bin/grep -q "$V2SHA" $BU || { echo "builder does not pin dn_mag2 $V2SHA (693p not applied?)"; bad=1; }
  [ -f $DNM2 ] && [ "$(sha $DNM2)" = "$V2SHA" ] || { echo "src/vcc2026/dn_mag2.py missing or not v2"; bad=1; }
  [ "$(sha "$DL")" = "$DLS" ] || { echo "LODO dump sha"; bad=1; }
  [ -s "$REF630" ] && [ -s "$REFZ" ] || { echo "ref meta missing"; bad=1; }
  [ "$(disk_gib)" -ge 35 ] || { echo "disk $(disk_gib) GiB < 35"; bad=1; }
  [ "$(fi_gib)" -ge 23 ] || { echo "free+inactive $(fi_gib) GiB < 23"; bad=1; }
  busy && { echo "another build or vcc prep is running"; bad=1; }
  return $bad; }
waitg(){ local r last=""; until r="$(gates)"; do [ "$r" != "$last" ] && log "waiting: $(echo "$r" | tr '\n' ';')"; last="$r"
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "STOP: gates not met by CAP $CAP"; exit 3; }; sleep 60; done; }
echo $$ > $W/ka693b.pid
POSTB="$(sha $BU)"
log "ka693b armed (pid $$): builder ${POSTB:0:16}, dn_mag2 ${V2SHA:0:16}, cap $CAP"

# ---- KAu2: all unset on the bumped builder (dn671 step 2, verbatim but for the tag and builder sha) ----
waitg
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM \
        VCC_RC_CALLSET_BLOCK VCC_DE_UPRN VCC_DE_DNMAG
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs KA_Q=1.75
  [ -e "$K/work/build_KA1015u2" ] && mv "$K/work/build_KA1015u2" "$K/work/build_KA1015u2.aborted-$(date -u +%Y%m%dT%H%M%SZ)"
  bash $K/run_ka_527.sh KA1015u2 "$REFZ" "$POSTB" ) >> "$W/KA1015u2.nohup" 2>&1
KL=$K/work/build_KA1015u2/runKA1015u2.log
log "KAu2: $(/usr/bin/grep -E 'driver rc=|KA PASS|KA FAIL|builder after restore' $KL 2>/dev/null | tr '\n' ' ' | cut -c1-300)"
if /usr/bin/grep -q "driver rc=0" $KL 2>/dev/null && /usr/bin/grep -q "KA PASS (identical)" $KL; then
  log "KAU2 PASS (== ec3a5ae6 on builder ${POSTB:0:16})"; rm -f $K/work/build_KA1015u2/ka527-KA1015u2.h5ad
else log "KAU2 FAIL (see $KL)"; exit 1; fi

# ---- DNv2: dn671's DN build on the v2 tree; the ledger judges the prediction ----
[ -e "$S/$STEM.vcc" ] && { log "DNV2 FAIL: a package already exists at $S/$STEM.vcc"; exit 1; }
waitg
NEWE="$(sha $EMP)"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REF630")"
EXPECT="de_dnmag"; [ -n "$INC_EM" ] && [ "$NEWE" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
AUDIT="$REPO/$W/$TAG.rc_audit.jsonl"; rm -f "$AUDIT"
t0=$(date +%s)
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs
  unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK
  export VCC_DE_UPRN=b_k100_q20r3p10 VCC_DE_DNMAG=d_m200 VCC_RC_AUDIT_DUMP="$AUDIT"
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REF630" dnmag631v2 "$STEM" "$POSTB" "$NEWE" "$EXPECT" ) >> "$W/build_$TAG.nohup" 2>&1
log "DNv2 build: $(( ($(date +%s) - t0) / 60 )) min; build stage $(python3 -c 'import json,sys; t=json.load(open(sys.argv[1])); t=t.get("stages", t); v=t["build"]; print(v.get("seconds", v.get("wall_s")) if isinstance(v, dict) else v)' "$E/work/build_$TAG/timings.json" 2>/dev/null || echo '?') s (v1: 9,727.9 s)"
[ -s "$S/$STEM.vcc" ] || { log "DNV2 FAIL: no package"; exit 1; }
NOV="$(uv run python scripts/check_submission_novelty.py check "$S/$STEM.vcc" 2>&1)"
printf '%s\n' "$NOV" >> "$LOG"
if echo "$NOV" | /usr/bin/grep -q "ALREADY BEEN SCORED" && echo "$NOV" | /usr/bin/grep -q "$DN_ENTRY"; then
  log "DNV2 PASS: the v2 tree rebuilt DN's scored prediction ($(echo "$NOV" | /usr/bin/grep -m1 'matched' | cut -c1-140)), entry $DN_ENTRY"
  rm -f "$S/$STEM.vcc" "$S/$STEM.h5ad"
else
  log "DNV2 FAIL: the rebuilt prediction is not DN's ($(echo "$NOV" | tail -2 | tr '\n' ' ' | cut -c1-200))"; exit 1
fi
