#!/usr/bin/env bash
# 704 (science session [bf003b], 2026-10-14): build the CONTINGENCY arm of experiments/704-alpha125-contingency/PREREGISTRATION.md,
# only if DN (671, 10-15 s1) reads NULL/HARM, and only when the coordinator launches it after looking at DN's falsifiers.
#   nohup bash experiments/704-alpha125-contingency/build704.sh > experiments/704-alpha125-contingency/work/build704.nohup 2>&1 &
# ARM 704a125 = 630arm's exact env (dn671's step-1 env without the DN knob) + VCC_ALPHA=1.2275 (carrier x1.25 under BULK_LOAD),
#     built by run_arm_generic.sh against 630arm's meta, EXPECT alpha (+ de_sign_emitter_sha256 if the pin differs).
# Gates: builder and emitter == HEAD, LODO dump sha, disk >= 35 GiB, free+inactive >= 23 GiB, no other build or vcc prep.
# After: "[a704] PASS: exactly", CLIFF OK, zero-gene audit <= 6 gated zero-total cells is the prereg falsifier (production's
#        audit line must read 0 in 0 blocks, as every arm's), novelty rc 0, package sha. Terminal: "READY a704 <sha>" |
#        "NOT READY <reason>", then "END"; READY is appended to experiments/701-kup-callset/work/READY-science.log as
#        "<ts> READY a704 <sha> <pkg>". Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/704-alpha125-contingency; W=$X/work; mkdir -p "$W"; LOG=$W/build704.log
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py; E=experiments/428-zero-gene-audit; S=submissions
TAG=a704
STEM="446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100-a125"
REF630="$REPO/$S/446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
CAP="${A704_CAP:-202610160600}"
READYLOG=experiments/701-kup-callset/work/READY-science.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null \
        || ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > "$W/build704.pid"
TERMINAL=0; trap 'rc=$?; [ "$TERMINAL" = 1 ] || { log "NOT READY unexpected exit rc $rc"; log "END"; }' EXIT
done_(){ TERMINAL=1; log "$1"; log "END"; exit 0; }
log "build704 armed (pid $$): a704 = 630arm + VCC_ALPHA=1.2275; cap $CAP"
[ -s "$REF630" ] || done_ "NOT READY 630arm meta missing"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"][:16])' "$REF630")" = 8b075f7d81c006a0 ] \
  || done_ "NOT READY 630arm's meta is not the scored 630arm (output 8b075f7d)"
[ -e "$S/$STEM.vcc" ] && done_ "NOT READY a package already exists at $S/$STEM.vcc"
until { eq_head $BU && eq_head $EMP && [ "$(sha "$DL")" = "$DLS" ] && [ "$(disk_gib)" -ge 35 ] && [ "$(fi_gib)" -ge 23 ] && ! busy; }; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && done_ "NOT READY gates not met by $CAP"
  sleep 60; done
BUS="$(sha $BU)"; EMS="$(sha $EMP)"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REF630")"
EXPECT="alpha"; [ -n "$INC_EM" ] && [ "$EMS" != "$INC_EM" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
log "gates PASS; builder ${BUS:0:16} emitter ${EMS:0:16}; EXPECT vs 630arm: $EXPECT"
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs \
         VCC_DE_UPRN=b_k100_q20r3p10 VCC_ALPHA=1.2275
  unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK \
        VCC_DE_DNMAG VCC_BULK_LOAD_CTX VCC_RC_AUDIT_DUMP
  bash $E/run_arm_generic.sh "$TAG" 1.75 "$REF630" alpha704 "$STEM" "$BUS" "$EMS" "$EXPECT" ) >> "$W/build_$TAG.nohup" 2>&1
RL=$E/work/build_$TAG/run$TAG.log
log "run log: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|zero-gene audit|package sha|novelty rc|DONE$|stopping|GATE" "$RL" 2>/dev/null | tr '\n' ' ' | cut -c1-900)"
[ "$(sha $BU)" = "$BUS" ] && [ "$(sha $EMP)" = "$EMS" ] || done_ "NOT READY builder or emitter moved during the build"
/usr/bin/grep -q "^\[$TAG\] PASS: exactly" "$RL" 2>/dev/null || done_ "NOT READY config diff is not exactly $EXPECT"
/usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" 2>/dev/null || done_ "NOT READY cliff check"
ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || done_ "NOT READY zero-gene audit: ${ZL:-missing}"
NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"; [ "$NOV" = "rc=0" ] || done_ "NOT READY novelty ${NOV:-missing}"
PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
[ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] && [ "$(sha "$S/$STEM.vcc")" = "$PK" ] || done_ "NOT READY no package or sha differs"
printf '%s READY a704 %s %s\n' "$(date -u +%FT%TZ)" "$PK" "$S/$STEM.vcc" >> "$READYLOG"
done_ "READY a704 $PK"
