#!/usr/bin/env bash
# 707 (science session [d80a10], 2026-10-15): build the strong-tier direction arm of
# experiments/707-strong-tier-direction/PREREGISTRATION.md (coordinator: 10-17 s2 or 10-18 s1).
#   launched by experiments/707-strong-tier-direction/sci1016.sh (the science chain, detached)
# ARM = ALLSW (658allsw) with ONE change: VCC_BULK_LOAD = the cap-6 dump with the top z_src tertile's rows blended toward
#       the production donors' unit mean, norm-preserving, beta 0.25 (swap707.py; KA-beta0 reproduces ALLSW's dump byte for
#       byte; carriers untouched; P2g tercile edges bitwise unchanged).
#       Re-levelled as 702 / 498: r1 at ALLSW's factors -> ntclevel417 -> r2 at factors x EV (rho within 0.5 % of 0.9793).
# EXPECT vs ALLSW's meta: r1 bulk_load,bulk_load_sha256; r2 + bulk_ctx_scale (+ de_sign_emitter_sha256 if the pin moved).
# START: not before START (default 10-15 18:00Z) and only once none of today's chains is alive (K-UP's build701, the
#       rehearsals run690 / rabc694, ka693b, def_day, R7a's wrapper); then the usual gates (builder and emitter == HEAD,
#       no other build, disk >= 35 GiB, free+inactive >= 23 GiB). CAP for starting a build: 10-16 18:00Z.
# Terminal: "<ts> READY <tag> <pkg sha>" | "<ts> NOT READY <reason>", then "<ts> END". READY lines are also appended to
# experiments/701-kup-callset/work/READY-science.log (the coordinator's arm gate). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/707-strong-tier-direction; W=$X/work; mkdir -p "$W"; E=experiments/428-zero-gene-audit; S=submissions
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
TAG=sw707; DUMP=$X/work/dump_h1x2_cd4x15_lodo05hi6_sw025.npz; PINNED=77e9de2de8d04947de39e9eb4c43fd2bccd2f4029666ac60dd9bc4cae43d4f45
DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"
LOG=$W/build_$TAG.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
REFA="$REPO/$S/658allsw-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5-zs250-hk005.h5ad.meta.json"
QV="tercile:2.5,2.0,1.75"
CAP="${B707_CAP:-202610181200}"; START="${B707_START:-202610152100}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null \
        || ps -axo args= | /usr/bin/grep -vE 'claude|grep' | /usr/bin/grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > "$W/build_$TAG.pid"
TERMINAL=0; trap 'rc=$?; [ "$TERMINAL" = 1 ] || { log "NOT READY unexpected exit rc $rc"; log "END"; }' EXIT
done_(){ TERMINAL=1; log "$1"; log "END"; exit 0; }
log "build707 armed (pid $$): $TAG = ALLSW + VCC_BULK_LOAD=$DUMP (${DSHA:0:16}); start $START; cap $CAP"
[ -s "$REFA" ] || done_ "NOT READY ALLSW meta missing"
[ "$DSHA" = "$PINNED" ] || done_ "NOT READY dump sha ${DSHA:0:16} != pinned ${PINNED:0:16}"
F0="$(python3 -c 'import json,sys; s=json.load(open(sys.argv[1]))["config"]["bulk_ctx_scale"]; print(",".join(f"{k}={s[k]:.4f}" for k in "ABC"))' "$REFA")"
INC_EM="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$REFA")"
gates(){ local bad=0
  eq_head $BU || { echo "builder != HEAD"; bad=1; }
  eq_head $EMP || { echo "emitter != HEAD"; bad=1; }
  [ "$(disk_gib)" -ge 35 ] || { echo "disk $(disk_gib) GiB < 35"; bad=1; }
  [ "$(fi_gib)" -ge 23 ] || { echo "free+inactive $(fi_gib) GiB < 23"; bad=1; }
  busy && { echo "another build or vcc prep is running"; bad=1; }
  return $bad; }
build(){ local t=$1 fact=$2 nm=$3 stem=$4 expect=$5
  until R="$(gates)"; do
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && return 9
    sleep 60; done
  local B EM; B="$(sha $BU)"; EM="$(sha $EMP)"
  [ -n "$INC_EM" ] && [ "$EM" != "$INC_EM" ] && expect="$expect,de_sign_emitter_sha256"
  log "$t: gates PASS; builder ${B:0:16} emitter ${EM:0:16}; factors $fact; EXPECT $expect"
  ( export VCC_BULK_LOAD="$REPO/$DUMP:$DSHA" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$fact" \
           VCC_DE_MAG=dn_k250_zs VCC_DE_MDE_SIGN=clean4:0.05
    unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM \
          VCC_RC_CALLSET_BLOCK VCC_DE_UPRN VCC_DE_DNMAG VCC_BULK_LOAD_CTX
    bash $E/run_arm_generic.sh "$t" "$QV" "$REFA" "$nm" "$stem" "$B" "$EM" "$expect" ) >> "$W/build_$t.nohup" 2>&1
  RL=$E/work/build_$t/run$t.log
  log "$t: $(/usr/bin/grep -E "^\[$t\] (PASS|FAIL|CLIFF)|zero-gene audit|package sha|novelty rc|stopping" "$RL" 2>/dev/null | tr '\n' ' ' | cut -c1-700)"
  [ "$(sha $BU)" = "$B" ] && [ "$(sha $EMP)" = "$EM" ] || { log "$t: builder or emitter moved during the build"; return 1; }
  /usr/bin/grep -q "^\[$t\] PASS: exactly" "$RL" 2>/dev/null; }
level(){ uv run python experiments/417-bulk-sharpen/ntclevel417.py "$1" --bundle data/challenge --contexts A,B,C --sd 0.01 >> "$LOG" 2>&1; }
STEM0="707$TAG-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6sw025load"
# ---- start conditions: START reached and none of today's chains alive (bash + script path, never pgrep -f)
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
# 10-16 07:2xZ ([d80a10], coordinator-approved): R7a resumes right after run690 and before sw707, so its waiter is a chain
# here too (alive from its arming until R7a's wrapper ends; no start race with sci1016's step 6).
CHAINS="experiments/701-kup-callset/build701.sh run690.sh rabc694.sh experiments/693-dnmag-speed/ka693b.sh scripts/def_day.sh experiments/552-r7-reserve/work/r7a_access_1015b.sh experiments/705-lodo-strength-safety/build705.sh after690.sh experiments/552-r7-reserve/work/r7a_after690.sh"
while :; do
  NOW="$(date -u +%Y%m%d%H%M)"
  [ "$NOW" -ge "$CAP" ] && done_ "NOT READY start conditions not met by $CAP"
  if [ "$NOW" -ge "$START" ]; then
    LIVE=""; for c in $CHAINS; do alive "$c" && LIVE="$LIVE $c"; done
    [ -z "$LIVE" ] && break
    [ "${LASTLIVE:-}" = "$LIVE" ] || log "waiting for:$LIVE"; LASTLIVE="$LIVE"
  fi
  sleep 60
done
log "start conditions hold"
# ---- round 1 at ALLSW's factors
T1="${TAG}r1"; S1="$STEM0-ctxlevel-mde010coarse-zr5-zs250-hk005"
build "$T1" "$F0" "$T1" "$S1" "bulk_load,bulk_load_sha256"; rc=$?
[ $rc = 9 ] && done_ "NOT READY gates not met by $CAP (r1)"
[ $rc = 0 ] || done_ "NOT READY r1 build failed or config diff not exact"
H1="$E/work/build_$T1/$T1.h5ad"; level "$H1" || done_ "NOT READY r1 ntclevel failed"
EV="$(/usr/bin/grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' "$LOG" | tail -1 | awk -F'= ' '{print $2}')"
[ -n "$EV" ] || done_ "NOT READY r1 ntclevel printed no factors"
F1="$(python3 -c "import sys; f0={k:float(v) for k,v in (p.split('=') for p in sys.argv[1].split(','))}; ev={k:float(v) for k,v in (p.split('=') for p in sys.argv[2].split(','))}; print(','.join(f'{k}={f0[k]*ev[k]:.4f}' for k in 'ABC'))" "$F0" "$EV")"
log "r1 rho: $(uv run python scripts/def_day_helpers.py rho "${H1%.h5ad}.ntclevel.json" 2>&1 | tail -1); EV $EV -> r2 factors $F1"
rm -f "$H1" "$REPO/$S/$S1.vcc"
# ---- round 2 (and at most one more) at factors x EV until every context's rho is within 0.5 % of 0.9793
inband(){ python3 -c "import sys,re; t=sys.argv[1]; v=[float(x) for x in re.findall(r'[ABC] ([0-9.]+)', t)]; sys.exit(0 if len(v)==3 and all(abs(x/0.9793-1)<=0.005 for x in v) else 1)" "$1"; }
FCUR="$F1"
for ROUND in 2 3; do
  T2="$TAG"; [ "$ROUND" = 3 ] && T2="${TAG}r3"
  S2="$STEM0-ctxlevelr-mde010coarse-zr5-zs250-hk005"; [ "$ROUND" = 3 ] && S2="$S2-r3"
  [ -e "$S/$S2.vcc" ] && done_ "NOT READY a package already exists at $S/$S2.vcc"
  build "$T2" "$FCUR" "$T2" "$S2" "bulk_load,bulk_load_sha256,bulk_ctx_scale"; rc=$?
  [ $rc = 9 ] && done_ "NOT READY gates not met by $CAP (round $ROUND)"
  [ $rc = 0 ] || done_ "NOT READY round $ROUND build failed or config diff not exact"
  RL=$E/work/build_$T2/run$T2.log
  /usr/bin/grep -q "^\[$T2\] CLIFF OK" "$RL" || done_ "NOT READY round $ROUND cliff check"
  ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
  echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || done_ "NOT READY zero-gene audit: ${ZL:-missing}"
  NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"; [ "$NOV" = "rc=0" ] || done_ "NOT READY novelty ${NOV:-missing}"
  H2="$E/work/build_$T2/$T2.h5ad"; level "$H2" || done_ "NOT READY round $ROUND ntclevel failed"
  RHO="$(uv run python scripts/def_day_helpers.py rho "${H2%.h5ad}.ntclevel.json" 2>&1 | tail -1)"
  log "round $ROUND rho: $RHO (must be within 0.5 % of 0.9793)"
  PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
  [ -s "$S/$S2.vcc" ] && [ -n "$PK" ] && [ "$(sha "$S/$S2.vcc")" = "$PK" ] || done_ "NOT READY no package or sha differs"
  if inband "$RHO"; then
    printf '%s READY %s %s %s\n' "$(date -u +%FT%TZ)" "$TAG" "$PK" "$S/$S2.vcc" >> "$READYLOG"
    done_ "READY $TAG $PK"
  fi
  [ "$ROUND" = 3 ] && done_ "NOT READY rho outside the band after round 3 ($RHO)"
  EV="$(/usr/bin/grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' "$LOG" | tail -1 | awk -F'= ' '{print $2}')"
  FCUR="$(python3 -c "import sys; f0={k:float(v) for k,v in (p.split('=') for p in sys.argv[1].split(','))}; ev={k:float(v) for k,v in (p.split('=') for p in sys.argv[2].split(','))}; print(','.join(f'{k}={f0[k]*ev[k]:.4f}' for k in 'ABC'))" "$FCUR" "$EV")"
  log "round $ROUND outside the band: EV $EV -> round 3 factors $FCUR"
  rm -f "$H2" "$REPO/$S/$S2.vcc"
done
