#!/usr/bin/env bash
# 690 chain_x1x2 (readiness [401d9b], 2026-10-22): carries the two 6a replacement candidates to a judged READY without a session.
# The coordinator's rulings: X1 GO (SIXA d98c1c48), X2 = C0 (SIXA 1ebd3a1f), one builder at a time, 6a (d) for both:
# relaxed CD4 from the start; a postbuild FAIL whose ONLY row is 554's with every lambda >= 0.85 is re-gated on the finished build
# with DEF_DAY_ALLOW_LOW_DOSE=1 (regate_r0.sh, as R0); any other FAIL or no safe build retires the candidate (no retry).
#   X1 phase: wait for X1's def_day (experiments/692-def-r0/work/x1, launched by x1_690.sh) and its r0ready692 (x1-ready) to end;
#             NOT READY on a postbuild STOP -> regate_r0.sh into x1g-ready (it refuses anything but a 554-only FAIL);
#             READY -> gate6a_x1.py --candidate x1. Logs "X1 RESULT: <READY dir> 6A PASS|FAIL" or "X1 RETIRED: <why>".
#   X2 phase: waits for experiments/670-coord-1015/work/GO-X2 (the coordinator: LZ PASS and a hook that covers D/E/F's panel);
#             NOGO-X2 at any time before the launch cancels it. Then:
#             preflight: the hook (committed and clean) on r0's stack dump into a fresh reference dir must exit 0 (it refused
#             D/E/F at 07:0xZ: 136 targets missing from 811's banks); build gates as x1_690; def_day in experiments/692-def-r0/work/x2
#             (x2def) with R0's recipe (asserted equal) + relaxed CD4 + DEF_DAY_DUMP_HOOK=dumphook811c0.sh, ALLOW_LOW_DOSE unset;
#             r0ready692 -> x2-ready; the same 554 re-gate path into x2g-ready; gate6a --candidate x2 --hook-ref <reference>.
#             Logs "X2 RESULT: …" or "X2 CANCELLED/RETIRED: …".
# Never submits and never arms: the coordinator arms the candidate that ships (X2 if eligible and passing, else X1, else R0).
# Cap CAP_X1X2 (default 202610240000). Terminal line, always: "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/chain_x1x2.sh > experiments/690-def-readiness-1015/work/chain_x1x2.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015
if [ -z "${CHAIN_X1X2_FROZEN:-}" ]; then
  FZ="$E/work/.frozen-$$-chain_x1x2.sh"; cp "$E/chain_x1x2.sh" "$FZ" || exit 1
  CHAIN_X1X2_FROZEN=1 exec bash "$FZ" "$@"
fi
LOG="${LOG_X1X2:-$E/work/chain_x1x2.log}"; CAP="${CAP_X1X2:-202610240000}"
B=data/challenge-def; R0W=experiments/692-def-r0/work/r0; CW=experiments/670-coord-1015/work
HOOK2=experiments/811-lolo-model/dumphook811c0.sh
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"; rm -f "$E/work/.frozen-$$-chain_x1x2.sh"' EXIT
now(){ date -u +%Y%m%d%H%M; }
sha8(){ shasum -a 256 "$1" | cut -c1-8; }
dd_alive(){ ps -axo comm=,args= | awk -v p="scripts/def_day.sh $B $1 " '$1 ~ /bash$/ && index($0, p)' | /usr/bin/grep -q .; }
rr_alive(){ ps -axo comm=,args= | awk -v p="r0ready692.sh $B $1 " '$1 ~ /bash$/ && index($0, p)' | /usr/bin/grep -q .; }
echo $$ > "$E/work/chain_x1x2.pid"
log "chain_x1x2 armed (pid $$): X1 to its judged READY, then X2 on GO-X2 ($CW/GO-X2; NOGO-X2 cancels); cap $CAP"

# settle <work> <stem> <ready dir> <regate dir>: wait for def_day + r0ready692, re-gate a postbuild STOP; echo the READY dir or "".
settle(){
  local W=$1 N=$2 RD=$3 GD=$4 tag=$5
  until ! dd_alive "$W" && /usr/bin/grep -q "builder after restore:" "$W/def_day.log" 2>/dev/null; do
    [ "$(now)" -ge "$CAP" ] && { log "$tag: def_day not finished by the cap"; return 1; }; sleep 60; done
  log "$tag def_day finished: $(/usr/bin/grep -E ' rho:|CONVERGED|NOT converged|STOP:|final:' "$W/def_day.log" | tail -4 | tr '\n' ' ' | cut -c1-500)"
  until ! rr_alive "$W" && /usr/bin/grep -qE "READY $N [0-9a-f]{64}|NOT READY" "$RD/$N.ready.log" 2>/dev/null; do
    [ "$(now)" -ge "$CAP" ] && { log "$tag: r0ready692 not finished by the cap"; return 1; }; sleep 30; done
  if /usr/bin/grep -qE "READY $N [0-9a-f]{64}" "$RD/$N.ready.log" && ! /usr/bin/grep -q "NOT READY" "$RD/$N.ready.log"; then
    log "$tag READY in $RD: $(tail -1 "$RD/$N.ready.log" | cut -c1-120)"; SETTLED="$RD"; return 0
  fi
  log "$tag NOT READY in $RD: $(/usr/bin/grep 'NOT READY' "$RD/$N.ready.log" | tail -1 | cut -c1-300)"
  if /usr/bin/grep -q "STOP: postbuild FAILED" "$W/def_day.log"; then
    log "$tag: postbuild STOP -> regate_r0.sh into $GD (it refuses anything but a 554-only FAIL with every lambda >= 0.85)"
    REGATE_W="$W" REGATE_N="$N" REGATE_OUT="$GD" REGATE_KA_SHA=meta bash "$E/regate_r0.sh" > "$E/work/regate_$tag.nohup" 2>&1
    if /usr/bin/grep -qE "READY $N [0-9a-f]{64}" "$GD/$N.ready.log" 2>/dev/null && ! /usr/bin/grep -q "NOT READY" "$GD/$N.ready.log"; then
      log "$tag READY after the 554 re-gate in $GD: $(tail -1 "$GD/$N.ready.log" | cut -c1-120)"; SETTLED="$GD"; return 0
    fi
    log "$tag RETIRED: the re-gate did not write READY: $(/usr/bin/grep 'NOT READY' "$GD/$N.ready.log" 2>/dev/null | tail -1 | cut -c1-300)"; return 1
  fi
  log "$tag RETIRED: NOT READY on something other than a postbuild STOP (6a (d): no retry)"; return 1
}

# ---- X1
SETTLED=""
if settle experiments/692-def-r0/work/x1 x1def experiments/692-def-r0/work/x1-ready experiments/692-def-r0/work/x1g-ready X1; then
  uv run python "$E/gate6a_x1.py" --candidate x1 --work experiments/692-def-r0/work/x1 --stem x1def --ready-dir "$SETTLED" > "$E/work/gate6a_x1.out" 2>&1
  log "X1 RESULT: READY dir $SETTLED; $(tail -1 "$E/work/gate6a_x1.out" | cut -c1-600)"
fi

# ---- X2
last=""
until [ -f "$CW/GO-X2" ] || [ -f "$CW/NOGO-X2" ]; do
  [ "$last" = w ] || log "X2: waiting for $CW/GO-X2 (or NOGO-X2)"; last=w
  [ "$(now)" -ge "$CAP" ] && { log "X2 CANCELLED: no GO-X2 by the cap"; exit 0; }; sleep 60; done
[ -f "$CW/NOGO-X2" ] && { log "X2 CANCELLED: NOGO-X2 ($(head -c 200 "$CW/NOGO-X2"))"; exit 0; }
log "X2: GO-X2 ($(head -c 200 "$CW/GO-X2"))"
git ls-files --error-unmatch "$HOOK2" >/dev/null 2>&1 && [ -z "$(git status --porcelain -- "$HOOK2" experiments/811-lolo-model/hook811.py)" ] \
  || { log "X2 CANCELLED: $HOOK2 or hook811.py is not committed and clean"; exit 1; }
REF="$E/work/def-day/hook811c0_ref-$(date -u +%H%M%S)"; mkdir -p "$REF"
bash "$HOOK2" "$B" "$R0W" "$R0W/dump_h1x2_cd4x1_s0.25.npz" "$R0W/identity.json" "$REF" > "$REF/hook.log" 2>&1; HR=$?
[ "$HR" = 0 ] && [ -s "$REF/dump_h1x2_cd4x1_s0.25.npz" ] \
  || { log "X2 CANCELLED: the hook fails its preflight on r0's stack dump (rc $HR): $(/usr/bin/grep -E 'Error|refuse' "$REF/hook.log" | tail -1 | cut -c1-300)"; exit 1; }
log "X2 preflight: the hook ($(sha8 "$HOOK2"), hook811 $(sha8 experiments/811-lolo-model/hook811.py)) on r0's stack dump -> $REF ($(sha8 "$REF/dump_h1x2_cd4x1_s0.25.npz"))"
WX=experiments/692-def-r0/work/x2; N=x2def; OUTX=experiments/692-def-r0/work/x2-ready
for p in "$WX" "$OUTX" experiments/692-def-r0/work/x2g-ready "submissions/$N.vcc"; do [ -e "$p" ] && { log "X2 CANCELLED: $p exists (fresh dirs only)"; exit 1; }; done
[ "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -c1-8)" = e99bbd4a ] || { log "X2 CANCELLED: HEAD's builder is not e99bbd4a"; exit 1; }
for fs in docs/october-decisions.json:04be2af3 scripts/def_day.sh:11ff102e src/vcc2026/rc_emit_k.py:0f56209d scripts/def_day_helpers.py:5354e449; do
  [ "$(sha8 "${fs%%:*}")" = "${fs##*:}" ] || { log "X2 CANCELLED: ${fs%%:*} is not ${fs##*:}"; exit 1; }
done
gate(){
  [ -f "$CW/NOGO-X2" ] && { echo "NOGO-X2"; return 2; }
  [ "$(sha8 scripts/make_submission_pm_cd4.py)" = e99bbd4a ] || { echo "builder on disk is not e99bbd4a"; return 1; }
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  [ "$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh|day708\.sh|regate_r0\.sh|r0ready692\.sh)/' | wc -l | tr -d ' ')" = 0 ] || { echo "a day-path chain, day708, a re-gate or a 692 is running"; return 1; }
  [ "$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc (prep|submit)( |$)' | tr -d ' ')" = 0 ] || { echo "a vcc prep / submit is running"; return 1; }
  local f mp d; f="$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)"
  mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"; d="$(df -g /System/Volumes/Data | awk 'NR==2{print $4}')"
  [ "$f" -ge 20 ] && [ "${mp:-0}" -ge 40 ] && [ "$d" -ge 35 ] || { echo "free+inactive $f GiB, mp ${mp:-?} %, disk $d GiB"; return 1; }
  echo "free+inactive $f GiB, mp $mp %, disk $d GiB"
}
last=""; while :; do r="$(gate)"; g=$?; [ "$g" = 0 ] && break
  [ "$g" = 2 ] && { log "X2 CANCELLED: NOGO-X2 before the launch"; exit 0; }
  [ "$r" = "$last" ] || log "X2 waiting: $r"; last="$r"
  [ "$(now)" -ge "$CAP" ] && { log "X2 CANCELLED: build gate not open by the cap"; exit 0; }; sleep 60; done
log "X2 gate PASS: $r"
RJ="$E/work/def-day/x2-recipe.json"
env -u OCTOBER_DECISIONS uv run python scripts/october_recipe.py command --bundle "$B" --work "$WX" --name "$N" --json > "$RJ" || { log "X2 CANCELLED: october_recipe failed"; exit 1; }
[ "$(python3 -c 'import json,sys; print(sorted(json.load(open(sys.argv[1]))["env"].items()))' "$RJ")" = \
  "$(python3 -c 'import json,sys; print(sorted(json.load(open(sys.argv[1]))["env"].items()))' "$E/work/def-day/r0-recipe-relaxed.json")" ] \
  || { log "X2 CANCELLED: the recipe env differs from R0's"; exit 1; }
ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in json.load(open(sys.argv[1]))["env"].items()))' "$RJ")"
ARGS="$(python3 -c 'import json,shlex,sys; print(" ".join(shlex.quote(x) for x in json.load(open(sys.argv[1]))["args"]))' "$RJ")"
[ "$ARGS" = "$B $WX $N 1.75 0.5 noprog 4" ] || { log "X2 CANCELLED: unexpected recipe args $ARGS"; exit 1; }
CMD="env -u OCTOBER_DECISIONS -u DEF_DAY_ALLOW_LOW_DOSE DEF_DAY_CD4_THRESHOLDS=relaxed DEF_DAY_DUMP_HOOK=$HOOK2 $ENVS nohup bash scripts/detach.sh bash scripts/def_day.sh $ARGS > $WX.nohup 2>&1 &"
log "X2 RUN: $CMD"; eval "$CMD"
sleep 20; /usr/bin/grep -q "start: bundle" "$WX/def_day.log" 2>/dev/null || { log "X2 RETIRED: def_day did not start (see $WX.nohup)"; exit 1; }
mkdir -p "$OUTX"
OUT692="$OUTX" nohup bash scripts/detach.sh bash experiments/692-def-r0/r0ready692.sh "$B" "$WX" "$N" 202610270000 > "$OUTX/r0ready692.nohup" 2>&1 &
sleep 10; log "X2 launched: def_day in $WX; r0ready692 -> $OUTX"
SETTLED=""
if settle "$WX" "$N" "$OUTX" experiments/692-def-r0/work/x2g-ready X2; then
  uv run python "$E/gate6a_x1.py" --candidate x2 --work "$WX" --stem "$N" --ready-dir "$SETTLED" --hook-ref "$REF/dump_h1x2_cd4x1_s0.25.npz" > "$E/work/gate6a_x2.out" 2>&1
  log "X2 RESULT: READY dir $SETTLED; $(tail -1 "$E/work/gate6a_x2.out" | cut -c1-600)"
fi
