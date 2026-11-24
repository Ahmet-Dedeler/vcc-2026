#!/usr/bin/env bash
# 690 X1 build (readiness [401d9b], 2026-10-22; the coordinator's GO, SIXA d98c1c48, after day708's decisions 06:09:50Z):
#   X1 = R0's day recipe (the frozen decisions 04be2af3, re-resolved, + DEF_DAY_CD4_THRESHOLDS=relaxed, which every D/E/F day
#        recipe carries from the start) + day708's FINAL: DEF_DAY_H1_STACK=x3, DEF_DAY_CD4_STACK=x1.5 (weights h3c1.5) and
#        DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh (swap707 ON). The 705 rider stays ON as in R0.
#   Built WITHOUT DEF_DAY_ALLOW_LOW_DOSE. 6a (d), ruled before the build: if postbuild's ONLY failing row is 554's with every
#   lambda >= 0.85, the finished build is re-gated with the override exactly as R0 was (a judge waiver); any lambda < 0.85, any
#   other real FAIL, or no safe build retires X1 (no plan B for a replacement).
# Fresh work dir experiments/692-def-r0/work/x1 (name x1def), r0ready692 on it with OUT692=experiments/692-def-r0/work/x1-ready
# (package submissions/x1def.vcc). r0 stays on disk (X1's 6a diff is judged against r0's meta and npz). Never submits: the
# coordinator arms the submit once X1 is the candidate that ships.
# Gates as run690 / watch691 AUTO_R0: builder == HEAD (e99bbd4a) and the other frozen shas, no builder / operator fit, no
# day-path chain, free+inactive >= 20 GiB and memory_pressure free >= 40 %, disk >= 35 GiB, the hook committed and clean.
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/x1_690.sh > experiments/690-def-readiness-1015/work/x1_690.nohup 2>&1 &
#   DRY=1 bash experiments/690-def-readiness-1015/x1_690.sh     (resolve + gates + print; launches nothing)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
DRY="${DRY:-0}"
B=data/challenge-def; WX=experiments/692-def-r0/work/x1; N=x1def; OUTX=experiments/692-def-r0/work/x1-ready
HOOK=experiments/707-strong-tier-direction/dumphook707.sh
LOG=experiments/690-def-readiness-1015/work/x1_690.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
die(){ log "X1 STOP: $*"; exit 1; }
sha8(){ shasum -a 256 "$1" | cut -c1-8; }
log "x1_690 (pid $$, DRY=$DRY): X1 = R0 recipe + relaxed CD4 + H1 x3 + CD4 x1.5 + hook707 in $WX ($N); 692 -> $OUTX"
# gates
[ "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -c1-8)" = e99bbd4a ] || die "HEAD's builder is not e99bbd4a"
[ "$(sha8 scripts/make_submission_pm_cd4.py)" = e99bbd4a ] || die "builder on disk is not e99bbd4a"
for fs in docs/october-decisions.json:04be2af3 scripts/def_day.sh:11ff102e src/vcc2026/rc_emit_k.py:0f56209d scripts/def_day_helpers.py:5354e449; do
  [ "$(sha8 "${fs%%:*}")" = "${fs##*:}" ] || die "${fs%%:*} is not ${fs##*:}"
done
git ls-files --error-unmatch "$HOOK" >/dev/null 2>&1 && [ -z "$(git status --porcelain -- "$HOOK")" ] || die "$HOOK is not committed and clean"
[ -e "$WX" ] && die "$WX exists (a build goes in a FRESH dir)"
[ -e "$OUTX" ] && die "$OUTX exists"
[ -e "submissions/$N.vcc" ] && die "submissions/$N.vcc exists"
[ -s experiments/692-def-r0/work/r0/r0def_final.h5ad.meta.json ] || die "r0's final meta is missing (keep r0: X1 is judged against it)"
gate(){
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  [ "$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh|day708\.sh)/' | wc -l | tr -d ' ')" = 0 ] || { echo "a day-path chain or day708 is running"; return 1; }
  local f mp d; f="$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)"
  mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"; d="$(df -g /System/Volumes/Data | awk 'NR==2{print $4}')"
  [ "$f" -ge 20 ] || { echo "free+inactive $f GiB < 20"; return 1; }
  [ "${mp:-0}" -ge 40 ] || { echo "memory_pressure free ${mp:-?} % < 40"; return 1; }
  [ "$d" -ge 35 ] || { echo "disk $d GiB < 35"; return 1; }
  echo "no builder, no chain, free+inactive $f GiB, memory_pressure free $mp %, disk $d GiB"
}
r="$(gate)" || die "gate: $r"
log "gate PASS: $r"
# recipe: R0's, re-resolved from the live decisions file, + the X1 switches
RJ=experiments/690-def-readiness-1015/work/def-day/x1-recipe.json
env -u OCTOBER_DECISIONS uv run python scripts/october_recipe.py command --bundle "$B" --work "$WX" --name "$N" --json > "$RJ" || die "october_recipe command failed"
ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in json.load(open(sys.argv[1]))["env"].items()))' "$RJ")"
ARGS="$(python3 -c 'import json,shlex,sys; print(" ".join(shlex.quote(x) for x in json.load(open(sys.argv[1]))["args"]))' "$RJ")"
case "$ARGS" in "$B $WX $N 1.75 0.5 noprog 4") ;; *) die "unexpected recipe args: $ARGS";; esac
R0E="$(python3 -c 'import json,sys; print(sorted(json.load(open(sys.argv[1]))["env"].items()))' experiments/690-def-readiness-1015/work/def-day/r0-recipe-relaxed.json)"
X1E="$(python3 -c 'import json,sys; print(sorted(json.load(open(sys.argv[1]))["env"].items()))' "$RJ")"
[ "$R0E" = "$X1E" ] || die "the recipe env differs from R0's: $X1E"
XTRA="DEF_DAY_CD4_THRESHOLDS=relaxed DEF_DAY_H1_STACK=x3 DEF_DAY_CD4_STACK=x1.5 DEF_DAY_DUMP_HOOK=$HOOK"
CMD="env -u OCTOBER_DECISIONS -u DEF_DAY_ALLOW_LOW_DOSE $XTRA $ENVS nohup bash scripts/detach.sh bash scripts/def_day.sh $ARGS > $WX.nohup 2>&1 &"
log "recipe env == R0's (${#R0E} chars); X1 adds: $XTRA (hook $(sha8 "$HOOK"))"
CMD692="OUT692=$OUTX nohup bash scripts/detach.sh bash experiments/692-def-r0/r0ready692.sh $B $WX $N 202610270000 > $OUTX/r0ready692.nohup 2>&1 &"
if [ "$DRY" = 1 ]; then log "DRY: $CMD"; log "DRY: mkdir -p $OUTX; $CMD692"; log "END"; exit 0; fi
log "RUN: $CMD"; eval "$CMD"
sleep 20; /usr/bin/grep -q "start: bundle" "$WX/def_day.log" 2>/dev/null || die "def_day in $WX did not start (see $WX.nohup)"
/usr/bin/grep -E "H1 stack|CD4 stack|690g|relaxed" "$WX/def_day.log" | head -6 >> "$LOG"
mkdir -p "$OUTX"; log "RUN: $CMD692"; eval "$CMD692"
sleep 10; /usr/bin/grep -q "r0ready692 armed" "$OUTX/$N.run.log" 2>/dev/null || die "r0ready692 on $WX did not arm (see $OUTX/r0ready692.nohup)"
log "launched: def_day in $WX; r0ready692 -> $OUTX/$N.ready.log, $OUTX/$N.run.log (package submissions/$N.vcc)"
log "END"
