#!/usr/bin/env bash
# 658b: build 524dr2 (experiments/658-prereg-1007-builds/PREREGISTRATION-524dr2.md). It runs recipe_dr524.sh's procedure on
# the applied post-10-06 tree, in its own work dir and stem. It resolves the October recipe at start, then branches against
# 524dr's recipe (experiments/524-october-recipe-dr/work/recipe524.json):
#   N0 (nothing new): build, then a LOCAL --exact content check vs 524dr's package; NEVER a slot (identical content)
#   N1 (only DE switches new: mde_sign): build; pre-slot bulk identity vs 524dr (bulkid658, sums) must hold; gate524
#   N2 / N3 (tercilec and/or a wired H1 dose; N3 = with a DE switch too): build; gate524
# Any other recipe change (a dropped switch, another knob) stops: amend 658b first. Never submits; unarmed.
#
# The LAST line of the log is always ONE terminal line:
#   <ts> READY 524dr2 <package sha256>
#   <ts> NOT READY <reason>
# The exact exports and def_day arguments are recorded in work/dr2/env.json.
#
# usage: DR2_BUILDER_SHA=<post-10-06 builder sha> DR2_EMITTER_SHA=<post-10-06 emitter sha> \
#          [DR2_FLOOR=YYYYmmddHHMM] [DR2_CAP=YYYYmmddHHMM] bash experiments/658-prereg-1007-builds/recipe_dr524dr2.sh
#   DRY=1: resolve the recipe, print the branch and the plan; build nothing.
# Named recipe_dr524dr2.sh (not chain5NN.sh) so no build gate's chain pattern mistakes it for itself (AGENTS 6.2); its own
# "other build" detector excludes its pid. bash 3.2 / BSD grep safe; waits capped by DR2_CAP.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/658-prereg-1007-builds; W=$X/work/dr2; LOG=$W/recipe_dr524dr2.log; R453=experiments/453-abc-dress-rehearsal/work/run
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/524dr2-october-recipe-abc; NAME=524dr2
PREV_RECIPE="${DR2_PREV_RECIPE:-experiments/524-october-recipe-dr/work/recipe524.json}"   # override for tests only
PREV_VCC=submissions/524dr-october-recipe-abc.vcc
DRY="${DRY:-0}"; FLOOR="${DR2_FLOOR:-202610070000}"; CAP="${DR2_CAP:-202610081200}"
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"          # a DRY run never writes the real run's log, recipe or env.json
LOG=$W/${PFX}recipe_dr524dr2.log
mkdir -p "$W"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"' EXIT
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
builders(){
  ps -axo pid,comm,args | awk -v me=$$ '$1 != me && $2 ~ /^bash$/ { for (i = 3; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524\.sh|build_allsw658\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
echo $$ > "$W/${PFX}recipe_dr524dr2.pid"
log "start (pid $$, DRY=$DRY): 524dr2 = 524's procedure on the applied tree; floor $FLOOR cap $CAP"
[ -s "$PREV_RECIPE" ] || finish "NOT READY 524dr's recipe ($PREV_RECIPE) is missing: nothing to branch against"
if [ "$DRY" = 1 ]; then WANT_B="${DR2_BUILDER_SHA:-$(shasum -a 256 $BU | cut -d' ' -f1)}"; WANT_EM="${DR2_EMITTER_SHA:-$(shasum -a 256 $EMP | cut -d' ' -f1)}"
else
  [ -n "${DR2_BUILDER_SHA:-}" ] && [ -n "${DR2_EMITTER_SHA:-}" ] || finish "NOT READY DR2_BUILDER_SHA / DR2_EMITTER_SHA not given (the post-10-06 shas)"
  WANT_B="$DR2_BUILDER_SHA"; WANT_EM="$DR2_EMITTER_SHA"
  log "waiting for >= $FLOOR, no other build, def_day/helpers/builder/emitter == HEAD (builder/emitter == the given shas), free+inactive >= 25 GiB, disk >= 35 GiB"
  until [ "$(date -u +%Y%m%d%H%M)" -ge "$FLOOR" ] && ! builders && eq_head $DD && eq_head $HP && eq_head $BU && eq_head $EMP \
        && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$WANT_B" ] && [ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$WANT_EM" ] \
        && [ "$(fi_gib)" -ge 25 ] && [ "$(free_gib)" -ge 35 ]; do
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && finish "NOT READY prerequisites not met by $CAP UTC"
    sleep 60; done
  log "released: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16), builder $(shasum -a 256 $BU | cut -c1-16), emitter $(shasum -a 256 $EMP | cut -c1-16)"
fi

# ---- resolve the recipe NOW (after 10-05's mde_sign and 10-06's step-9 decisions)
RW="$W/run"; [ "$DRY" = 1 ] && RW="$W/run-dry"
cp docs/october-decisions.json "$W/${PFX}october-decisions.at-resolve.json"
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$RW" --name "$NAME" --allow-pending --json \
  > "$W/${PFX}recipe524.json" 2>> "$LOG" || finish "NOT READY the recipe resolver failed (rc 3 = an installed script lacks a knob or a dependency is violated)"
log "recipe: $(tr -d '\n ' < "$W/${PFX}recipe524.json" | cut -c1-700)"
BR="$(uv run python $X/branch658.py "$W/${PFX}recipe524.json" "$PREV_RECIPE")"
case "$BR" in STOP*) finish "NOT READY ${BR#STOP }";; N0*|N1*|N2*|N3*) ;; *) finish "NOT READY could not classify the recipe ($BR)";; esac
BRANCH="${BR%% *}"; NEW="${BR#* }"
log "branch $BRANCH (new vs 524dr: $NEW)"
H1KEY="$(printf '%s' "$NEW" | sed -nE 's/.*h1_dose\(([A-Z0-9_]+)\).*/\1/p')"
if [ -n "$H1KEY" ] && ! /usr/bin/grep -qE "\"$H1KEY\"[[:space:]]*:" experiments/524-october-recipe-dr/gate524.py; then
  finish "NOT READY $H1KEY (an H1-dose switch) is in the recipe but gate524.py has no SWITCH_KEYS row for it (658b 1: add it first)"
fi
ENVS="$(uv run python -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' "$W/${PFX}recipe524.json")"
ARGS="$(uv run python -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' "$W/${PFX}recipe524.json")"
case "$ARGS" in "$BUNDLE $RW $NAME "*) ;; *) finish "NOT READY resolved args do not start with the bundle/work/name asked for: $ARGS";; esac
uv run python - "$W/${PFX}env.json" "$W/${PFX}recipe524.json" "$BRANCH" "$NEW" "$WANT_B" "$WANT_EM" <<'PY'
import json, sys
o, rp, br, new, b, em = sys.argv[1:]
r = json.load(open(rp))
json.dump({"exports": r["env"], "def_day_args": r["args"], "pending_treated_off": r.get("pending_treated_off"), "branch": br,
           "new_vs_524dr": new, "builder_sha": b, "emitter_sha": em,
           "operator": "experiments/453-abc-dress-rehearsal/work/run/whitener_union_pruned_def.npz (66610be7)"}, open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json"
if [ "$DRY" = 1 ]; then
  log "DRY plan: builders now: $(builders && echo busy || echo none); disk $(free_gib) GiB; free+inactive $(fi_gib) GiB; 524dr package $( [ -s "$PREV_VCC" ] && echo present || echo MISSING)"
  log "DRY def_day: $ENVS bash $DD $ARGS"
  rm -rf "$W/run-dry"
  finish "NOT READY DRY=1 (plan only; nothing built; branch $BRANCH)"
fi
case "$BRANCH" in N0|N1) [ -s "$PREV_VCC" ] || finish "NOT READY branch $BRANCH needs 524dr's package for its check ($PREV_VCC missing)";; esac

# ---- the build: exactly recipe_dr524.sh's def_day call
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || finish "NOT READY 453's day operator is not 66610be7"
[ -e "$W/run" ] && mv "$W/run" "$W/run.aborted-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "def_day: $ENVS bash $DD $ARGS"
eval "env $ENVS bash $DD $ARGS" > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || finish "NOT READY def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
log "def_day: $(/usr/bin/grep -E "LODO key ON|per-target sharpening|same-line mode|programme|CONVERGED at|NOT converged|final:|postbuild PASS" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-900)"
eq_head $BU || finish "NOT READY builder differs from HEAD after def_day (its exit trap should have restored it)"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/${NAME}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || finish "NOT READY relabel failed"
[ "$(fi_gib)" -ge 23 ] || finish "NOT READY free+inactive $(fi_gib) GiB < 23 -- not prepping; rerun the prep by hand"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || finish "NOT READY prep dry run failed"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || finish "NOT READY prep failed"
PKGSHA="$(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
log "package $STEM.vcc sha256 $PKGSHA"

if [ "$BRANCH" = N0 ]; then
  # same recipe on the applied tree: the content must be 524dr's. Compare both packages through the same restore path.
  mkdir -p "$W/n0"; mv "$STEM.vcc" "$W/n0/"           # out of submissions/: identical content must never take a slot
  [ "$(free_gib)" -ge 30 ] || finish "NOT READY N0: disk $(free_gib) GiB < 30 for two restored copies; run the check by hand"
  rm -f "$STEM.h5ad" "$W/run/"*.h5ad
  A="$W/n0/r524dr.h5ad"; Bh="$W/n0/r524dr2.h5ad"; rm -f "$A" "$Bh"
  uv run python scripts/restore_prediction_h5ad.py "$PREV_VCC" "$A" >> "$LOG" 2>&1 || { rm -f "$A"; finish "NOT READY N0: could not restore 524dr"; }
  uv run python scripts/restore_prediction_h5ad.py "$W/n0/$(basename "$STEM").vcc" "$Bh" >> "$LOG" 2>&1 || { rm -f "$A" "$Bh"; finish "NOT READY N0: could not restore 524dr2"; }
  uv run python $X/bulkid658.py "$Bh" "$A" --exact --out "$W/n0/bulkid_n0.json" >> "$LOG" 2>&1; BID=$?
  rm -f "$A" "$Bh"
  [ "$BID" = 0 ] && finish "NOT READY N0 tree check IDENTICAL: no new switch, the applied tree carries the October path bit for bit (no slot)"
  finish "NOT READY N0 tree check DIFFERENT (bulkid rc $BID, $W/n0/bulkid_n0.json): the applied tree changed the October path -- diagnose"
fi

uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
if [ "$BRANCH" = N1 ]; then
  # only DE switches are new: the bulk must be 524dr's exactly (the primary read's falsifier, checked before the slot)
  [ "$(free_gib)" -ge 20 ] || finish "NOT READY N1: disk $(free_gib) GiB < 20 for 524dr's restored copy; run the check by hand"
  A="$W/r524dr.h5ad"; rm -f "$A"
  uv run python scripts/restore_prediction_h5ad.py "$PREV_VCC" "$A" >> "$LOG" 2>&1 || { rm -f "$A"; finish "NOT READY N1: could not restore 524dr"; }
  uv run python $X/bulkid658.py "$STEM.h5ad" "$A" --out "$W/bulkid_n1.json" >> "$LOG" 2>&1; BID=$?
  rm -f "$A"
  [ "$BID" = 0 ] || finish "NOT READY N1: the bulk moved vs 524dr with only DE switches new (bulkid rc $BID): tree or determinism defect -- diagnose"
  log "N1 bulk identity vs 524dr: IDENTICAL"
fi
uv run python experiments/524-october-recipe-dr/gate524.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate524 rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate524 failed (see the gate lines above)"
rm -f "$STEM.h5ad" "$W/run/"*.h5ad     # the package exists, its sha is logged and the gates passed
finish "READY 524dr2 $PKGSHA"
