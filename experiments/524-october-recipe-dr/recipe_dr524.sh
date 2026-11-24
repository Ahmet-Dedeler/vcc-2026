#!/usr/bin/env bash
# 524: the October recipe dress rehearsal (the 10-06 slot-1 candidate; PREREGISTRATION.md).
# At >= 10-05 02:00Z (after the 10-05 reads and the coordinator's `october_recipe.py decide` calls) it resolves
# `october_recipe.py command --allow-pending --json` on 180's def-bundle (A/B/C relabelled D/E/F), saves it as
# work/recipe524.json, clones 453's day operator (66610be7), runs def_day.sh with EXACTLY the resolved exports and
# arguments, relabels to A/B/C, preps (dry run, then real), checks novelty and runs gate524.py --write vs 453dr.
# Serial: no other build running; def_day.sh, its helpers, the builder and the emitter equal HEAD. Never submits.
# Named recipe_dr524.sh (not chain5NN.sh) so no build gate's `chain5..` pattern can mistake it for itself (AGENTS §6.2).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/524-october-recipe-dr; W=$X/work; LOG=$W/recipe_dr524.log; R453=experiments/453-abc-dress-rehearsal/work/run
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/524dr-october-recipe-abc
FLOOR=${RECIPE524_FLOOR:-202610050200}; CAP=${RECIPE524_CAP:-202610052000}
mkdir -p $W
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
builders(){
  # bash launchers, as an argv word ending in the script name (a mention inside one argument does not count),
  # plus anything scripts/def_pipeline.py other_builders() would call a builder (python/uv/time). 526.
  ps -axo comm,args | awk '$1 ~ /^bash$/ { for (i = 2; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
echo $$ > $W/recipe_dr524.pid
# 2026-10-02: floor 30 -> 25 GiB. def_day passes --min-free-gib 23 to def_pipeline, which counts free+inactive+speculative+
# purgeable (never less than this vm_stat free+inactive), so 25 keeps a 2 GiB margin above the build-stage refusal; at
# 30 the gate sat at 26-29 GiB with five Claude sessions open and 524dr could miss its 20:00Z cap (council [4d0207]).
log "armed (pid $$): waiting for >= $FLOOR UTC, no other build, def_day/helpers/builder/emitter == HEAD, free+inactive >= 25 GiB"
until [ "$(date -u +%Y%m%d%H%M)" -ge "$FLOOR" ] && ! builders && eq_head $DD && eq_head $HP && eq_head $BU && eq_head $EMP \
      && [ "$(fi_gib)" -ge 25 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "prerequisites not met by $CAP UTC -- stopping"; log "END"; exit 1; }
  sleep 60; done
log "released: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16), builder $(shasum -a 256 $BU | cut -c1-16), emitter $(shasum -a 256 $EMP | cut -c1-16)"
cp docs/october-decisions.json $W/october-decisions.at-resolve.json
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$W/run" --name 524dr --allow-pending --json \
  > $W/recipe524.json 2>> $LOG; RC=$?
[ $RC -eq 0 ] || { log "resolver rc=$RC (3 = an installed script lacks a knob or a dependency is violated) -- stopping"; log "END"; exit 1; }
log "recipe: $(tr -d '\n ' < $W/recipe524.json | cut -c1-700)"
# the exports and the arguments, exactly as resolved
ENVS="$(uv run python -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); assert not r["problems"], r["problems"]; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' $W/recipe524.json)" \
  || { log "recipe524.json has problems -- stopping"; log "END"; exit 1; }
ARGS="$(uv run python -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' $W/recipe524.json)"
case "$ARGS" in "$BUNDLE $W/run 524dr "*) ;; *) log "resolved args do not start with the bundle/work/name asked for: $ARGS -- stopping"; log "END"; exit 1;; esac
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || { log "453's day operator is not 66610be7 -- stopping"; log "END"; exit 1; }
[ -e "$W/run" ] && mv "$W/run" "$W/run.aborted-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "def_day: $ENVS bash $DD $ARGS"
eval "env $ENVS bash $DD $ARGS" > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-400)"; log "END"; exit 1; }
log "def_day: $(/usr/bin/grep -E "LODO key ON|per-target sharpening|same-line mode|programme|CONVERGED at|NOT converged|final:|postbuild PASS" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-900)"
eq_head $BU || { log "builder differs from HEAD after def_day (its exit trap should have restored it) -- stopping"; log "END"; exit 1; }
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/524dr_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; log "END"; exit 1; }
[ "$(fi_gib)" -ge 23 ] || { log "free+inactive $(fi_gib) GiB < 23 -- not prepping; rerun the prep by hand"; log "END"; exit 1; }
# prep overwrites the package before the gate. Keep the previous one so a failed re-run
# does not destroy a package whose gates524.log still matches (526).
[ -e "$STEM.vcc" ] && { cp -c "$STEM.vcc" "$STEM.vcc.prev" 2>/dev/null || cp "$STEM.vcc" "$STEM.vcc.prev"; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || { log "prep dry run failed"; log "END"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || { log "prep failed"; log "END"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate524.py --write >> "$LOG" 2>&1; G=$?; log "gate524 rc=$G"
if [ "$G" -ne 0 ]; then
  [ -e "$STEM.vcc.prev" ] && mv "$STEM.vcc.prev" "$STEM.vcc"
  log "END"
  exit 1
fi
rm -f "$STEM.vcc.prev" "$W/run/"*.h5ad
log "END recipe_dr524 (READY iff gate524 rc=0)"
