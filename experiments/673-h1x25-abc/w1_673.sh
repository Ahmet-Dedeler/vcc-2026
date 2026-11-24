#!/usr/bin/env bash
# 673 W1 build (readiness [1efb7d] for the coordinator's 673 prereg 5c4f4bc5 + 1182d686; tag w1; for 10-21 s1).
# W1 = R-ABC2's resolved recipe (695: the frozen October recipe 04be2af3) + DEF_DAY_H1_STACK=x2.5 (690i: the dump subshell's
# VCC_H1_UNION_G = 0.7; x2, the decided default, is 0.56), no hook, on def_day on A/B/C (180's def-bundle) like 695.
# The export is added by THIS driver; docs/october-decisions.json (h1_stack DECIDED x2) is only read.
# Inert x2 KA: the frozen decisions resolve to NO H1 export (x2 = def_day's default), and R-KA2 (877e7c28) rebuilt R-ABC2's def_day
# output byte for byte from exactly that recipe; test690i I2/I5/I7 cover the shell paths. This driver re-asserts that the resolved
# recipe equals R-ABC2's exactly (no DEF_DAY_H1_STACK in it) before adding x2.5.
# Starts after C0 (699's c0811_699.sh) has exited (one builder at a time) and the build gates open; no start after CAP (10-21 12:00Z).
# gate673: content-aware diff vs R-ABC2 within {bulk_load, bulk_load_sha256, bulk_ctx_scale} containing bulk_load (the dump is
# renamed), the dump records VCC_H1_UNION_G 0.7, all of 699's rows. READY "w1" to READY-science.log; evidence build_w1/runw1.log.
# usage: nohup bash scripts/detach.sh bash experiments/673-h1x25-abc/w1_673.sh >> experiments/673-h1x25-abc/work/w1_673.nohup 2>&1 &
#        DRY=1 bash experiments/673-h1x25-abc/w1_673.sh   (plan only)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/673-h1x25-abc; W=$X/work
if [ -z "${W1_673_FROZEN:-}" ]; then
  mkdir -p "$W"; FZ="$W/.frozen-$$-w1_673.sh"; cp "$REPO/$X/w1_673.sh" "$FZ" || exit 1
  W1_673_FROZEN=1 exec bash "$FZ" "$@"
fi
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
DEC=docs/october-decisions.json
FROZEN="$DEC:04be2af3289dabd0 $DD:11ff102e652ae315 $BU:e99bbd4a63ad834c $EMP:0f56209d57e8385a $HP:5354e44900194ab2"
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/673w1-october-recipe-abc; NAME=w1
R453=experiments/453-abc-dress-rehearsal/work/run; OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
R2_RECIPE=experiments/695-rabc2/work/recipe695.json; R2_RUN=experiments/695-rabc2/work/run
H1=x2.5; H1G=0.7; DUMP25=dump_h1x25_cd4x1_s0.25.npz
DNM2SHA=b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87
RUNLOG_ARM=experiments/428-zero-gene-audit/work/build_w1/runw1.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
DRY="${DRY:-0}"; CAP="${W1_CAP:-202610211200}"
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"
LOG=$W/${PFX}w1_673.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; rm -f "$W/.frozen-$$-w1_673.sh"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"; rm -f "$W/.frozen-$$-w1_673.sh"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
frozen_ok(){
  local e f want; for e in $FROZEN; do f="${e%%:*}"; want="${e##*:}"
    [ "$(sha "$f" | cut -c1-16)" = "$want" ] || { echo "$f is $(sha "$f" | cut -c1-16), not the frozen $want"; return 1; }; done
}
builders(){
  ps -axo pid,comm,args | awk -v me=$$ '$1 != me && $2 ~ /^bash$/ { for (i = 3; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524\.sh|recipe_dr524dr2\.sh|build_allsw658\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
gate(){
  local a; for a in c0811_699.sh c1811_698.sh rka2_690.sh; do ps -axo comm=,args= | awk -v p="$a" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}' && { echo "$a alive (runs first)"; return 1; }; done
  [ -e "$W/HOLD" ] && { echo "HOLD file ($(head -c 200 "$W/HOLD"))"; return 1; }
  builders && { echo "another build is running"; return 1; }
  local f r; for f in $DD $HP $BU $EMP; do eq_head $f || { echo "$f != HEAD"; return 1; }; done
  r="$(frozen_ok)" || { echo "$r"; return 1; }
  /usr/bin/grep -q "$DNM2SHA" $BU && /usr/bin/grep -q "from vcc2026 import dn_mag2 as _dn_mag" $BU || { echo "builder does not pin dn_mag2 ${DNM2SHA:0:16}"; return 1; }
  local fg mp; fg="$(fi_gib)"; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$fg" -ge 20 ] || { echo "free+inactive $fg GiB < 20"; return 1; }
  [ "${mp:-0}" -ge 40 ] || { echo "memory_pressure free ${mp:-?} % < 40"; return 1; }
  [ "$(disk_gib)" -ge 50 ] || { echo "disk $(disk_gib) GiB < 50"; return 1; }
  echo "GO, no build, scripts == HEAD == frozen (builder $(sha $BU | cut -c1-16), def_day $(sha $DD | cut -c1-16)), free+inactive $fg GiB, pressure $mp %, disk $(disk_gib) GiB"
}
mkdir -p "$W"; echo $$ > "$W/${PFX}w1_673.pid"
log "start (pid $$, DRY=$DRY): W1 = R-ABC2 + DEF_DAY_H1_STACK=$H1 (dump VCC_H1_UNION_G $H1G), on def_day on A/B/C; after C0; cap $CAP"
if [ "$DRY" != 1 ]; then
  last=""
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
    capped && finish "NOT READY gates not met by cap $CAP ($r)"
    sleep 60; done
  log "gate PASS: $r"
fi
RW="$W/run"; [ "$DRY" = 1 ] && RW="$W/run-dry"
cp "$DEC" "$W/${PFX}october-decisions.at-resolve.json"
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$RW" --name "$NAME" --json \
  > "$W/${PFX}recipe673.json" 2>> "$LOG" || finish "NOT READY the recipe resolver failed"
CMP="$(uv run python - "$W/${PFX}recipe673.json" "$R2_RECIPE" "$BUNDLE" "$RW" "$NAME" <<'PY'
import json, sys
r, p2 = (json.load(open(x)) for x in sys.argv[1:3])
bundle, rw, name = sys.argv[3:6]
bad = []
if r.get("problems") or r.get("pending_treated_off"):
    bad.append(f"problems {r.get('problems')} pending {r.get('pending_treated_off')}")
if r["args"][:3] != [bundle, rw, name]:
    bad.append(f"args head {r['args'][:3]}")
if r["args"][3:] != p2["args"][3:]:
    bad.append(f"args {r['args'][3:]} vs R-ABC2 {p2['args'][3:]}")
if r["env"] != p2["env"]:
    bad.append(f"resolved env != R-ABC2's: {sorted(set(r['env'].items()) ^ set(p2['env'].items()))}")
if "DEF_DAY_H1_STACK" in r["env"]:
    bad.append("the frozen recipe already carries DEF_DAY_H1_STACK (x2 must resolve to no export)")
print("STOP " + "; ".join(bad) if bad else "OK")
PY
)"
[ "$CMP" = OK ] || finish "NOT READY recipe check: ${CMP#STOP }"
log "inert-x2 KA: the frozen decisions (h1_stack DECIDED x2) resolve to R-ABC2's recipe exactly, with no DEF_DAY_H1_STACK export (R-KA2 877e7c28 rebuilt R-ABC2 byte for byte from it)"
python3 - "$W/${PFX}recipe673.json" "$H1" <<'PY2' || finish "NOT READY could not add the H1 dose to the recipe"
import json, sys
p, h1 = sys.argv[1:3]; r = json.load(open(p))
r["env"]["DEF_DAY_H1_STACK"] = h1; r["added_by_673"] = {"DEF_DAY_H1_STACK": h1}
json.dump(r, open(p, "w"), indent=1)
PY2
log "recipe: R-ABC2's + DEF_DAY_H1_STACK=$H1 (the only addition)"
ENVS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' "$W/${PFX}recipe673.json")"
ARGS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' "$W/${PFX}recipe673.json")"
python3 - "$W/${PFX}env.json" "$W/${PFX}recipe673.json" "$(sha $BU)" "$(sha $EMP)" "$(sha $DD)" <<'PY'
import json, sys
o, rp, b, em, dd = sys.argv[1:]
r = json.load(open(rp))
json.dump({"exports": r["env"], "def_day_args": r["args"], "changed_vs_reference": ["DEF_DAY_H1_STACK"], "builder_sha": b,
           "emitter_sha": em, "def_day_sha": dd, "operator": "experiments/453-abc-dress-rehearsal/work/run/whitener_union_pruned_def.npz (66610be7)",
           "reference": "R-ABC2 (695; package 54e6c5cf; board 0.2739849, nnHQfcZHsT3as9EafA3i): experiments/695-rabc2/work/recipe695.json"}, open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json"
if [ "$DRY" = 1 ]; then
  log "DRY gates now: $(gate || true)"
  log "DRY def_day: $ENVS bash $DD $ARGS"
  rm -rf "$W/run-dry"
  finish "NOT READY DRY=1 (plan only; nothing built)"
fi
[ -e "$STEM.vcc" ] && finish "NOT READY a package already exists at $STEM.vcc"
[ "$(sha "$R453/whitener_union_pruned_def.npz")" = "$OPSHA" ] || finish "NOT READY 453's day operator is not 66610be7"
[ -e "$W/run" ] && mv "$W/run" "$W/run.aborted-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "def_day: $ENVS bash $DD $ARGS"
t0=$(date +%s)
eval "env $ENVS bash $DD $ARGS" > "$W/def_day.nohup" 2>&1
log "def_day rc $? in $(( ($(date +%s) - t0) / 60 )) min"
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || finish "NOT READY def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
log "def_day: $(/usr/bin/grep -E "LODO key ON|LODO strength rider|H1 stack|per-target sharpening|DE rider|CONVERGED at|NOT converged|final:|postbuild PASS|690h" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-1500)"
eq_head $BU || finish "NOT READY builder differs from HEAD after def_day (its exit trap should have restored it)"
r="$(frozen_ok)" || finish "NOT READY after def_day: $r"
G="$(uv run python -c 'import json,sys,numpy as np; print(json.loads(str(np.load(sys.argv[1], allow_pickle=False)["env"])).get("VCC_H1_UNION_G"))' "$W/run/$DUMP25" 2>/dev/null)"
[ "$G" = "$H1G" ] || finish "NOT READY the x2.5 dump ($W/run/$DUMP25) records VCC_H1_UNION_G '$G', not $H1G"
log "[w1] H1 dose check: $DUMP25 records VCC_H1_UNION_G $G"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/${NAME}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || finish "NOT READY relabel failed"
[ "$(fi_gib)" -ge 20 ] || finish "NOT READY free+inactive $(fi_gib) GiB < 20 -- not prepping; rerun the prep by hand"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || finish "NOT READY prep dry run failed"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || finish "NOT READY prep failed"
PKGSHA="$(sha "$STEM.vcc")"
log "package $STEM.vcc sha256 $PKGSHA"
rm -f "$W/run/"*.h5ad          # the package exists and its sha is logged
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate673.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; GR=$?
log "gate673 rc=$GR"
[ "$GR" = 0 ] || finish "NOT READY gate673 failed (see the gate lines above)"
rm -f "$STEM.h5ad"
mkdir -p "$(dirname "$RUNLOG_ARM")"; cp "$LOG" "$RUNLOG_ARM"
printf '%s READY w1 %s %s\n' "$(date -u +%FT%TZ)" "$PKGSHA" "$STEM.vcc" >> "$READYLOG"
log "evidence copied to $RUNLOG_ARM; READY line appended to $READYLOG"
finish "READY w1 $PKGSHA"
