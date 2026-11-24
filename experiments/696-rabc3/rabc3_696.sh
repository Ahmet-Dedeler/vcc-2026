#!/usr/bin/env bash
# 696 R-ABC3 (readiness [4592ee]; 10-19 s2 Branch A, tag rabc3). PREREGISTRATION.md. R-ABC2's recipe + DEF_DAY_RC_CALLSET_CPM=3000
# (K 3,400, only if read711 says ADOPT K 3,400), built exactly as R-ABC2 (rabc2_695.sh, from which this is derived). Starts only after
# GO-rabc3 (readiness touches it on the coordinator's ADOPT message) and the queue / memory gates; R-ABC2 must be READY (its
# package is the bulk judge). The recipe must be R-ABC2's env + exactly DEF_DAY_RC_CALLSET_CPM=3000, same args; bulk identity vs R-ABC2
# (K moves only the DE call set) before READY; gate696. READY "rabc3" to READY-science.log; evidence to build_rabc3/runrabc3.log.
# usage: nohup bash scripts/detach.sh bash experiments/696-rabc3/rabc3_696.sh >> experiments/696-rabc3/work/rabc3_696.nohup 2>&1 &
#        DRY=1 bash experiments/696-rabc3/rabc3_696.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/696-rabc3; W=$X/work
if [ -z "${RABC3_696_FROZEN:-}" ]; then
  mkdir -p "$W"; FZ="$W/.frozen-$$-rabc3_696.sh"; cp "$REPO/$X/rabc3_696.sh" "$FZ" || exit 1
  RABC3_696_FROZEN=1 exec bash "$FZ" "$@"
fi
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/696rabc3-october-recipe-abc; NAME=rabc3
R453=experiments/453-abc-dress-rehearsal/work/run; OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
PREV_RECIPE=experiments/695-rabc2/work/recipe695.json                   # R-ABC2's resolved recipe (the reference)
PREV_VCC=submissions/695rabc2-october-recipe-abc.vcc
# R-ABC2's package sha, read from its own READY line at run time (it does not exist when this driver is written)
PREV_SHA="$(/usr/bin/grep -oE "READY rabc2 [0-9a-f]{64}" experiments/695-rabc2/work/rabc2_695.log 2>/dev/null | tail -1 | awk '{print $3}')"
DNM2SHA=b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87
RUNLOG_ARM=experiments/428-zero-gene-audit/work/build_rabc3/runrabc3.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
DRY="${DRY:-0}"; CAP="${RABC3_CAP:-202610191500}"
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"
LOG=$W/${PFX}rabc3_696.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; rm -f "$W/.frozen-$$-rabc3_696.sh"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"; rm -f "$W/.frozen-$$-rabc3_696.sh"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
builders(){
  ps -axo pid,comm,args | awk -v me=$$ '$1 != me && $2 ~ /^bash$/ { for (i = 3; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524\.sh|recipe_dr524dr2\.sh|build_allsw658\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
gate(){
  [ -e "$W/GO-rabc3" ] || { echo "no GO-rabc3 (readiness touches it after the last verdict)"; return 1; }
  [ -e "$W/HOLD" ] && { echo "HOLD file ($(head -c 200 "$W/HOLD"))"; return 1; }
  builders && { echo "another build is running"; return 1; }
  local f; for f in $DD $HP $BU $EMP; do eq_head $f || { echo "$f != HEAD"; return 1; }; done
  /usr/bin/grep -q "$DNM2SHA" $BU && /usr/bin/grep -q "from vcc2026 import dn_mag2 as _dn_mag" $BU || { echo "builder does not pin dn_mag2 ${DNM2SHA:0:16}"; return 1; }
  local fg mp; fg="$(fi_gib)"; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$fg" -ge 20 ] || { echo "free+inactive $fg GiB < 20"; return 1; }
  [ "${mp:-0}" -ge 40 ] || { echo "memory_pressure free ${mp:-?} % < 40"; return 1; }
  [ "$(disk_gib)" -ge 50 ] || { echo "disk $(disk_gib) GiB < 50"; return 1; }
  echo "GO, no build, scripts == HEAD (builder $(sha $BU | cut -c1-16), def_day $(sha $DD | cut -c1-16)), free+inactive $fg GiB, pressure $mp %, disk $(disk_gib) GiB"
}
mkdir -p "$W"; echo $$ > "$W/${PFX}rabc3_696.pid"
log "start (pid $$, DRY=$DRY): R-ABC3 = R-ABC2 + K 3,400 (Branch A of 10-19 s2), on def_day on A/B/C; cap $CAP"
if [ "$DRY" != 1 ]; then
  last=""
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
    capped && finish "NOT READY gates not met by cap $CAP ($r)"
    sleep 60; done
  log "gate PASS: $r"
fi

RW="$W/run"; [ "$DRY" = 1 ] && RW="$W/run-dry"
cp docs/october-decisions.json "$W/${PFX}october-decisions.at-resolve.json"
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$RW" --name "$NAME" --json \
  > "$W/${PFX}recipe696.json" 2>> "$LOG" || finish "NOT READY the recipe resolver failed (a PENDING switch, or a missing knob / dependency)"
log "recipe: $(tr -d '\n ' < "$W/${PFX}recipe696.json" | cut -c1-900)"
CMP="$(uv run python - "$W/${PFX}recipe696.json" "$PREV_RECIPE" "$BUNDLE" "$RW" "$NAME" <<'PY'
import json, sys
r, p = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
bundle, rw, name = sys.argv[3:6]
allowed = {"DEF_DAY_RC_CALLSET_CPM": {"3000"}}   # 696: R-ABC3 = R-ABC2 + K 3,400 only
bad = []
if r.get("problems") or r.get("pending_treated_off"):
    bad.append(f"problems {r.get('problems')} pending {r.get('pending_treated_off')}")
if r["args"][:3] != [bundle, rw, name]:
    bad.append(f"args head {r['args'][:3]}")
if r["args"][3:] != p["args"][3:]:
    bad.append(f"args {r['args'][3:]} vs R-ABC2 {p['args'][3:]}")
for k, v in p["env"].items():
    if r["env"].get(k) != v:
        bad.append(f"{k}={r['env'].get(k)!r} (R-ABC2 {v!r})")
new = {k: v for k, v in r["env"].items() if k not in p["env"]}
for k, v in new.items():
    if v not in allowed.get(k, set()):
        bad.append(f"unexpected export {k}={v!r} (allowed {sorted(allowed)})")
if not bad and not new:
    print("NOTHING K 3,400 is not in the recipe: R-ABC3 is not built (Branch A needs read711 ADOPT K 3,400)")
else:
    print("STOP " + "; ".join(bad) if bad else "OK " + json.dumps(new, sort_keys=True))
PY
)"
case "$CMP" in
  OK*) log "recipe check: R-ABC2's env + ${CMP#OK }; args $NAME $(python3 -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["args"][3:]))' "$W/${PFX}recipe696.json")";;
  NOTHING*) finish "NOT READY nothing new: ${CMP#NOTHING }";;
  *) finish "NOT READY recipe is not R-ABC's plus adopted switches only: ${CMP#STOP }";;
esac
ADOPTED="$(python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); p=json.load(open(sys.argv[2])); print(",".join(sorted(k for k in r["env"] if k not in p["env"])))' "$W/${PFX}recipe696.json" "$PREV_RECIPE")"
ENVS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' "$W/${PFX}recipe696.json")"
ARGS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' "$W/${PFX}recipe696.json")"
python3 - "$W/${PFX}env.json" "$W/${PFX}recipe696.json" "$(sha $BU)" "$(sha $EMP)" "$(sha $DD)" "$ADOPTED" <<'PY'
import json, sys
o, rp, b, em, dd, ad = sys.argv[1:]
r = json.load(open(rp))
json.dump({"exports": r["env"], "def_day_args": r["args"], "adopted": ad.split(",") if ad else [], "builder_sha": b,
           "emitter_sha": em, "def_day_sha": dd, "operator": "experiments/453-abc-dress-rehearsal/work/run/whitener_union_pruned_def.npz (66610be7)",
           "reference": "R-ABC (entry PooEbFClgTLMQoXiOCQP, 0.2743318): experiments/694-rabc/work/recipe694.json"}, open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json (adopted: $ADOPTED)"
if [ "$DRY" = 1 ]; then
  log "DRY gates now: $(gate || true)"
  log "DRY def_day: $ENVS bash $DD $ARGS"
  rm -rf "$W/run-dry"
  finish "NOT READY DRY=1 (plan only; nothing built)"
fi
[ -n "$PREV_SHA" ] && [ -s "$PREV_VCC" ] && [ "$(sha "$PREV_VCC")" = "$PREV_SHA" ] || finish "NOT READY R-ABC2's package is missing or not its READY sha (the bulk judge)"
[ -e "$STEM.vcc" ] && finish "NOT READY a package already exists at $STEM.vcc"
[ "$(sha "$R453/whitener_union_pruned_def.npz")" = "$OPSHA" ] || finish "NOT READY 453's day operator is not 66610be7"
[ -e "$W/run" ] && mv "$W/run" "$W/run.aborted-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "def_day: $ENVS bash $DD $ARGS"
t0=$(date +%s)
eval "env $ENVS bash $DD $ARGS" > "$W/def_day.nohup" 2>&1
log "def_day rc $? in $(( ($(date +%s) - t0) / 60 )) min"
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || finish "NOT READY def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
log "def_day: $(/usr/bin/grep -E "LODO key ON|LODO strength rider|690g: dump|per-target sharpening|DE rider|CONVERGED at|NOT converged|final:|postbuild PASS|690h" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-1500)"
eq_head $BU || finish "NOT READY builder differs from HEAD after def_day (its exit trap should have restored it)"
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
if [ "$ADOPTED" = "DEF_DAY_RC_CALLSET_CPM" ]; then   # 696: always (K is the only allowed addition)
  # K alone moves only the DE call set: the bulk must be R-ABC's exactly (pds / mse identical on the board)
  [ "$(disk_gib)" -ge 20 ] || finish "NOT READY disk $(disk_gib) GiB < 20 for R-ABC's restored copy"
  A="$W/r695rabc2.h5ad"; rm -f "$A"
  uv run python scripts/restore_prediction_h5ad.py "$PREV_VCC" "$A" >> "$LOG" 2>&1 || { rm -f "$A"; finish "NOT READY could not restore R-ABC2"; }
  uv run python experiments/658-prereg-1007-builds/bulkid658.py "$STEM.h5ad" "$A" --out "$W/bulkid_rabc3.json" >> "$LOG" 2>&1; BID=$?
  rm -f "$A"
  [ "$BID" = 0 ] || finish "NOT READY K-only: the bulk differs from R-ABC2's (bulkid rc $BID) -- diagnose"
  log "K-only bulk identity vs R-ABC2: IDENTICAL"
fi
uv run python $X/gate696.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate696 rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate696 failed (see the gate lines above)"
rm -f "$STEM.h5ad"
mkdir -p "$(dirname "$RUNLOG_ARM")"; cp "$LOG" "$RUNLOG_ARM"
printf '%s READY rabc3 %s %s\n' "$(date -u +%FT%TZ)" "$PKGSHA" "$STEM.vcc" >> "$READYLOG"
log "evidence copied to $RUNLOG_ARM; READY line appended to $READYLOG"
finish "READY rabc3 $PKGSHA"
