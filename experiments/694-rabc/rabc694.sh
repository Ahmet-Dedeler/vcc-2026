#!/usr/bin/env bash
# 694 R-ABC (2026-10-15, session [17adab]; asked by [400b6e] for the 10-16 s1 slot, s2 fallback). PREREGISTRATION.md.
# The October recipe exactly as resolved (`october_recipe.py command`, NO --allow-pending: every switch is decided), run by
# def_day.sh on 180's def-bundle (A/B/C relabelled D/E/F; same-panel branch: 453's operator 66610be7 pre-placed), exactly as
# 524dr2 was built (658's recipe_dr524dr2.sh). Then relabel453 -> vcc prep -> novelty -> bulk identity vs the board-scored
# 524dr2 (bulkid658 sums: the riders keep every per-(context, target, gene) total, so pds / mse must be 524dr2's) -> gate694.
# Never submits: the coordinator arms arm1015 on this log.
#
# Waits, in order (AGENTS 6.2; each capped by RABC_CAP):
#   1. run690 (R-NEW + R-KA) has EXITED (process, then its log tail is recorded); ka693b has exited with KAU2 PASS and
#      DNV2 PASS (the build-level known answers of the v2 tree: an all-unset rebuild == ec3a5ae6, and DN's scored prediction)
#   2. no $W/HOLD file (the coordinator's pause), no other build, def_day / helpers / builder / emitter == HEAD, the builder
#      pins dn_mag2 b9d6cfcf, free+inactive >= 25 GiB, disk >= 50 GiB
# Stops NOT READY if the resolved recipe is not 524dr2's env plus exactly DEF_DAY_DE_UPRN=b_k100_q20r3p10 and
# DEF_DAY_DE_DNMAG=d_m200 with the same arguments.
#
# Log = ready log = run log: experiments/694-rabc/work/rabc694.log. gate694 prints its "[rabc] PASS ..." / "[rabc] CLIFF OK ..."
# lines RAW into it (arm1015 greps them anchored at line start); "novelty rc=<n>" and "package <path> sha256 <sha>" are
# logged. The LAST line is always ONE terminal line:
#   <ts> READY rabc <package sha256>
#   <ts> NOT READY <reason>
# usage: nohup bash experiments/694-rabc/rabc694.sh > experiments/694-rabc/work/rabc694.nohup 2>&1 &
#        DRY=1 bash experiments/694-rabc/rabc694.sh    (no waits: resolve + check the recipe, print the plan, build nothing)
# Named rabc694.sh (not chainNNN.sh) so no build gate's chain pattern mistakes it for a builder. bash 3.2 / BSD grep safe.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/694-rabc; W=$X/work
# Run from a frozen per-launch copy, so an edit of this file never corrupts a live run (bash reads scripts by offset)
if [ -z "${RABC694_FROZEN:-}" ]; then
  mkdir -p "$W"; FZ="$W/.frozen-$$-rabc694.sh"; cp "$REPO/$X/rabc694.sh" "$FZ" || exit 1
  RABC694_FROZEN=1 exec bash "$FZ" "$@"
fi
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/694rabc-october-recipe-abc; NAME=rabc
R453=experiments/453-abc-dress-rehearsal/work/run; OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
PREV_RECIPE=experiments/658-prereg-1007-builds/work/dr2/recipe524.json      # 524dr2's resolved recipe (10-06)
PREV_VCC=submissions/524dr2-october-recipe-abc.vcc; PREV_SHA=93e49b8186bbaeeeae6696df8a20e16ebb9f214d1698a494de0d422a4aae106f
DNM2SHA=b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87
L690=experiments/690-def-readiness-1015/work/run690.log; L693B=experiments/693-dnmag-speed/work/ka693b.log
DRY="${DRY:-0}"; CAP="${RABC_CAP:-202610161700}"
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"          # a DRY run never writes the real run's log, recipe or env.json
LOG=$W/${PFX}rabc694.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
live(){ ps -axo comm=,args= | awk -v pat="$1" '$1 ~ /bash$/ && $0 ~ pat' | /usr/bin/grep -q .; }   # $1 = an awk regex on args
builders(){
  ps -axo pid,comm,args | awk -v me=$$ '$1 != me && $2 ~ /^bash$/ { for (i = 3; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524\.sh|recipe_dr524dr2\.sh|build_allsw658\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
gate(){
  [ -e "$W/HOLD" ] && { echo "HOLD file ($(head -c 200 "$W/HOLD"))"; return 1; }
  builders && { echo "another build is running"; return 1; }
  local f; for f in $DD $HP $BU $EMP; do eq_head $f || { echo "$f != HEAD"; return 1; }; done
  /usr/bin/grep -q "$DNM2SHA" $BU && /usr/bin/grep -q "from vcc2026 import dn_mag2 as _dn_mag" $BU || { echo "builder does not pin dn_mag2 ${DNM2SHA:0:16}"; return 1; }
  [ "$(fi_gib)" -ge 25 ] || { echo "free+inactive $(fi_gib) GiB < 25"; return 1; }
  [ "$(disk_gib)" -ge 50 ] || { echo "disk $(disk_gib) GiB < 50"; return 1; }
  echo "no build, scripts == HEAD (builder $(sha $BU | cut -c1-16), emitter $(sha $EMP | cut -c1-16)), free+inactive $(fi_gib) GiB, disk $(disk_gib) GiB"
}
mkdir -p "$W"; echo $$ > "$W/${PFX}rabc694.pid"
log "start (pid $$, DRY=$DRY): R-ABC = the resolved October recipe on def_day on A/B/C; cap $CAP"

if [ "$DRY" != 1 ]; then
  # ---- 1. run690 and ka693b first (the coordinator's order: ka693b, then the two rehearsals, then R-ABC)
  P690='(\.frozen-[0-9]+-run690\.sh|690-def-readiness-1015/run690\.sh)'; P693B='693-dnmag-speed/ka693b\.sh'
  if live "$P690" || live "$P693B"; then
    log "waiting for run690 (R-NEW + R-KA) and ka693b to exit"
    while live "$P690" || live "$P693B"; do capped && finish "NOT READY run690 / ka693b still running at cap $CAP"; sleep 60; done
  fi
  log "run690 exited; its tail: $(tail -3 "$L690" 2>/dev/null | tr '\n' ' ' | cut -c1-500)"
  KAV="$(/usr/bin/grep -E 'KAU2 (PASS|FAIL)|DNV2 (PASS|FAIL)' "$L693B" 2>/dev/null | tr '\n' ' ' | cut -c1-400)"
  log "ka693b verdicts: ${KAV:-none}"
  /usr/bin/grep -q "KAU2 PASS" "$L693B" && /usr/bin/grep -q "DNV2 PASS" "$L693B" \
    || finish "NOT READY ka693b did not pass both known answers (KAU2 + DNV2): the v2 tree is not proven -- diagnose first"
  # ---- 2. the build gates
  last=""
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
    capped && finish "NOT READY gates not met by cap $CAP ($r)"
    sleep 60; done
  log "gate PASS: $r"
fi

# ---- resolve the recipe NOW (plain `command`: a PENDING switch would make it fail, which is a stop)
RW="$W/run"; [ "$DRY" = 1 ] && RW="$W/run-dry"
cp docs/october-decisions.json "$W/${PFX}october-decisions.at-resolve.json"
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$RW" --name "$NAME" --json \
  > "$W/${PFX}recipe694.json" 2>> "$LOG" || finish "NOT READY the recipe resolver failed (a PENDING switch, or a missing knob / dependency)"
log "recipe: $(tr -d '\n ' < "$W/${PFX}recipe694.json" | cut -c1-700)"
CMP="$(uv run python - "$W/${PFX}recipe694.json" "$PREV_RECIPE" "$BUNDLE" "$RW" "$NAME" <<'PY'
import json, sys
r, p = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
bundle, rw, name = sys.argv[3:6]
want_new = {"DEF_DAY_DE_UPRN": "b_k100_q20r3p10", "DEF_DAY_DE_DNMAG": "d_m200"}
bad = []
if r.get("problems") or r.get("pending_treated_off"):
    bad.append(f"problems {r.get('problems')} pending {r.get('pending_treated_off')}")
if r["args"][:3] != [bundle, rw, name]:
    bad.append(f"args head {r['args'][:3]}")
if r["args"][3:] != p["args"][3:]:
    bad.append(f"args {r['args'][3:]} vs 524dr2 {p['args'][3:]}")
for k, v in p["env"].items():
    if r["env"].get(k) != v:
        bad.append(f"{k}={r['env'].get(k)!r} (524dr2 {v!r})")
new = {k: v for k, v in r["env"].items() if k not in p["env"]}
if new != want_new:
    bad.append(f"new exports {new} (want exactly {want_new})")
print("STOP " + "; ".join(bad) if bad else "OK " + json.dumps(new, sort_keys=True))
PY
)"
case "$CMP" in OK*) log "recipe check: 524dr2's env + ${CMP#OK }; args ${NAME} $(python3 -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["args"][3:]))' "$W/${PFX}recipe694.json")";;
  *) finish "NOT READY recipe is not 524dr2's plus exactly UPRN + DNMAG: ${CMP#STOP }";; esac
ENVS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' "$W/${PFX}recipe694.json")"
ARGS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' "$W/${PFX}recipe694.json")"
python3 - "$W/${PFX}env.json" "$W/${PFX}recipe694.json" "$(sha $BU)" "$(sha $EMP)" "$(sha $DD)" <<'PY'
import json, sys
o, rp, b, em, dd = sys.argv[1:]
r = json.load(open(rp))
json.dump({"exports": r["env"], "def_day_args": r["args"], "builder_sha": b, "emitter_sha": em, "def_day_sha": dd,
           "operator": "experiments/453-abc-dress-rehearsal/work/run/whitener_union_pruned_def.npz (66610be7)",
           "reference": "524dr2 (entry dVHLws6KdMxudq0Yn5ap, 0.2730889): experiments/658-prereg-1007-builds/work/dr2/recipe524.json"},
          open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json"
if [ "$DRY" = 1 ]; then
  log "DRY gates now: $(gate || true)"
  log "DRY 524dr2 package $( [ -s "$PREV_VCC" ] && echo present || echo MISSING); ka693b verdicts now: $(/usr/bin/grep -E 'KAU2 (PASS|FAIL)|DNV2 (PASS|FAIL)' "$L693B" 2>/dev/null | tr '\n' ' ')"
  log "DRY def_day: $ENVS bash $DD $ARGS"
  rm -rf "$W/run-dry"
  finish "NOT READY DRY=1 (plan only; nothing built)"
fi
[ -s "$PREV_VCC" ] && [ "$(sha "$PREV_VCC")" = "$PREV_SHA" ] || finish "NOT READY 524dr2's package is missing or not 93e49b81 (the bulk judge)"
[ -e "$STEM.vcc" ] && finish "NOT READY a package already exists at $STEM.vcc"

# ---- the build: exactly recipe_dr524dr2.sh's def_day call
[ "$(sha "$R453/whitener_union_pruned_def.npz")" = "$OPSHA" ] || finish "NOT READY 453's day operator is not 66610be7"
[ -e "$W/run" ] && mv "$W/run" "$W/run.aborted-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "def_day: $ENVS bash $DD $ARGS"
t0=$(date +%s)
eval "env $ENVS bash $DD $ARGS" > "$W/def_day.nohup" 2>&1
log "def_day rc $? in $(( ($(date +%s) - t0) / 60 )) min"
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || finish "NOT READY def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
log "def_day: $(/usr/bin/grep -E "LODO key ON|per-target sharpening|same-line mode|DE rider|CONVERGED at|NOT converged|final:|postbuild PASS|690h" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-1200)"
eq_head $BU || finish "NOT READY builder differs from HEAD after def_day (its exit trap should have restored it)"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/${NAME}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || finish "NOT READY relabel failed"
[ "$(fi_gib)" -ge 23 ] || finish "NOT READY free+inactive $(fi_gib) GiB < 23 -- not prepping; rerun the prep by hand"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || finish "NOT READY prep dry run failed"
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || finish "NOT READY prep failed"
PKGSHA="$(sha "$STEM.vcc")"
log "package $STEM.vcc sha256 $PKGSHA"
rm -f "$W/run/"*.h5ad          # the package exists and its sha is logged; $STEM.h5ad stays for the bulk check
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"

# ---- F1 before the slot: the bulk must be 524dr2's exactly (per-(context, target, gene) sums)
[ "$(disk_gib)" -ge 20 ] || finish "NOT READY disk $(disk_gib) GiB < 20 for 524dr2's restored copy; run the bulk check by hand"
A="$W/r524dr2.h5ad"; rm -f "$A"
uv run python scripts/restore_prediction_h5ad.py "$PREV_VCC" "$A" >> "$LOG" 2>&1 || { rm -f "$A"; finish "NOT READY could not restore 524dr2"; }
uv run python experiments/658-prereg-1007-builds/bulkid658.py "$STEM.h5ad" "$A" --out "$W/bulkid_rabc.json" >> "$LOG" 2>&1; BID=$?
rm -f "$A"
[ "$BID" = 0 ] || finish "NOT READY F1: the bulk differs from 524dr2's (bulkid rc $BID, $W/bulkid_rabc.json): the riders or the tree moved a total -- diagnose"
log "F1 bulk identity vs 524dr2: IDENTICAL (per-(context, target, gene) sums)"
uv run python $X/gate694.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate694 rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate694 failed (see the gate lines above)"
rm -f "$STEM.h5ad"             # the package exists, its sha is logged and every gate passed
finish "READY rabc $PKGSHA"
