#!/usr/bin/env bash
# 698 C1 (readiness [4592ee]; tag c1811; for 10-20 s1; requested by the coordinator [2d6f46] 10-18, PLAN f0aaf5f9). PREREGISTRATION.md;
# the lab's prediction is experiments/811-lolo-model/PREREG-C1-ABC.md.
# C1 = R-ABC2H's exact recipe (697: the frozen October recipe + DEF_DAY_DUMP_HOOK=dumphook707.sh) with the hook path replaced by the
# lab's UNAPPLIED experiments/811-lolo-model/dumphook811.sh (swap707.py unchanged, then 811's c0_l25), on def_day on A/B/C like 697.
# docs/october-decisions.json is FROZEN (04be2af3) and only read: the hook is added by THIS driver, never by the decisions file.
# Starts only after readiness writes $W/GO-c1811 (after the lab's READY), with lines:
#   hook <path>        the committed, clean hook script
#   sha <sha256>       its sha (the lab's READY note)
#   ka_env <VAR=val>   the hook's lambda-0 switch: with it the hook must reproduce 697's dumphook707 output byte for byte on 697's
#                      own stack dump (2dfe327d -> af27efbd)   |   ka_skip <reason>   (only if the lab gives no switch; logged)
# Then the gates (no build, frozen shas, memory, disk). gate698: content-aware diff vs R-ABC2H = exactly bulk_load_sha256 (REAL) +
# bulk_ctx_scale, plus 690g's hook lines. READY "c1811" to READY-science.log; evidence to build_c1811/runc1811.log.
# No GO by CAP (default 10-19 19:00Z, so the build ends before 10-20 00:00Z) = NOT READY. Pause with $W/HOLD.
# usage: nohup bash scripts/detach.sh bash experiments/698-c1811/c1811_698.sh >> experiments/698-c1811/work/c1811_698.nohup 2>&1 &
#        DRY=1 bash experiments/698-c1811/c1811_698.sh     (plan only; needs a GO file, or DRY_HOOK=<path> to plan against)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/698-c1811; W=$X/work
if [ -z "${C1811_698_FROZEN:-}" ]; then
  mkdir -p "$W"; FZ="$W/.frozen-$$-c1811_698.sh"; cp "$REPO/$X/c1811_698.sh" "$FZ" || exit 1
  C1811_698_FROZEN=1 exec bash "$FZ" "$@"
fi
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
DEC=docs/october-decisions.json
# the frozen October recipe (coordinator, 10-18, after R-KA2 877e7c28): 16-hex prefixes
FROZEN="$DEC:04be2af3289dabd0 $DD:11ff102e652ae315 $BU:e99bbd4a63ad834c $EMP:0f56209d57e8385a $HP:5354e44900194ab2"
BUNDLE=experiments/180-def-readiness/work/def-bundle; STEM=submissions/698c1811-october-recipe-abc; NAME=c1811
R453=experiments/453-abc-dress-rehearsal/work/run; OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
R2H=experiments/697-rabc2h/work; R2H_RECIPE=$R2H/recipe697.json; R2_RECIPE=experiments/695-rabc2/work/recipe695.json
H707=experiments/707-strong-tier-direction/dumphook707.sh
KA_IN=$R2H/run/dump_h1x2_cd4x1_s0.25.npz; KA_IN_SHA=2dfe327de44733301f6f8ddea8a5dd3e52d5474eea10c132f31ef346ba91f3c8
KA_OUT_SHA=af27efbd7e8aaf18501b1400dba999b5c773f41bf9a13c69b3513ebc6c8d28a2   # 697's dumphook707 output on that dump
DNM2SHA=b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87
RUNLOG_ARM=experiments/428-zero-gene-audit/work/build_c1811/runc1811.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
DRY="${DRY:-0}"; CAP="${C1811_CAP:-202610191900}"; GO=$W/GO-c1811
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"
LOG=$W/${PFX}c1811_698.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; rm -f "$W/.frozen-$$-c1811_698.sh"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"; rm -f "$W/.frozen-$$-c1811_698.sh"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
eq_head(){ [ "$(sha "$1")" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
gofield(){ /usr/bin/grep -E "^$1 " "$GO" 2>/dev/null | tail -1 | cut -d' ' -f2-; }
frozen_ok(){
  local e f want; for e in $FROZEN; do f="${e%%:*}"; want="${e##*:}"
    [ "$(sha "$f" | cut -c1-16)" = "$want" ] || { echo "$f is $(sha "$f" | cut -c1-16), not the frozen $want"; return 1; }; done
}
builders(){
  ps -axo pid,comm,args | awk -v me=$$ '$1 != me && $2 ~ /^bash$/ { for (i = 3; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524\.sh|recipe_dr524dr2\.sh|build_allsw658\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
gate(){
  local a; for a in rabc2h_697.sh rka2_690.sh; do ps -axo comm=,args= | awk -v p="$a" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}' && { echo "$a alive (runs first)"; return 1; }; done
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
mkdir -p "$W"; echo $$ > "$W/${PFX}c1811_698.pid"
log "start (pid $$, DRY=$DRY): C1 = R-ABC2H with the hook path 707 -> 811, on def_day on A/B/C; GO file $GO; cap $CAP"

# ---- 1. GO (readiness writes it after the lab's READY)
if [ "$DRY" = 1 ] && [ ! -e "$GO" ]; then
  HOOK="${DRY_HOOK:-experiments/811-lolo-model/dumphook811.sh}"; HSHA="(dry: none)"; KAENV=""; KASKIP="dry run"
  log "DRY: no GO file; planning against $HOOK"
else
  last=""
  until [ -s "$GO" ]; do
    [ "$last" != w ] && log "waiting: no $GO (readiness writes it after the lab's READY)"; last=w
    capped && finish "NOT READY no GO by cap $CAP (the lab did not report READY in time)"
    sleep 60; done
  HOOK="$(gofield hook)"; HSHA="$(gofield sha)"; KAENV="$(gofield ka_env)"; KASKIP="$(gofield ka_skip)"
  log "GO: hook $HOOK sha $HSHA ka_env '${KAENV}' ka_skip '${KASKIP}'"
  [ -n "$HOOK" ] && [ -f "$HOOK" ] || finish "NOT READY GO names no hook file ($HOOK)"
  [ "$(sha "$HOOK")" = "$HSHA" ] || finish "NOT READY $HOOK is $(sha "$HOOK"), not the GO sha $HSHA"
  git ls-files --error-unmatch "$HOOK" >/dev/null 2>&1 || finish "NOT READY $HOOK is not committed (def_day requires a committed hook)"
  [ -z "$(git status --porcelain -- "$HOOK")" ] || finish "NOT READY $HOOK has uncommitted changes"
  [ -n "$KAENV" ] || [ -n "$KASKIP" ] || finish "NOT READY GO has neither ka_env nor ka_skip"
  case "$KAENV" in ""|[A-Za-z_]*=*) ;; *) finish "NOT READY ka_env must be VAR=value (got '$KAENV')";; esac
  log "hook committed and clean at $(git log -1 --format=%h -- "$HOOK"); $(git log -1 --format=%s -- "$HOOK" | cut -c1-160)"
fi
[ "$HOOK" != "$H707" ] || finish "NOT READY the GO hook is dumphook707 itself (that is R-ABC2H, not C1)"

# ---- 2. gates
if [ "$DRY" != 1 ]; then
  last=""
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
    capped && finish "NOT READY gates not met by cap $CAP ($r)"
    sleep 60; done
  log "gate PASS: $r"
fi

# ---- 3. the hook's lambda-0 known answer on 697's own stack dump (production's A/B/C theta dump)
if [ -n "$KAENV" ]; then
  [ "$(sha "$KA_IN")" = "$KA_IN_SHA" ] || finish "NOT READY 697's stack dump is not ${KA_IN_SHA:0:16}"
  KW="$W/${PFX}ka811"; rm -rf "$KW"; mkdir -p "$KW/in"; cp -c "$KA_IN" "$KW/in/"
  ( env "$KAENV" bash "$HOOK" "$BUNDLE" "$KW" "$KW/in/$(basename "$KA_IN")" "$R2H/run/identity.json" "$KW/hd" ) > "$KW/hook.log" 2>&1
  KRC=$?; KOUT="$KW/hd/$(basename "$KA_IN")"
  KSHA="$( [ -s "$KOUT" ] && sha "$KOUT" || echo none)"
  log "KA-off ($KAENV) on 697's stack dump ${KA_IN_SHA:0:16}: hook rc $KRC, out ${KSHA:0:16} vs dumphook707's ${KA_OUT_SHA:0:16}"
  [ "$KRC" = 0 ] && [ "$KSHA" = "$KA_OUT_SHA" ] || finish "NOT READY KA-off FAILED: with $KAENV the hook does not reproduce dumphook707 byte for byte (see $KW/hook.log)"
  log "[c1811] KA-off PASS: $KAENV reproduces dumphook707's output byte for byte on the A/B/C theta dump"
  rm -rf "$KW"
else
  log "KA-off NOT run: ${KASKIP:-no ka_env}"
fi
# 10-18 (coordinator): KA-on -- the hook as def_day will run it (no env) on the same dump must write the lab's C1 dump (GO "ka_on <sha>");
# after def_day, the dump the build used must be that file (checked in step 5)
KAON="$(gofield ka_on)"
if [ -n "$KAON" ]; then
  KW="$W/${PFX}ka811on"; rm -rf "$KW"; mkdir -p "$KW/in"; cp -c "$KA_IN" "$KW/in/"
  ( bash "$HOOK" "$BUNDLE" "$KW" "$KW/in/$(basename "$KA_IN")" "$R2H/run/identity.json" "$KW/hd" ) > "$KW/hook.log" 2>&1
  KRC=$?; KOUT="$KW/hd/$(basename "$KA_IN")"
  KSHA="$( [ -s "$KOUT" ] && sha "$KOUT" || echo none)"
  log "KA-on (as def_day runs it) on 697's stack dump ${KA_IN_SHA:0:16}: hook rc $KRC, out ${KSHA:0:16} vs the lab's ${KAON:0:16}"
  [ "$KRC" = 0 ] && [ "$KSHA" = "$KAON" ] || finish "NOT READY KA-on FAILED: the hook does not write the lab's C1 dump $KAON (see $KW/hook.log)"
  log "[c1811] KA-on PASS: the hook writes the lab's C1 dump ${KAON:0:16} byte for byte on the A/B/C theta dump"
  KAON_FILE="$KOUT"   # kept: step 5 compares the dump def_day used against it by content (each npz carries its run's own path)
fi

# ---- 4. the recipe: frozen decisions -> R-ABC2's recipe exactly; + the 811 hook; must then equal R-ABC2H's recipe except the hook path
RW="$W/run"; [ "$DRY" = 1 ] && RW="$W/run-dry"
cp "$DEC" "$W/${PFX}october-decisions.at-resolve.json"
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$RW" --name "$NAME" --json \
  > "$W/${PFX}recipe698.json" 2>> "$LOG" || finish "NOT READY the recipe resolver failed"
CMP="$(uv run python - "$W/${PFX}recipe698.json" "$R2_RECIPE" "$R2H_RECIPE" "$BUNDLE" "$RW" "$NAME" "$HOOK" <<'PY'
import json, sys
r, p2, p2h = (json.load(open(x)) for x in sys.argv[1:4])
bundle, rw, name, hook = sys.argv[4:8]
bad = []
if r.get("problems") or r.get("pending_treated_off"):
    bad.append(f"problems {r.get('problems')} pending {r.get('pending_treated_off')}")
if r["args"][:3] != [bundle, rw, name]:
    bad.append(f"args head {r['args'][:3]}")
if r["args"][3:] != p2["args"][3:]:
    bad.append(f"args {r['args'][3:]} vs R-ABC2 {p2['args'][3:]}")
if r["env"] != p2["env"]:
    bad.append(f"resolved env != R-ABC2's: {sorted(set(r['env'].items()) ^ set(p2['env'].items()))}")
want = dict(p2h["env"]); want["DEF_DAY_DUMP_HOOK"] = hook
got = dict(r["env"]); got["DEF_DAY_DUMP_HOOK"] = hook
if got != want:
    bad.append(f"recipe + 811 hook != R-ABC2H's recipe with the hook replaced: {sorted(set(got.items()) ^ set(want.items()))}")
if p2h["env"].get("DEF_DAY_DUMP_HOOK") != "experiments/707-strong-tier-direction/dumphook707.sh":
    bad.append(f"R-ABC2H's recipe hook is {p2h['env'].get('DEF_DAY_DUMP_HOOK')}")
print("STOP " + "; ".join(bad) if bad else "OK")
PY
)"
[ "$CMP" = OK ] || finish "NOT READY recipe check: ${CMP#STOP }"
python3 - "$W/${PFX}recipe698.json" "$HOOK" <<'PY2' || finish "NOT READY could not add the hook to the recipe"
import json, sys
p, hook = sys.argv[1:3]; r = json.load(open(p))
assert "DEF_DAY_DUMP_HOOK" not in r["env"], "the frozen recipe carries a hook"
r["env"]["DEF_DAY_DUMP_HOOK"] = hook; r["added_by_698"] = {"DEF_DAY_DUMP_HOOK": hook}
json.dump(r, open(p, "w"), indent=1)
PY2
log "recipe check: the frozen decisions resolve to R-ABC2's recipe exactly; + DEF_DAY_DUMP_HOOK=$HOOK == R-ABC2H's recipe with the hook path replaced"
ENVS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(r["env"].items())))' "$W/${PFX}recipe698.json")"
ARGS="$(python3 -c 'import json,shlex,sys; r=json.load(open(sys.argv[1])); print(" ".join(shlex.quote(a) for a in r["args"]))' "$W/${PFX}recipe698.json")"
python3 - "$W/${PFX}env.json" "$W/${PFX}recipe698.json" "$(sha $BU)" "$(sha $EMP)" "$(sha $DD)" "$HOOK" "$( [ -f "$HOOK" ] && sha "$HOOK" || echo none)" <<'PY'
import json, sys
o, rp, b, em, dd, hk, hs = sys.argv[1:]
r = json.load(open(rp))
json.dump({"exports": r["env"], "def_day_args": r["args"], "changed_vs_reference": ["DEF_DAY_DUMP_HOOK"], "hook": hk, "hook_sha": hs,
           "builder_sha": b, "emitter_sha": em, "def_day_sha": dd,
           "operator": "experiments/453-abc-dress-rehearsal/work/run/whitener_union_pruned_def.npz (66610be7)",
           "reference": "R-ABC2H (697; package a43a64d0; board read 10-19 s2): experiments/697-rabc2h/work/recipe697.json"}, open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json"
if [ "$DRY" = 1 ]; then
  log "DRY gates now: $(gate || true)"
  log "DRY def_day: $ENVS bash $DD $ARGS"
  rm -rf "$W/run-dry"
  finish "NOT READY DRY=1 (plan only; nothing built)"
fi

# ---- 5. build, package, gate
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
r="$(frozen_ok)" || finish "NOT READY after def_day: $r"
# a dump's file sha carries its run dir (the `env` string), so these compare CONTENT with def_metadiff.npz_same (arrays exact,
# strings equal once each run's own dir is normalised): 695's and 697's stack dumps differ in sha and are content-identical.
DB="$(basename "$KA_IN")"
npzeq(){ uv run python -c 'import sys; sys.path.insert(0, "scripts"); from pathlib import Path; from def_metadiff import npz_same
ok, why = npz_same(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])); print(why); sys.exit(0 if ok else 1)' "$@"; }
r="$(npzeq "$W/run/$DB" "$W/run" "$KA_IN" "$R2H/run")" || finish "NOT READY def_day's stack dump differs in content from 697's (the KAs ran on 697's): $r"
log "[c1811] stack dump == 697's in content ($r)"
if [ -n "$KAON" ]; then
  r="$(npzeq "$W/run/dumphook/$DB" "$W/run" "$KAON_FILE" "$R2H/run")" \
    || finish "NOT READY the dump def_day used differs in content from the lab's C1 dump ${KAON:0:16}: $r"
  log "[c1811] dump check PASS: the hook's dump def_day used == the lab's C1 dump ${KAON:0:16} in content ($r)"
  rm -rf "$W/${PFX}ka811on"
fi
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
uv run python $X/gate698.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate698 rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate698 failed (see the gate lines above)"
rm -f "$STEM.h5ad"
mkdir -p "$(dirname "$RUNLOG_ARM")"; cp "$LOG" "$RUNLOG_ARM"
printf '%s READY c1811 %s %s\n' "$(date -u +%FT%TZ)" "$PKGSHA" "$STEM.vcc" >> "$READYLOG"
log "evidence copied to $RUNLOG_ARM; READY line appended to $READYLOG"
finish "READY c1811 $PKGSHA"
