#!/usr/bin/env bash
# 691 — the Oct 22 release watcher (2026-10-15, session [17adab], D/E/F readiness).
#
# From START it polls `vcc --json datasets list` every EVERY seconds and compares it with a baseline taken when it was
# armed. A NEW or CHANGED available dataset is the release. Then, with nobody awake:
#   1. download it (`vcc datasets download <id> -d <DL>/<ts>`; the CLI resumes and verifies the server's checksum);
#   2. unpack it as untrusted data into a FRESH bundle dir (release691.py unpack: no traversal, no symlinks, data file
#      types only), recording every file's size and sha256;
#   3. run the driver's builder-free stages (preflight, identity, cd4-extract, census; no pin, ~15 s, < 4 GiB). If
#      preflight trips a guard, the other three still run so the identity read is on disk;
#   4. summarize (release691.py summarize): contexts, SAME-PANEL vs NEW-PANEL against A/B/C, gene axis identical,
#      depth and spread, identity calls; and write NEXT.txt with the exact day command for that branch;
#   5. only with AUTO_R0=1, and only if preflight and the summary are clean and the gates pass (no HOLD, builder == HEAD,
#      no builder, no day-path chain, free+inactive >= 26 GiB, disk >= 60 GiB): build R0 with def_day.sh and the resolved
#      recipe (plain `october_recipe.py command`, NO --allow-pending: a PENDING switch stops R0) in a fresh work dir.
#      It never runs `vcc prep` (real) or `vcc submit`: those stay by hand (DEF-RUNBOOK §2's tail, AGENTS §14).
# It touches data/challenge never (A/B/C stays where every A/B/C script expects it).
#
# Knobs (tests): VCC691 (the vcc binary), START691 / CAP691 (UTC YYYYmmddHHMM), EVERY691 (s), DL691, BUNDLE691, W691.
# Terminal line, always: "<ts> END". Signals for sessions: <W>/RELEASE-DETECTED, <W>/summary.json, <W>/NEXT.txt.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 1
E=experiments/691-release-watcher
W="${W691:-$E/work}"; LOG="$W/watch691.log"
# 10-16: START default 10-19 00:00Z (was 10-21 22:00Z) and CAP 10-27 00:00Z (was 10-23 18:00Z): Arc gave no day or hour ("next week")
START="${START691:-202610190000}"; CAP="${CAP691:-202610270000}"; EVERY="${EVERY691:-300}"
VCC="${VCC691:-vcc}"; DL="${DL691:-data/def-download}"; BUNDLE="${BUNDLE691:-data/challenge-def}"
AUTO_R0="${AUTO_R0:-0}"; R0W="${R0W691:-experiments/692-def-r0/work/r0}"
UPGRADE="${UPGRADE691:-uv tool upgrade vcc-cli}"; DLWAIT="${DLWAIT691:-60}"; DLHOURS="${DLHOURS691:-4}"; UPG=0
OPS=experiments/453-abc-dress-rehearsal/work/run; OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
mkdir -p "$W"
# Run from a frozen per-launch copy, so an edit of this file can never corrupt a live run (bash reads scripts by offset;
# AGENTS 6.2). The copy's name keeps "watch691.sh" for the single-instance guard and 691's relauncher.
if [ -z "${WATCH691_FROZEN:-}" ]; then
  for f in "$W"/.frozen-*-watch691.sh; do [ -e "$f" ] || continue; p="${f##*/.frozen-}"; p="${p%%-*}"; kill -0 "$p" 2>/dev/null || rm -f "$f"; done
  FZ="$W/.frozen-$$-watch691.sh"; cp "$REPO/experiments/691-release-watcher/watch691.sh" "$FZ" || exit 1
  WATCH691_FROZEN=1 exec bash "$FZ" "$@"
fi
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"' EXIT
now(){ date -u +%Y%m%d%H%M; }
rel(){ uv run python -I "$E/release691.py" "$@"; }
# macOS notifications only for the real CLI: a stub (rehearsals, tests) must never tell Ahmet that D/E/F has landed
notify(){ [ "$VCC" = vcc ] || return 0; osascript -e "display notification \"$1\" with title \"VCC D/E/F release\"" >/dev/null 2>&1 || true; }
# review 720 F4: one standing-approved upgrade of the vcc CLI (2025's release day needed one); the CLI is a uv tool, not the
# project venv, so no numerics KA is involved (AGENTS 11.0 lists numpy/scipy/anndata/h5py/cell-eval2)
upgrade_once(){
  [ "$UPG" = 1 ] && return 0; UPG=1
  log "UPGRADE (once, $1): vcc before: $("$VCC" --version 2>&1 | head -1)"
  $UPGRADE >> "$LOG" 2>&1; log "upgrade rc $?"; "$VCC" skill install >> "$LOG" 2>&1 || true
  log "UPGRADE: vcc after: $("$VCC" --version 2>&1 | head -1)"
  log "UPGRADE: the day's real vcc prep now runs on a new CLI: re-run DEF-RUNBOOK §6.7's prep KA before the first D/E/F submit (720 R4)"
  notify "vcc CLI upgraded on the day: re-run the prep KA"
}
# 720 R1: never two watchers on one work dir (both could pass the R0-dir checks and start def_day seconds apart)
if [ -s "$W/watch691.pid" ]; then
  _op="$(cat "$W/watch691.pid")"
  if [ "$_op" != "$$" ] && ps -p "$_op" -o args= 2>/dev/null | /usr/bin/grep -q "watch691.sh"; then
    log "STOP: another watch691 (pid $_op) is alive on $W; not starting a second one"; trap - EXIT; exit 1
  fi
fi
echo $$ > "$W/watch691.pid"
log "watch691 armed (pid $$): start $START, cap $CAP, every ${EVERY}s, vcc '$VCC', bundle -> $BUNDLE, AUTO_R0=$AUTO_R0, AUTO_ALLOW_DEPTH691=${AUTO_ALLOW_DEPTH691:-0}, AUTO_ALLOW_WIDE691=${AUTO_ALLOW_WIDE691:-0}"

# baseline: today's list (taken once; a re-arm keeps the first one)
if [ ! -s "$W/baseline.json" ]; then
  "$VCC" --json datasets list > "$W/baseline.json" 2>> "$LOG" || { log "STOP: baseline list failed"; exit 1; }
fi
log "baseline: $(python3 -c 'import json,sys; print([(d["id"], d.get("filename")) for d in json.load(open(sys.argv[1]))["datasets"]])' "$W/baseline.json")"
# RESUME (a re-launch after a reboot or a kill, AGENTS 6.2 step 6): a release already detected skips the poll; a bundle
# that matches its unpack record skips the download; the cheap stages always re-run.
RESUME=0; [ -s "$W/RELEASE-DETECTED" ] && RESUME=1
if [ "$RESUME" = 0 ]; then
[ -e "$BUNDLE" ] && { log "STOP: $BUNDLE already exists (move it aside; the watcher writes a fresh bundle only)"; exit 1; }

until [ "$(now)" -ge "$START" ]; do sleep 60; done
log "polling"
FAILS=0; IDS=""
while :; do
  [ "$(now)" -ge "$CAP" ] && { log "STOP: no release by CAP $CAP"; exit 3; }
  if "$VCC" --json datasets list > "$W/now.json" 2> "$W/now.err"; then
    [ "$FAILS" -gt 0 ] && log "list answers again after $FAILS failures"
    FAILS=0
    IDS="$(rel diff "$W/baseline.json" "$W/now.json")"
    [ -n "$IDS" ] && break
  else
    FAILS=$((FAILS + 1))
    # Cold Turkey or the network: log the first failure and then once an hour (12 x 5 min)
    [ $((FAILS % 12)) = 1 ] && log "list failed ($FAILS in a row): $(head -c 300 "$W/now.err" | tr '\n' ' ')"
    # a refusal that names the client version is not the network: upgrade once (720 F4)
    [ "$FAILS" -ge 3 ] && /usr/bin/grep -qiE "upgrade|unsupported|426|client version" "$W/now.err" && upgrade_once "list refused: $(head -c 120 "$W/now.err" | tr '\n' ' ')"
  fi
  sleep "$EVERY"
done
TS="$(date -u +%Y%m%dT%H%MZ)"
log "RELEASE DETECTED: $(echo $IDS | tr '\n' ' ')"; printf '%s %s\n' "$TS" "$(echo $IDS | tr '\n' ' ')" > "$W/RELEASE-DETECTED"
cp "$W/now.json" "$W/release-list.json"
notify "dataset(s) $(echo $IDS | tr '\n' ' ') appeared; downloading"
else
  TS="$(awk '{print $1}' "$W/RELEASE-DETECTED")"; IDS="$(cut -d' ' -f2- "$W/RELEASE-DETECTED" | tr ' ' '\n' | sed '/^$/d')"
  cp "$W/release-list.json" "$W/now.json" 2>/dev/null || { log "STOP: RESUME without release-list.json"; exit 1; }
  log "RESUME: release $(echo $IDS | tr '\n' ' ')detected at $TS; resuming"
fi
log "vcc $("$VCC" --version 2>&1 | head -1); PyPI vcc-cli $(curl -s -m 10 https://pypi.org/pypi/vcc-cli/json | python3 -c 'import json,sys; print(json.load(sys.stdin)["info"]["version"])' 2>/dev/null || echo '?')"

# 720 F5: choose by name among several new ids; ambiguous -> a person decides (NEXT.txt), nothing is downloaded
if [ -n "${PICK691:-}" ]; then   # a person's choice after an ambiguous release (re-launch with PICK691=<id>)
  echo "$IDS" | /usr/bin/grep -qx "$PICK691" && { ID="$PICK691"; PRC=0; log "PICK691: $ID chosen by hand"; } || { ID=""; PRC=2; log "PICK691=$PICK691 is not among the new ids"; }
else
  ID="$(rel pick "$W/now.json" $IDS 2>> "$LOG")"; PRC=$?
fi
if [ "$PRC" != 0 ] || [ -z "$ID" ]; then
  { echo "# 691 NEXT ($(date -u +%FT%TZ)): AMBIGUOUS release, nothing downloaded. New datasets:"; echo "$IDS" | sed 's/^/#   /'
    echo "# Pick the final one, then re-launch the watcher with PICK691=<id> in front of the same command (it resumes), or by hand:"
    echo "vcc datasets download <id> -d $DL/manual && uv run python -I $E/release691.py unpack <zip> $BUNDLE $W/release.json"; } > "$W/NEXT.txt"
  log "STOP: ambiguous release (pick rc $PRC); NEXT.txt lists the candidates"; notify "ambiguous D/E/F release: read NEXT.txt"; exit 1
fi
log "release dataset chosen: $ID (of $(echo $IDS | tr '\n' ' '))"
FN="$(python3 -c 'import json,sys; print([d for d in json.load(open(sys.argv[1]))["datasets"] if d["id"] == sys.argv[2]][0]["filename"])' "$W/now.json" "$ID")"
SKIPDL=0
if [ "$RESUME" = 1 ] && [ -e "$BUNDLE" ]; then
  if rel verify "$BUNDLE" "$W/release.json" >> "$LOG" 2>&1; then SKIPDL=1; log "RESUME: bundle verified; skipping download and unpack"
  else log "STOP: RESUME found $BUNDLE but it does not match release.json (move it aside, then re-launch)"; exit 1; fi
fi
[ "$RESUME" = 1 ] && [ -e "$BUNDLE.unzipped" ] && [ "$SKIPDL" = 0 ] && { mv "$BUNDLE.unzipped" "$BUNDLE.unzipped.aside-$(date -u +%H%M%S)"; log "RESUME: moved an interrupted unpack aside"; }
if [ "$SKIPDL" = 0 ]; then
mkdir -p "$DL/$TS"
# 720 F4: retry with backoff for DLHOURS (Cold Turkey's morning blocks lasted hours), one CLI upgrade after the 3rd failure
ok=0; i=0; d=$DLWAIT; T0DL=$(date +%s)
while :; do
  i=$((i + 1))
  if "$VCC" datasets download "$ID" -d "$DL/$TS" >> "$W/download.log" 2>&1; then ok=1; break; fi
  log "download attempt $i failed: $(tail -2 "$W/download.log" | tr '\n' ' ' | cut -c1-200)"; notify "download attempt $i failed"
  [ "$i" = 3 ] && upgrade_once "3 download failures"
  [ $(( $(date +%s) - T0DL )) -ge $(( DLHOURS * 3600 )) ] && break
  [ "$(now)" -ge "$CAP" ] && break
  sleep "$d"; d=$(( d * 2 )); [ "$d" -gt 900 ] && d=900
done
[ "$ok" = 1 ] || { log "STOP: download failed $i times over $(( ($(date +%s) - T0DL) / 60 )) min (see download.log)"; notify "download FAILED"; exit 1; }
log "download ok after $i attempt(s)"
ZIP="$DL/$TS/$FN"
[ -s "$ZIP" ] || { log "STOP: $ZIP missing after the download"; exit 1; }
log "downloaded $ZIP ($(stat -f %z "$ZIP") bytes, sha256 $(shasum -a 256 "$ZIP" | cut -c1-16))"
rel unpack "$ZIP" "$BUNDLE" "$W/release.json" >> "$LOG" 2>&1 || { log "STOP: unpack refused (see above)"; notify "unpack REFUSED"; exit 1; }
fi
log "df: $(df -g /System/Volumes/Data | awk 'NR==2{print $4}') GiB free"

# builder-free stages; a tripped preflight guard does not stop the identity read
LW="$W/light"; rm -rf "$LW"; mkdir -p "$LW"
DRV=(uv run python scripts/def_pipeline.py --bundle "$BUNDLE" --work "$LW" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$LW/frozen_inputs.json" --name light691)
"${DRV[@]}" --stages preflight > "$LW/preflight.out" 2>&1; PFRC=$?
log "preflight rc $PFRC: $(grep -E '^\s*FAIL' "$LW/logs/preflight.log" 2>/dev/null | tr -s ' ' | tr '\n' ';' | cut -c1-600)"
# 690 depth contingency: the DECISION-CARD §1 overrides that are pre-decided, only when armed (AUTO_ALLOW_DEPTH691 /
# AUTO_ALLOW_WIDE691) and only in their exact cases; preflight must then pass with them, and def_day gets the same env
OVR=""
if [ "$PFRC" != 0 ] && OVR="$(rel overrides "$LW/preflight.json" 2>> "$LOG")" && [ -n "$OVR" ]; then
  log "AUTO OVERRIDE (pre-decided, DECISION-CARD §1): $OVR; re-running preflight with it"; notify "D/E/F preflight override: $OVR"
  env $OVR "${DRV[@]}" --stages preflight > "$LW/preflight-ovr.out" 2>&1; PFRC=$?
  log "preflight with the override rc $PFRC"
  [ "$PFRC" = 0 ] || OVR=""
fi
"${DRV[@]}" --stages identity,cd4-extract,census > "$LW/stages.out" 2>&1; STRC=$?
log "identity/cd4-extract/census rc $STRC"
rel summarize "$BUNDLE" "$LW" "$W/summary.json" >> "$LOG" 2>&1; SMRC=$?
BRANCH="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["branch"])' "$W/summary.json" 2>/dev/null || echo UNKNOWN)"
log "summary rc $SMRC, branch $BRANCH"

# NEXT.txt: the exact day command for this branch (plain command: PENDING refuses, by design)
# R0_ALLOW_PENDING691=1 is for REHEARSALS ONLY (690's R-SAME): a PENDING switch then counts as off. Never on the day.
AP=(); [ "${R0_ALLOW_PENDING691:-0}" = 1 ] && { AP=(--allow-pending); log "REHEARSAL: R0_ALLOW_PENDING691=1 (PENDING treated off)"; }
CMD="$(uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$R0W" --name r0def ${AP[@]+"${AP[@]}"} 2>&1)"; CMDRC=$?
{
  echo "# 691 NEXT ($(date -u +%FT%TZ)); release $ID ($FN); branch $BRANCH; preflight rc $PFRC; light stages rc $STRC; summary rc $SMRC"
  echo "# Read $W/summary.json (identity reading rule, DEF-RUNBOOK top block) and the 690 checklist before building."
  [ "$BRANCH" = SAME-PANEL ] && echo "mkdir -p $R0W && cp -c $OPS/whitener_union_pruned_def.npz $OPS/whitener_union_pruned_def.json $R0W/   # SAME-PANEL branch (690): 453's operator 66610be7"
  [ "$PFRC" != 0 ] && echo "# preflight FAILED: read the guard table (DECISION-CARD §1) before any override"
  [ -n "$OVR" ] && echo "# AUTO OVERRIDE applied (pre-decided case): export $OVR before the def_day line below"
  echo "$CMD"
} > "$W/NEXT.txt"
log "NEXT.txt written (recipe command rc $CMDRC)"
notify "D/E/F unpacked: branch $BRANCH, preflight rc $PFRC; see experiments/691-release-watcher/work/NEXT.txt"

if [ "$AUTO_R0" != 1 ]; then log "AUTO_R0 off: stopping after NEXT.txt"; exit 0; fi
[ "$PFRC" = 0 ] && [ "$STRC" = 0 ] && [ "$SMRC" = 0 ] && [ "$CMDRC" = 0 ] || { log "AUTO_R0: not starting (preflight $PFRC, stages $STRC, summary $SMRC, recipe $CMDRC): a human reads NEXT.txt"; exit 0; }
if [ -e "$R0W" ]; then
  # RESUME: what the previous run left in the R0 dir decides (checked BEFORE the gates: a pinned builder is exactly
  # what an interrupted R0 leaves, and the gates would only wait on it)
  if /usr/bin/grep -q "builder after restore:" "$R0W/def_day.log" 2>/dev/null; then
    if /usr/bin/grep -q "STOP:" "$R0W/def_day.log"; then
      log "AUTO_R0: R0 def_day in $R0W STOPPED: $(/usr/bin/grep 'STOP:' "$R0W/def_day.log" | tail -1 | cut -c1-200). A person reruns it in a fresh dir (checklist §C2)"
      notify "R0 def_day STOPPED: read $R0W/def_day.log"; exit 1
    fi
    log "AUTO_R0: $R0W already holds a finished def_day ($(/usr/bin/grep -E 'DONE|STOP' "$R0W/def_day.log" | tail -1 | cut -c1-200)); 692 reads it; nothing to do"; exit 0
  fi
  if ps -axo comm=,args= | awk -v w="$R0W" '$1 ~ /bash$/ && /def_day\.sh/ && index($0, w) > 0' | /usr/bin/grep -q .; then
    log "AUTO_R0: a def_day for $R0W is still running; nothing to do"; exit 0
  fi
  _db="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  _ours="$(cat "$R0W"/pin-*.json 2>/dev/null | python3 -c 'import json,sys
t = sys.stdin.read()
dec = json.JSONDecoder(); i = 0; out = []
while i < len(t):
    t2 = t[i:].lstrip()
    if not t2: break
    o, n = dec.raw_decode(t2); i += len(t[i:]) - len(t2) + n
    out.append(o.get("builder_sha256_after", ""))
print(" ".join(x for x in out if x))' 2>/dev/null)"
  case " $_ours " in
    *" $_db "*)
      log "STOP: $R0W holds an INTERRUPTED def_day and the builder on disk is ITS pin (${_db:0:16}). Recovery (AGENTS 6.2):"
      log "      uv run python scripts/def_pin_panel_artifacts.py --restore --record $R0W/pin-artifacts.json"
      log "      then re-launch this watcher: it moves $R0W aside and rebuilds R0 (692 waits on the same path)"
      exit 1 ;;
  esac
  mv "$R0W" "$R0W.aborted-$(date -u +%H%M%S)"
  log "AUTO_R0: an interrupted def_day in $R0W (not holding the builder) was moved aside (AGENTS 6.2 step 4); rebuilding"
fi
# gates (as in 690's run690.sh), waiting up to CAP
gate(){
  [ -e "$W/HOLD" ] && { echo "HOLD file"; return 1; }
  local hb db; hb="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"; db="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  [ "$hb" = "$db" ] || { echo "builder != HEAD"; return 1; }
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99)" = 0 ] || { echo "a builder is running"; return 1; }
  [ "$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh)/' | wc -l | tr -d ' ')" = 0 ] || { echo "a day-path chain is running"; return 1; }
  local f; f="$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)"
  # 10-16: the AGENTS §9 rule (free+inactive >= 20 GiB AND memory_pressure free >= 40 %; no swap clause, no swap file since
  # 10-15). 26 GiB kept run690 waiting 10-16 with ~7 GB held by other apps; def_day level builds peak at 20.2 GiB RSS.
  local mp; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$f" -ge "${FI_MIN691:-20}" ] || { echo "free+inactive $f GiB < ${FI_MIN691:-20}"; return 1; }
  [ "${mp:-0}" -ge "${MP_MIN691:-40}" ] || { echo "memory_pressure free ${mp:-?} % < ${MP_MIN691:-40}"; return 1; }
  [ "$(df -g /System/Volumes/Data | awk 'NR==2{print $4}')" -ge 60 ] || { echo "disk < 60 GiB"; return 1; }
}
last=""; until r="$(gate)"; do [ "$r" != "$last" ] && log "AUTO_R0 waiting: $r"; last="$r"; [ "$(now)" -ge "$CAP" ] && { log "STOP: AUTO_R0 gate not open by CAP"; exit 3; }; sleep 60; done
mkdir -p "$R0W"
if [ "$BRANCH" = SAME-PANEL ]; then
  [ "$(shasum -a 256 $OPS/whitener_union_pruned_def.npz | cut -d' ' -f1)" = "$OPSHA" ] || { log "STOP: 453's operator is not 66610be7"; exit 1; }
  cp -c "$OPS/whitener_union_pruned_def.npz" "$OPS/whitener_union_pruned_def.json" "$R0W/" || { log "STOP: operator copy failed"; exit 1; }
  log "AUTO_R0: SAME-PANEL branch: pre-placed 453's operator 66610be7"
fi
uv run python scripts/october_recipe.py command --bundle "$BUNDLE" --work "$R0W" --name r0def ${AP[@]+"${AP[@]}"} --json > "$W/r0-recipe.json" 2>> "$LOG" \
  || { log "STOP: recipe --json failed"; exit 1; }
ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in json.load(open(sys.argv[1]))["env"].items()))' "$W/r0-recipe.json")"
ARGS="$(python3 -c 'import json,shlex,sys; print(" ".join(shlex.quote(a) for a in json.load(open(sys.argv[1]))["args"]))' "$W/r0-recipe.json")"
log "AUTO_R0: def_day.sh starting: $ENVS $ARGS"
notify "R0 build started (def_day, ~2 h)"
t0=$(date +%s)
[ -n "$OVR" ] && log "AUTO_R0: with the pre-decided override $OVR"
eval "env $OVR $ENVS bash scripts/def_day.sh $ARGS" > "$R0W.nohup" 2>&1; rc=$?
log "AUTO_R0: def_day rc $rc in $(( ($(date +%s) - t0) / 60 )) min; $(grep -E 'DONE|STOP' "$R0W/def_day.log" 2>/dev/null | tail -2 | tr '\n' ' ' | cut -c1-400)"
# 690h: the DNMAG delivery report goes into NEXT.txt (a LOW context is recorded, not rebuilt: the coordinator's rule)
if /usr/bin/grep -q "690h: dnmag" "$R0W/def_day.log" 2>/dev/null; then
  { echo "# 690h DNMAG delivery report (REPORT only; a LOW context is recorded here, never rebuilt):"
    /usr/bin/grep -E "690h: (dnmag|WARNING)" "$R0W/def_day.log" | sed 's/^/#   /'; } >> "$W/NEXT.txt"
  log "690h: $(/usr/bin/grep -oE 'dnmag LOW_CONTEXTS .*' "$R0W/def_day.log" | tail -1)"
fi
if /usr/bin/grep -q "STOP:" "$R0W/def_day.log" 2>/dev/null || [ "$rc" != 0 ]; then
  log "AUTO_R0: R0 def_day STOPPED (rc $rc): a person reruns it in a fresh dir (checklist §C2)"; notify "R0 def_day STOPPED (rc $rc): read $R0W/def_day.log"
else
  notify "R0 def_day DONE: 692 preps it and writes READY; arm1015 R0 submits"
fi
