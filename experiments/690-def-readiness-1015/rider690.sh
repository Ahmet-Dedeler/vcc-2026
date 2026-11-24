#!/usr/bin/env bash
# 690 rider rehearsal (10-15 21:45Z, readiness [4592ee]) for the two October switches that may clear late. Arm it ONLY for
# the ones whose board read cleared:
#   RIDER_LSTR=1  705's lodo_strength theta 0.25 (DEF_DAY_LODO_STRENGTH; 705 reads CLEARED on 10-17 s1)
#   RIDER_K=1     K 3,400 (DEF_DAY_RC_CALLSET_CPM=3000; 711 ADOPT; rider690-K, RIDER_AFTER="rabc2_695.sh rabc3_696.sh")
#   RIDER_HOOK=1  707's day-of dump hook (DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh; 707 reads
#                 TRANSFERS >= +0.0010 on 10-17 s2 / 10-18 s1)
# Both set = one run that also covers the composition (the rider changes per-target amplitude before the dump, the hook
# the strong tier's direction after it).
# What runs: def_day.sh on the R-NEW bundle (work/bundles/new: foreign controls, a NEW 300-target panel, as on Oct 22) with the
# resolved October recipe (`october_recipe.py command --allow-pending`, as R-NEW) + DEF_DAY_ALLOW_DEPTH=1, plus each armed
# rider's export when the recipe does not already carry it (recorded in env.json). Judge = R-NEW's own final meta (same
# bundle, riders off; work/dr_new/dr690n_final.h5ad.meta.json). PASS needs:
#   - LSTR: 522 §8 (a) `def_strength_key522.py --ka` PASS on the frozen tree; def_day logs "LODO strength rider: theta 0.25"
#   - HOOK: def_day logs the hook on the stack dump and "dump replaced by the hook's" (so def_dumphook_check PASSed)
#   - DONE, postbuild all PASS, zero-gene audit 0, builder == HEAD after
#   - config diff vs R-NEW, content-aware (scripts/def_metadiff.py: sameline.path ignored, npz shas compared by content)
#     within {bulk_load, bulk_load_sha256, bulk_ctx_scale} (+ {pert_scale, pert_scale_sha256} with LSTR), and it includes
#     bulk_load_sha256 (the dump moved for real) and, with LSTR, bulk_load (renamed _s0.25; + bulk_sharpen_pt allowed)
# Gates: no 690 HOLD, builder == HEAD, no builder / day-path chain, free+inactive >= 20 GiB and memory_pressure free >= 40 %, disk >= 35 GiB; capped by
# RIDER_CAP (default 10-20 12:00Z). Never preps, never submits. Verdict line "RIDER690 PASS|FAIL ..."; terminal "<ts> END".
# usage: RIDER_LSTR=1 [RIDER_HOOK=1] nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/rider690.sh \
#          >> experiments/690-def-readiness-1015/work/rider690.nohup 2>&1 &
#        DRY=1 RIDER_LSTR=1 RIDER_HOOK=1 bash experiments/690-def-readiness-1015/rider690.sh   (resolve + plan only)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work
if [ -z "${RIDER690_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-rider690.sh"; cp "$REPO/$E/rider690.sh" "$FZ" || exit 1
  RIDER690_FROZEN=1 exec bash "$FZ" "$@"
fi
DRY="${DRY:-0}"; PFX=""; [ "$DRY" = 1 ] && PFX="dry-"
LOG=$W/${PFX}rider690.log; HOLD=$W/HOLD; CAP="${RIDER_CAP:-202610201200}"
LSTR="${RIDER_LSTR:-0}"; HOOK="${RIDER_HOOK:-0}"; KSW="${RIDER_K:-0}"; THETA=0.25; HOOKS=experiments/707-strong-tier-direction/dumphook707.sh; KCPM=3000
# 10-16: run only after these drivers have EXITED (the coordinator's order: R-ABC2 -> rider690 -> R-ABC3 -> rider690-K)
AFTER="${RIDER_AFTER:-rabc2_695.sh}"
B=$W/bundles/new; WD=$W/dr_rider; N=dr690r; REFM=$W/dr_new/dr690n_final.h5ad.meta.json
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'rm -f "$W/.frozen-$$-rider690.sh"; log "END"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
n_builders(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
chains(){ ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh)/' | wc -l | tr -d ' '; }
gate(){
  [ -e "$HOLD" ] && { echo "HOLD file"; return 1; }
  local hb db; hb="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"; db="$(sha scripts/make_submission_pm_cd4.py)"
  [ "$hb" = "$db" ] || { echo "builder on disk ${db:0:16} != HEAD ${hb:0:16}"; return 1; }
  [ "$(n_builders)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  [ "$(chains)" = 0 ] || { echo "a day-path chain is running"; return 1; }
  # skip this run's own process tree: the main pid and its pipeline subshells all carry ".frozen-<this pid>-" in their args
  local a; for a in $AFTER; do ps -axo comm=,args= | awk -v p="$a" -v me=".frozen-$$-" '$1 ~ /bash$/ && index($0, p) && !index($0, me) {n++} END {exit n ? 0 : 1}' && { echo "$a still alive (runs first)"; return 1; }; done
  local f d; f="$(fi_gib)"; d="$(disk_gib)"
  # 10-16: the AGENTS §9 rule (free+inactive >= 20 GiB AND memory_pressure free >= 40 %; no swap clause: since the 10-15
  # reboot there is no swap file). 26 never opened with ~7 GB held by other apps; R-ABC's level builds peaked at 20.2 GiB RSS.
  local mp; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$f" -ge "${FI_MIN690:-20}" ] || { echo "free+inactive ${f} GiB < ${FI_MIN690:-20}"; return 1; }
  [ "${mp:-0}" -ge "${MP_MIN690:-40}" ] || { echo "memory_pressure free ${mp:-?} % < ${MP_MIN690:-40}"; return 1; }
  [ "$d" -ge 35 ] || { echo "disk ${d} GiB < 35"; return 1; }
  echo "builder ${db:0:16} == HEAD, no builder, no chain, free+inactive ${f} GiB, disk ${d} GiB"
}
fail(){ log "RIDER690 FAIL: $*"; exit 1; }
log "rider690 (pid $$, DRY=$DRY): LSTR=$LSTR HOOK=$HOOK K=$KSW on $B; judge $REFM; after [$AFTER]; cap $CAP"
[ "$LSTR" = 1 ] || [ "$HOOK" = 1 ] || [ "$KSW" = 1 ] || fail "arm at least one rider (RIDER_LSTR=1, RIDER_HOOK=1 and/or RIDER_K=1)"
if [ "$HOOK" = 1 ]; then
  git ls-files --error-unmatch "$HOOKS" >/dev/null 2>&1 && [ -z "$(git status --porcelain -- "$HOOKS")" ] || fail "$HOOKS is not committed and clean"
fi
if [ "$DRY" != 1 ]; then
  [ -s "$REFM" ] || fail "no R-NEW final meta at $REFM (run690's R-NEW must have finished first)"
  last=""
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
    [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && fail "gates not met by $CAP ($r)"
    sleep 60
  done
  log "gate PASS: $r"
fi
if [ "$LSTR" = 1 ]; then
  KO="$W/${PFX}ka522"; rm -rf "$KO"; mkdir -p "$KO"
  uv run python experiments/522-october-bulk-key/def_strength_key522.py --ka --out "$KO/ka522.npz" >> "$LOG" 2>&1 || fail "522 §8 (a) def_strength_key522 --ka did not pass"
  log "522 §8 (a): def_strength_key522 --ka PASS"
fi
RJ="$W/${PFX}recipe-rider.json"
uv run python scripts/october_recipe.py command --bundle "$B" --work "$WD" --name "$N" --allow-pending --json > "$RJ" 2>> "$LOG" \
  || fail "october_recipe command failed"
ENVJ="$W/${PFX}env-rider.json"
python3 - "$RJ" "$ENVJ" "$LSTR" "$HOOK" "$THETA" "$HOOKS" "$KSW" "$KCPM" <<'PY' || fail "could not compose the env (see the log)"
import json, sys
rj, out, lstr, hook, theta, hooks, ksw, kcpm = sys.argv[1:]
r = json.load(open(rj)); env = dict(r["env"]); added = {}
# 10-17: a rider the RECIPE already carries (decided ON) must also be flagged, so its evidence is checked, not silently rehearsed
flags = {"DEF_DAY_LODO_STRENGTH": lstr, "DEF_DAY_DUMP_HOOK": hook, "DEF_DAY_RC_CALLSET_CPM": ksw}
unflagged = [k for k, f in flags.items() if k in env and f != "1"]
if unflagged:
    print(f"the recipe carries {unflagged} but their RIDER_* flags are not set: re-arm with them"); sys.exit(1)
want = {}
if lstr == "1":
    want["DEF_DAY_LODO_STRENGTH"] = theta
if hook == "1":
    want["DEF_DAY_DUMP_HOOK"] = hooks
if ksw == "1":
    want["DEF_DAY_RC_CALLSET_CPM"] = kcpm
for k, v in want.items():
    if k in env and env[k] != v:
        print(f"the recipe already sets {k}={env[k]!r}, not {v!r}"); sys.exit(1)
    if k not in env:
        env[k] = v; added[k] = v
env["DEF_DAY_ALLOW_DEPTH"] = "1"
json.dump({"env": env, "args": r["args"], "added_by_rider690": added, "pending_treated_off": r.get("pending_treated_off")},
          open(out, "w"), indent=1)
PY
ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in sorted(json.load(open(sys.argv[1]))["env"].items())))' "$ENVJ")"
ARGS="$(python3 -c 'import json,shlex,sys; print(" ".join(shlex.quote(a) for a in json.load(open(sys.argv[1]))["args"]))' "$ENVJ")"
log "env: $ENVS; args $ARGS; added by rider690: $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["added_by_rider690"])' "$ENVJ")"
if [ "$DRY" = 1 ]; then log "DRY: gates now: $(gate || true); R-NEW meta $( [ -s "$REFM" ] && echo present || echo 'not yet')"; exit 0; fi
[ -e "$WD" ] && mv "$WD" "$WD.aside-$(date -u +%H%M%S)"
mkdir -p "$WD"
t0=$(date +%s)
eval "env $ENVS bash scripts/def_day.sh $ARGS" > "$WD.nohup" 2>&1; rc=$?
DL="$WD/def_day.log"
log "def_day rc $rc in $(( ($(date +%s) - t0) / 60 )) min; $(/usr/bin/grep -E 'LODO strength rider|690g: dump|CONVERGED at|NOT converged|final:|postbuild PASS|DONE' "$DL" 2>/dev/null | tr '\n' ' ' | cut -c1-900)"
[ "$(sha scripts/make_submission_pm_cd4.py)" = "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)" ] || fail "builder != HEAD after def_day"
[ "$rc" = 0 ] && /usr/bin/grep -q "DONE. Next" "$DL" || fail "def_day did not finish (rc $rc)"
[ "$LSTR" = 1 ] && { /usr/bin/grep -q "LODO strength rider: theta $THETA" "$DL" || fail "no 'LODO strength rider: theta $THETA' in def_day.log"; }
if [ "$HOOK" = 1 ]; then
  /usr/bin/grep -q "690g: dump hook $HOOKS" "$DL" && /usr/bin/grep -q "690g: dump replaced by the hook's" "$DL" || fail "the dump hook did not run and validate"
fi
FM="$(ls "$WD"/${N}_final.h5ad.meta.json 2>/dev/null)"; [ -s "$FM" ] || fail "no final meta"
# the config diff vs R-NEW, content-aware (scripts/def_metadiff.py, 694 v2 rule): a new work dir changes the npz FILE shas
# (bulk_load_sha256, sameline.sha256) through the path inside them; only a content change counts. The riders must move
# the dump for real (bulk_load_sha256 is needed) and, with LSTR, the key (pert_scale_sha256).
# per rider (10-16): the hook moves the dump's content; θ renames it; K moves only the call set (the dump stays)
ALLOWK=""; NEEDK=""
[ "$HOOK" = 1 ] && { ALLOWK="$ALLOWK,bulk_load_sha256,bulk_ctx_scale"; NEEDK="$NEEDK,bulk_load_sha256"; }
[ "$KSW" = 1 ] && { ALLOWK="$ALLOWK,rc_callset_cpm"; NEEDK="$NEEDK,rc_callset_cpm"; }
# 10-16: θ rewrites only the DUMP's key (same LODO id, the strength file), so the meta shows it as bulk_load renamed _s0.25 (+sha),
# the tercilec record (bulk_sharpen_pt: carried_key_sha256) and the re-level; the level builds keep their pert_scale.
[ "$LSTR" = 1 ] && { ALLOWK="$ALLOWK,bulk_load,bulk_load_sha256,bulk_ctx_scale,pert_scale,pert_scale_sha256,bulk_sharpen_pt"; NEEDK="$NEEDK,bulk_load"; }
ALLOWK="${ALLOWK#,}"; NEEDK="${NEEDK#,}"
uv run python scripts/def_metadiff.py "$FM" "$WD" "$REFM" "$W/dr_new" --allow "$ALLOWK" --need "$NEEDK" --json "$W/rider690.metadiff.json" >> "$LOG" 2>&1 \
  || fail "config diff vs R-NEW (content-aware) is not within the riders' keys (see the [metadiff] lines above)"
python3 - "$WD" >> "$LOG" 2>&1 <<'PY' || fail "postbuild / zero-gene audit (see the lines above)"
import json, re, sys
from pathlib import Path
wd = Path(sys.argv[1])
pb = json.loads((wd / "postbuild.json").read_text())
bad_pb = [x.get("check") for x in pb.get("checks", []) if not x.get("pass")]
au = (wd / "audit428.txt").read_text().splitlines()
zero = bool(au) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", au[0]) is not None
print(f"[rider690] postbuild failing {bad_pb}; zero-gene audit {'0' if zero else 'NOT 0'}")
sys.exit(0 if not bad_pb and zero else 1)
PY
rm -f "$WD"/*.h5ad
log "RIDER690 PASS: LSTR=$LSTR HOOK=$HOOK rehearsed on the NEW-panel bundle; meta differs from R-NEW only where the riders act"
