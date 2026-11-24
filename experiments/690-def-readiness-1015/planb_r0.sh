#!/usr/bin/env bash
# 690 plan B for R0 (readiness [401d9b], 2026-10-22; asked by the coordinator after R0's lvl0 and lvl1 were both deleted past
# their cliffs: lvl0 E 1.5555, lvl1 D 1.0151 / F 1.0800). Run it ONLY if R0's def_day in experiments/692-def-r0/work/r0 ended
# with "STOP: no build with every rho <= 1; do not ship" (lvl2..lvl4 all past a cliff). OCT22-CHECKLIST §C2: that dir is spent;
# rerun in a FRESH dir and re-arm 692 + arm1015 on it.
#
# What it does (each step fail-closed; DRY=1 prints the commands and touches nothing):
#   1. checks: r0's def_day.log has that STOP line AND its trap line "builder after restore:", no def_day.sh runs, the builder on
#      disk == HEAD (e99bbd4a), the four other frozen files at their shas, r0b does not exist, disk >= 35 GiB;
#   2. stops wait708 (its log is no one's ready-log: OCT22-CHECKLIST 4b says stopping it is safe). It must be re-armed only AFTER
#      the coordinator's new arm1015 R0 has logged its ARMED line, or it would read the old arm's NOT SUBMITTING END as "R0 done";
#   3. launches def_day in experiments/692-def-r0/work/r0b: NEXT.txt's recipe + DEF_DAY_CD4_THRESHOLDS=relaxed (as R0) +
#      DEF_DAY_RHO_SCALE=$PLANB_RHO_SCALE (default 0.98: FINAL-PHASE row 2's pre-studied value; target rho* x 0.98 = 0.9597,
#      ~4 % below the cliff in rho, ~1 % in E's scale at its local exponent ~4.2) and max_levels $PLANB_MAXL (default 6).
#      Same name r0def. lvl0 and lvl1 replay r0's (deterministic; the first step is ntclevel's EV scale, not rho-scaled), so the
#      new information starts at lvl2: expect ~75 min per build, READY ~5-6 h after launch;
#   4. arms r0ready692 on r0b with OUT692=experiments/692-def-r0/work/r0b-ready (a fresh ready log: arm1015 stops on any
#      NOT READY line, and r0's ready log will carry one).
# Then (the coordinator's part, printed at the end): re-arm arm1015 R0 + name_guard on the new ready / run logs, then wait708.
#   PLANB_RHO_SCALE=0.98|off PLANB_MAXL=6 DRY=1 bash experiments/690-def-readiness-1015/planb_r0.sh
#
# PLANB_MODE=lowdose (added 10-22 05:30Z, after the real case): R0 CONVERGED at lvl3 (rho D 0.9762 / E 0.9764 / F 0.9749) and
# then STOPPED on ONE postbuild row, 554's emitter steering dose (D median prior_lambda 0.8975 < floor 0.9; E/F 1.0; D's MDE
# 0.954; missed calls 1/24/39). 554 RESULTS prescribes this case: a context at 0.85-0.9 is a ~10-15 % dose loss; "stop, look,
# and override with DEF_DAY_ALLOW_LOW_DOSE=1 if nothing else is wrong". The override exists only as a def_day env, so this
# mode reruns the SAME recipe in the fresh dir with DEF_DAY_ALLOW_LOW_DOSE=1 (rho-scale off, max_levels 4). The chain is
# deterministic (R-KA2 rebuilt R-ABC2 byte for byte in another dir), so r0b must reproduce r0's lvl3 output fa832988…:
# ka_r0b.sh (launched here, detached) touches experiments/670-coord-1015/work/NOGO-R0 on any other final output sha.
#   PLANB_MODE=lowdose DRY=1 bash experiments/690-def-readiness-1015/planb_r0.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
DRY="${DRY:-0}"; MODE="${PLANB_MODE:-rho}"
case "$MODE" in
  rho) RS="${PLANB_RHO_SCALE:-0.98}"; MAXL="${PLANB_MAXL:-6}"; LDE="" ;;
  lowdose) RS="${PLANB_RHO_SCALE:-off}"; MAXL="${PLANB_MAXL:-4}"; LDE="DEF_DAY_ALLOW_LOW_DOSE=1" ;;
  *) echo "PLANB_MODE must be rho or lowdose (got $MODE)" >&2; exit 2 ;;
esac
KA_SHA=fa83298813978cd8210d08f320f62d63d9501a630ca102e6310f1dfb95c07799   # r0's lvl3 = r0def_final output (05:17:30Z)
R0=experiments/692-def-r0/work/r0; R0B=experiments/692-def-r0/work/r0b; OUTB=experiments/692-def-r0/work/r0b-ready
LOG=experiments/690-def-readiness-1015/work/planb_r0.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
die(){ log "PLANB STOP: $*"; exit 1; }
run(){ if [ "$DRY" = 1 ]; then log "DRY: $*"; else log "RUN: $*"; eval "$@"; fi; }
sha8(){ shasum -a 256 "$1" | cut -c1-8; }

log "planb_r0 (pid $$, DRY=$DRY, mode $MODE): rho-scale $RS, max_levels $MAXL, fresh dir $R0B${LDE:+, $LDE}"
if [ "$RS" = off ]; then RSE=""; log "rho-scale off: the decided recipe's target rho* x1 (no DEF_DAY_RHO_SCALE in the env)"
else python3 -c 'import sys; s=float(sys.argv[1]); sys.exit(0 if 0.9<=s<=1.1 else 1)' "$RS" || die "PLANB_RHO_SCALE $RS outside [0.9, 1.1] (or off)"; RSE="DEF_DAY_RHO_SCALE=$RS"; fi
# 1. checks
if [ "$MODE" = rho ]; then
  /usr/bin/grep -q "STOP: no build with every rho <= 1" "$R0/def_day.log" || [ "$DRY" = 1 ] \
    || die "$R0/def_day.log has no 'no build with every rho <= 1' STOP: plan B (rho) is not for this case"
else
  /usr/bin/grep -q "STOP: postbuild FAILED" "$R0/def_day.log" || die "$R0/def_day.log has no 'postbuild FAILED' STOP: lowdose mode is not for this case"
  # the ONLY failing postbuild row must be 554's dose row, every failing dose floor a lambda in [0.85, 0.9), missed calls fine
  python3 - "$R0/postbuild.json" <<'PY' || die "postbuild.json fails something besides a 0.85-0.9 dose (lowdose mode refuses)"
import json, sys
d = json.load(open(sys.argv[1]))
bad = [c for c in d["checks"] if not c["pass"]]
assert len(bad) == 1 and "(554)" in bad[0]["check"], [c["check"] for c in bad]
g = d["dose_guard"]
for ctx, v in g["per_context"].items():
    for k in ("prior_lambda", "mde_lambda"):
        assert v[k] >= 0.85, (ctx, k, v[k])
    assert v["call_set_missed"] <= g["floors"]["call_set_missed"], (ctx, v["call_set_missed"])
print("lowdose precondition: the only failing row is 554's; every lambda >= 0.85")
PY
  [ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$R0/r0def_final.h5ad.meta.json")" = "$KA_SHA" ] \
    || die "r0's final meta does not record output $KA_SHA"
fi
n7="$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /708-def-panel-repricing\/day708\.sh/' | wc -l | tr -d ' ')"
[ "$n7" = 0 ] || [ "$DRY" = 1 ] || die "day708 is running: launch after its DECISIONS WRITTEN (the coordinator's order, 10-22)"
/usr/bin/grep -q "builder after restore:" "$R0/def_day.log" || [ "$DRY" = 1 ] || die "r0's def_day has not restored the builder yet"
n="$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /scripts\/def_day\.sh/' | wc -l | tr -d ' ')"
[ "$n" = 0 ] || [ "$DRY" = 1 ] || die "a def_day.sh is running"
[ "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -c1-8)" = e99bbd4a ] || die "HEAD's builder is not e99bbd4a"
[ "$(sha8 scripts/make_submission_pm_cd4.py)" = e99bbd4a ] || [ "$DRY" = 1 ] || die "builder on disk is not e99bbd4a (pinned? restore first)"
for fs in docs/october-decisions.json:04be2af3 scripts/def_day.sh:11ff102e src/vcc2026/rc_emit_k.py:0f56209d scripts/def_day_helpers.py:5354e449; do
  [ "$(sha8 "${fs%%:*}")" = "${fs##*:}" ] || die "${fs%%:*} is not ${fs##*:}"
done
[ -e "$R0B" ] && die "$R0B exists (a rerun goes in a FRESH dir)"
[ -e "$OUTB" ] && die "$OUTB exists"
d="$(df -g /System/Volumes/Data | awk 'NR==2{print $4}')"; [ "$d" -ge 35 ] || die "disk $d GiB < 35"
log "checks PASS (disk $d GiB)"
# 2. stop wait708
WP="$(ps -axo pid=,comm=,args= | awk '$2 ~ /bash$/ && /wait708\.sh/ {print $1}' | head -1)"
if [ -n "$WP" ]; then run "kill $WP"; log "wait708 (pid $WP) stopped; re-arm it only after the new arm1015 R0's ARMED line"; else log "no wait708 running"; fi
# 3. def_day in the fresh dir (the recipe re-resolved from the live decisions file, as NEXT.txt; the two day knobs added)
RJ=experiments/690-def-readiness-1015/work/def-day/r0b-recipe.json
env -u OCTOBER_DECISIONS uv run python scripts/october_recipe.py command --bundle data/challenge-def --work "$R0B" --name r0def --json > "$RJ" \
  || die "october_recipe command failed"
ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in json.load(open(sys.argv[1]))["env"].items()))' "$RJ")"
ARGS="$(python3 -c 'import json,shlex,sys; a=json.load(open(sys.argv[1]))["args"]; a[-1]=sys.argv[2]; print(" ".join(shlex.quote(x) for x in a))' "$RJ" "$MAXL")"
case "$ARGS" in "data/challenge-def $R0B r0def 1.75 0.5 noprog $MAXL") ;; *) die "unexpected recipe args: $ARGS";; esac
run "env -u OCTOBER_DECISIONS DEF_DAY_CD4_THRESHOLDS=relaxed $RSE $LDE $ENVS nohup bash scripts/detach.sh bash scripts/def_day.sh $ARGS > $R0B.nohup 2>&1 &"
[ "$DRY" = 1 ] || { sleep 20; /usr/bin/grep -q "re-level\|start: bundle" "$R0B/def_day.log" 2>/dev/null || die "def_day in $R0B did not start (see $R0B.nohup)"; }
# 3b. lowdose: the known-answer guard (r0b must rebuild r0's output byte for byte; else NOGO-R0 before the submit)
[ "$MODE" = lowdose ] && run "nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/ka_r0b.sh $KA_SHA > experiments/690-def-readiness-1015/work/ka_r0b.nohup 2>&1 &"
# 4. 692 on the fresh dir, with a fresh ready log
run "mkdir -p $OUTB"
run "OUT692=$OUTB nohup bash scripts/detach.sh bash experiments/692-def-r0/r0ready692.sh data/challenge-def $R0B r0def 202610270000 > $OUTB/r0ready692.nohup 2>&1 &"
log "launched. The coordinator's part (its scripts; not run here):"
log "  arm1015 R0 on ready-log $OUTB/r0def.ready.log, run-log $OUTB/r0def.run.log, pkg submissions/r0def.vcc, the same name and spec as 5181"
log "  (the old arm never submitted), plus name_guard with the matching spec; then re-arm wait708:"
log "  WAIT708_R710=POSITIVE nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/wait708.sh > experiments/690-def-readiness-1015/work/wait708.nohup 2>&1 &"
log "  (day708 reads r0's whitener_union_pruned_def.npz and cd4_compact/, which r0 still holds; keep r0)"
log "END"
