#!/usr/bin/env bash
# def_day.sh — the October D/E/F day in one gated script (2026-09-25, coordinator; generalises 448's run448.sh).
# 476 (2026-09-27): 485's review fixes merged onto 474 (481's patches). Fail closed on driver rc, dump rc, audit and prep;
# a bad secant step ships the nearest safe build instead of aborting; lvlC must sit inside the 1% band; past_cliff crashes
# are not "safe". Needs the driver with 485's dump-assert skip (def_pipeline.py): without it the driver exits 1 after
# every successful dump and this script's dump rc check stops the day.
#
#   bash scripts/def_day.sh <bundle> <work_dir> <name> [q=1.75] [sameline_w=auto] [sameline_flags=pt:noprog] [max_levels=4]
#   (sameline_w `auto`, builder >= 4feae27d: each matched context gets w from its bank's quality class, 2 strong / 0.5 weak,
#    written by def_sameline.py; a number forces one global w)
#   DEF_DAY_LODO=1 (460, default off): the stack dump is scaled by 415's LODO-V1 key re-fitted on the day's panel
#    (446bj2's key on A/B/C) instead of 410's; see stage 2b.
#   DEF_DAY_PROGRAMME_SOURCE=<bank> (461, default off; e.g. h1_2025_full = Arc's H1 from all public cells, 458d): every
#    context's block prior comes from that bank (def_sameline.py --programme-source). Needs sameline_flags WITHOUT noprog
#    (pass `pt`); the script refuses `noprog` with it, because the builder would then ignore the programme.
#   DEF_DAY_PRIOR_GATE=q<NN> (468, default off; e.g. q33): with DEF_DAY_PROGRAMME_SOURCE, each context's design-weakest
#    NN % of targets keep CLEAN4 instead of the programme (the builder's VCC_DE_PRIOR_GATE, builder >= 693fde15; 467 on H1).
#    Refused without DEF_DAY_PROGRAMME_SOURCE. Exported for the lvl0 and re-level builds only (never for the stack dump).
#   DEF_DAY_LODO_HI_PCT=<P> (497, default off; needs DEF_DAY_LODO=1): the LODO key's cap is the panel's own P-th
#    percentile of carried key / median (def_lodo_key.py --hi-pct; 543/546), pinned as lodo_lo0.5_pct_def.
#    P 96.64697522769474 is exactly cap 6 on A/B/C (495's artifact); the rounded 96.647 is not.
#   DEF_DAY_SHARPEN_PT=<tercile|stercile>:<q>,<q>,<q> (497, default off; e.g. stercile:2.5,2.0,1.75 = 479's S_P2): the
#    level builds sharpen each target's bulk with its tercile's q (VCC_BULK_SHARPEN_Q's per-target form, builder >=
#    d09a986c) instead of the scalar q. The stack dump is unchanged (q is a render knob).
#   DEF_DAY_DE_MAG=<variant> (502, default off; e.g. dn_k250_zs): the zs exact-U pass (VCC_DE_MAG) on the level builds.
#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).
#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.
#    The default restores, which is what the 452 rehearsal checks. Pins are re-applied by re-running the pin stages.
#
# Stages (each logged to <work>/def_day.log; a failed driver, audit, or prep stops the script):
#   1. the driver's light and pin stages with a bare env: preflight, identity, cd4-extract, census, the day's whitener,
#      the day's 410 key, repin-frozen, pin-artifacts, repin-panel-support (437/439);
#   2. the same-line artifact from the identity call (451: `scripts/def_sameline.py`). If any context is matched,
#      VCC_SAMELINE=<npz>:<sha>:<w>:<flags> is exported for every build below. Unmatched: nothing changes;
#   3. the stack dump, H1 x2 / CD4 x1 (438's October default: VCC_H1_UNION_G=0.56, CD4 at DEF_ENV's default);
#   4. lvl0, the un-levelled LOAD build at q, with SHIPPED_ENV's coarse MDE and zero repair (driver defaults);
#   5. the re-level loop (448): ntclevel417 (sd 0.01) on each build, then `relevel_iter.py` over every build so far,
#      until CONVERGED (every context within 1 % of rho* and none above 1) or max_levels levelled builds. A build with any
#      rho > 1 is never the final one. A proposal outside [0.7, 1.4], or a degenerate secant (relevel exit 2), ships the
#      nearest already-built safe artifact instead of stopping;
#   6. audits on the final build: ntclevel, zero-gene audit (stop unless the zero count is 0), and a `vcc prep --dry-run`
#      through the driver (stop unless that stage exits 0). lfc < -5 is logged and not gated (452 shipped with 7).
# It stops BEFORE the real `vcc prep` and `vcc submit`: those stay manual (DEF-RUNBOOK §2's tail, AGENTS §6.2, §14).
# The final artifact is <work>/<name>_final.h5ad (a hard link to the converged build). Its VCC_BULK_CTX_SCALE string is
# in <work>/def_day.final.env. The link is made only after the meta's q and scale match that string.
#
#   DEF_DAY_CD4_THRESHOLDS=relaxed (481 F1, default `default`): the CD4 extraction runs without 047's gates (275 per
#    condition, 290 union). Arc's other 300-target panel (VCC 2025 H1) covers 262/265/262, union 267, so the gates fire
#    on an Arc-style panel; even A/B/C cleared the union gate by one target. The census's own floors (union >= 120,
#    complete >= 100) still gate. Check the mapping first (244's mapping244.py), as DEF-RUNBOOK section 0.0 says.
#   DEF_DAY_INCUMBENT_META=<meta.json> (481 F5): the shipped A/B/C arm the final artifact's config is diffed against by
#    def_postbuild.py before the builder is restored (default: 446b's meta).
#
# Rehearse on a synthetic bundle first (e.g. experiments/238-essential-free-panel/work/ef-bundle-r).
# A rerun needs a fresh work dir: this script deletes a tag's h5ad before rebuilding it, and refuses to
# treat a pre-existing file as the build it just asked for.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
B="${1:?bundle}"; W="${2:?work dir}"; N="${3:?name}"; Q="${4:-1.75}"; SLW="${5:-auto}"; SLF="${6-pt:noprog}"; MAXL="${7:-4}"
mkdir -p "$W"; LOG="$W/def_day.log"
CTXS="$(python3 -c "import json, sys; print(','.join(json.load(open(sys.argv[1]))['contexts']))" "$B/manifest.json")" || exit 1
UNIT="$(python3 -c "import sys; print(','.join(c + '=1' for c in sys.argv[1].split(',')))" "$CTXS")"
CD4T="${DEF_DAY_CD4_THRESHOLDS:-default}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
die(){ log "STOP: $*"; exit 1; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "start: bundle $B work $W name $N q $Q sameline w $SLW flags $SLF; builder ${BSHA0:0:16} driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"
# A reboot skips this script's EXIT trap, so the builder can still be the LODO pin. A fresh work
# dir would record that pin as BSHA0 and the trap would then "restore" it (526). Refuse.
_HEADB="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
if [ "$BSHA0" != "$_HEADB" ]; then
  die "builder on disk ${BSHA0:0:16} is not HEAD ${_HEADB:0:16}. Restore before a fresh work dir: uv run python scripts/def_pin_panel_artifacts.py --restore --record <previous-work>/pin-artifacts.json"
fi
grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_SAMELINE (apply 451's patch first)"
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
case "$CD4T" in
  default) ;;
  relaxed) log "CD4 extract thresholds relaxed (DEF_DAY_CD4_THRESHOLDS=relaxed): check the mapping first (DEF-RUNBOOK 0.0)" ;;
  *) die "DEF_DAY_CD4_THRESHOLDS must be relaxed or unset (got $CD4T)" ;;
esac
# 497: the October forms of 495 (percentile cap) and 479 (per-target sharpening); both default off
if [ -n "${DEF_DAY_LODO_HI_PCT:-}" ]; then
  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_HI_PCT=$DEF_DAY_LODO_HI_PCT needs DEF_DAY_LODO=1"
  python3 -c "import sys; p = float(sys.argv[1]); sys.exit(0 if 50.0 < p < 100.0 else 1)" "$DEF_DAY_LODO_HI_PCT" 2>/dev/null \
    || die "DEF_DAY_LODO_HI_PCT must be a percentile in (50, 100) (got $DEF_DAY_LODO_HI_PCT)"
  grep -q -- "--hi-pct" scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi-pct (497)"
fi
# 502: the DE-channel riders (validated here, exported for the level builds after the dump)
if [ -n "${DEF_DAY_DE_MAG:-}" ]; then
  case "$DEF_DAY_DE_MAG" in *[!a-z0-9_]*) die "DEF_DAY_DE_MAG must be a variant name like dn_k250_zs (got $DEF_DAY_DE_MAG)" ;; esac
  grep -q "VCC_DE_MAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MAG"
fi
if [ -n "${DEF_DAY_MDE_SIGN:-}" ]; then
  case "$DEF_DAY_MDE_SIGN" in clean4:*) ;; *) die "DEF_DAY_MDE_SIGN must be clean4:<tau> (got $DEF_DAY_MDE_SIGN)" ;; esac
  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"
fi
SQ="$Q"
if [ -n "${DEF_DAY_SHARPEN_PT:-}" ]; then
  case "$DEF_DAY_SHARPEN_PT" in tercile:*,*,*|stercile:*,*,*) ;; *) die "DEF_DAY_SHARPEN_PT must be tercile:<q>,<q>,<q> or stercile:<q>,<q>,<q> (got $DEF_DAY_SHARPEN_PT)" ;; esac
  grep -q "def sharpen_q_by_tercile" scripts/make_submission_pm_cd4.py || die "the builder has no per-target sharpening (479)"
  SQ="$DEF_DAY_SHARPEN_PT"
  log "per-target sharpening for the level builds: VCC_BULK_SHARPEN_Q=$SQ (scalar q $Q not used)"
fi
restore(){
  local rc=$?
  if [ "${DEF_DAY_KEEP_PINS:-0}" = "1" ]; then
    log "DEF_DAY_KEEP_PINS=1: builder left pinned $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (pre-pin ${BSHA0:0:16}); commit the pins, then restore from $W/pin-artifacts.json"
    exit "$rc"
  fi
  if [ -f "$W/pin-artifacts.json" ]; then
    if ! uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1; then
      log "STOP: builder restore command failed"
      rc=1
    fi
  fi
  local now
  now="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  log "builder after restore: ${now:0:16} (want ${BSHA0:0:16})"
  if [ "$now" != "$BSHA0" ]; then
    log "STOP: restored builder sha does not match the pre-pin sha"
    rc=1
  fi
  exit "$rc"
}
trap restore EXIT

# 1. setup
"${DRV[@]}" --name "$N-setup" --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support \
    --cd4-thresholds "$CD4T" > "$W/driver-setup.log" 2>&1
if ! grep -q "^--- \[repin-panel-support\].*rc=0" "$W/driver-setup.log"; then
  grep -qE "condition coverage below|panel targets map across any condition" "$W/logs/cd4-extract.log" 2>/dev/null \
    && die "CD4 extraction gate fired ($(grep -E 'panel coverage|only [0-9]+ panel targets|coverage below' "$W/logs/cd4-extract.log" | head -2 | tr '\n' ' ')). Check the mapping (experiments/244-october-rehearsal/mapping244.py), then rerun with DEF_DAY_CD4_THRESHOLDS=relaxed in a fresh work dir"
  die "setup did not complete (see driver-setup.log)"
fi
log "setup: $(grep -E '^--- \[' "$W/driver-setup.log" | awk '{print $2, $3}' | tr '\n' ' ')"
[ -s "$W/identity.json" ] || die "no identity.json"
log "identity: $(python3 -c "import json; d=json.load(open('$W/identity.json')); print(d['call'], 'unanimous', d['call_unanimous_across_methods'])")"
# Context labels come from the preflight the setup stage just wrote. Hardcoding D,E,F ships a
# mismatch the moment the manifest uses other names: lvl0 has no scale and succeeds, then ntclevel
# looks for context_D.h5ad. isidentifier() rejects anything that is not a bare name.
CTXP="$(python3 -c 'import json,sys
cs=json.load(open(sys.argv[1]))["contexts"]
assert isinstance(cs, list) and len(cs)==3 and all(isinstance(c, str) and c.isidentifier() for c in cs), cs
print(",".join(cs))' "$W/preflight.json")" || die "could not read three context names from $W/preflight.json"
[ "$CTXP" = "$CTXS" ] || die "preflight contexts $CTXP differ from the manifest's $CTXS"
log "contexts $CTXS; unit factors $UNIT"

# 2. same-line
SL=""
if [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ]; then
  case ":$SLF:" in *:noprog:*) die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE needs sameline_flags without noprog (got $SLF)";; esac
  log "programme-source $DEF_DAY_PROGRAMME_SOURCE: every context's block prior from it (flags $SLF)"
fi
if [ -n "${DEF_DAY_PRIOR_GATE:-}" ]; then
  [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ] || die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE gates a programme: it needs DEF_DAY_PROGRAMME_SOURCE"
  grep -q "VCC_DE_PRIOR_GATE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_PRIOR_GATE (apply 468's patch first)"
  log "prior gate $DEF_DAY_PRIOR_GATE: each context's design-weakest targets keep CLEAN4"
fi
# DEF_DAY_FORCE="D=<CCLE sample>,..." forces the same-line call (rehearsals only; never on the day)
uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/sameline.npz" \
    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} ${DEF_DAY_PROGRAMME_SOURCE:+--programme-source "$DEF_DAY_PROGRAMME_SOURCE"} \
    >> "$LOG" 2>&1 || die "def_sameline.py failed"
[ -n "${DEF_DAY_FORCE:-}" ] && log "REHEARSAL: same-line call forced to $DEF_DAY_FORCE"
NCTX="$(uv run python scripts/def_day_helpers.py nctx "$W/sameline.npz")" || die "nctx could not read $W/sameline.npz"
case "$NCTX" in
  ''|*[!0-9]*) die "nctx printed '$NCTX'" ;;
esac
if [ "$NCTX" -gt 0 ]; then
  SL="$W/sameline.npz:$(shasum -a 256 "$W/sameline.npz" | cut -d' ' -f1):$SLW${SLF:+:$SLF}"
  log "same-line mode ON for $NCTX context(s): VCC_SAMELINE=${SL##*/}"
else
  if [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ]; then
    die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE but the same-line artifact switches on 0 contexts"
  fi
  log "same-line mode off (no unanimous registered line)"
fi

# 2b. DEF_DAY_LODO=1 (460; default off): the stack dump is scaled by 415's LODO-V1 key at clip (0.5, 3), re-fitted on
#     this panel, the day's operator and CD4 banks (`scripts/def_lodo_key.py`; its --ka reproduces A/B/C's 8c82bef3). A
#     second pin registers it as `lodo_lo0.5_hi3_def` next to the 410 key, and only the dump build selects it, exactly as
#     446bj2 did on A/B/C (config diff vs the default path: bulk_load, bulk_load_sha256). The restore trap undoes both pins.
LODO=""
if [ "${DEF_DAY_LODO:-0}" = "1" ]; then
  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"; LID=lodo_lo0.5_hi3_def
  [ -n "${DEF_DAY_LODO_HI_PCT:-}" ] && { LK="$W/pert_scale_lodo_def_lo0.5_pct.npz"; LID=lodo_lo0.5_pct_def; }
  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \
      ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"
  uv run python scripts/def_pin_panel_artifacts.py --whitener "$W/whitener_union_pruned_def.npz" \
      --whitener-metadata "$W/whitener_union_pruned_def.json" --whitener-id union_pruned_def --cd4-root "$W/cd4_compact" \
      --pert-scale "$LK" --pert-scale-id "$LID" --apply --out "$W/pin-lodo.json" >> "$LOG" 2>&1 \
      || die "LODO key pin failed"
  LODO="$LID"
  log "LODO key ON for the stack dump: $LODO $(shasum -a 256 "$LK" | cut -c1-16); cap $(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["hi"])' "${LK%.npz}.json")"
fi

# 3. the stack dump. The subshell keeps VCC_H1_UNION_G / VCC_BULK_DUMP / VCC_PERT_SCALE off the level builds.
#    Remove any previous npz first: a failed retry must not keep the last successful dump.
DUMP="$W/dump_h1x2_cd4x1.npz"
rm -f "$DUMP"
(
  export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"
  [ -n "$LODO" ] && export VCC_PERT_SCALE="$LODO"
  "${DRV[@]}" --name "$N-dump" --stages build > "$W/driver-dump.log" 2>&1
)
rc=$?
if [ "$rc" -ne 0 ] || [ ! -s "$DUMP" ]; then
  die "dump stage failed rc=$rc (see $W/driver-dump.log)"
fi
# The full builder log, not the driver's 25-line tail. time -l appends ~18 lines after the
# dump message, so the tail still contains it today; the stage log is the one that cannot drop it.
grep -q '\[bulk-dump\] wrote' "$W/logs/build.log" || die "dump stage log has no [bulk-dump] wrote line"
DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "dump ${DSHA:0:16}"
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"
[ -n "${DEF_DAY_DE_MAG:-}" ] && { export VCC_DE_MAG="$DEF_DAY_DE_MAG"; log "DE rider for the level builds: VCC_DE_MAG=$VCC_DE_MAG"; }
[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$DEF_DAY_MDE_SIGN"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }
[ -n "$SL" ] && export VCC_SAMELINE="$SL"
if [ -n "${DEF_DAY_PRIOR_GATE:-}" ]; then
  [ -n "$SL" ] || die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE, but same-line mode is off: the programme never reached the builder"
  export VCC_DE_PRIOR_GATE="$DEF_DAY_PRIOR_GATE"
fi

# 4-5. lvl0, then the re-level loop
POINTS=()
FACT="$UNIT"
build(){ local tag=$1 fact=$2 h5 rc pc
  if [ "$fact" = "$UNIT" ]; then unset VCC_BULK_CTX_SCALE; else export VCC_BULK_CTX_SCALE="$fact"; fi
  h5="$W/$N-$tag.h5ad"
  # The driver refuses to overwrite. Delete this tag's previous artifact so a retry rebuilds it,
  # then require the driver's own exit code. A pre-existing file is not evidence the new build ran.
  rm -f "$h5" "$h5.meta.json"
  "${DRV[@]}" --name "$N-$tag" --stages build > "$W/driver-$tag.log" 2>&1
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$h5" ]; then
    rm -f "$h5" "$h5.meta.json"
    die "$tag build failed rc=$rc: $(tail -3 "$W/driver-$tag.log" 2>/dev/null)"
  fi
  uv run python experiments/417-bulk-sharpen/ntclevel417.py "$h5" --bundle "$B" --contexts "$CTXS" --sd 0.01 >> "$LOG" 2>&1 \
    || die "ntclevel failed on $tag"
  POINTS+=("$W/$N-$tag.ntclevel.json:$fact")
  log "$tag ($fact) rho: $(uv run python scripts/def_day_helpers.py rho "$W/$N-$tag.ntclevel.json")"
  # past_cliff: 0 = a context is past its cliff (delete), 1 = safe, anything else = unreadable (stop).
  uv run python scripts/def_day_helpers.py past_cliff "$W/$N-$tag.ntclevel.json"
  pc=$?
  if [ "$pc" -eq 0 ]; then
    rm -f "$h5" "$h5.meta.json"; log "$tag has a context past its cliff: artifact deleted (never shippable)"
  elif [ "$pc" -ne 1 ]; then
    die "past_cliff failed on $tag rc=$pc; not treating the artifact as safe"
  fi; }
same_f(){ uv run python scripts/def_day_helpers.py same "$1" "$2"; }
build lvl0 "$FACT"
FINAL=""
for L in $(seq 1 "$MAXL"); do
  OUT="$(uv run python experiments/448-def-rehearsal-q175/relevel_iter.py "${POINTS[@]}" 2>&1)"
  rc=$?
  echo "$OUT" >> "$LOG"
  if [ "$rc" -ne 0 ]; then
    log "relevel_iter rc=$rc; shipping the nearest safe build"
    break
  fi
  NEXT="$(echo "$OUT" | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2-)"
  if [ -z "$NEXT" ]; then
    log "relevel_iter printed no factors; shipping the nearest safe build"
    break
  fi
  if echo "$OUT" | grep -q "^CONVERGED"; then
    # the converged factors are one of the builds already made: find it
    for P in "${POINTS[@]}"; do same_f "${P##*:}" "$NEXT" && [ -s "${P%%.ntclevel.json:*}.h5ad" ] && FINAL="${P%%.ntclevel.json:*}"; done
    if [ -z "$FINAL" ]; then
      # each context converged in a different build: build the combination once, and take it only
      # when the new measurement is inside the same 1% band relevel required of the pieces.
      log "CONVERGED per context at $NEXT, but no single build carries all of them: building it"
      build lvlC "$NEXT"
      if [ -s "$W/$N-lvlC.h5ad" ] && uv run python scripts/def_day_helpers.py inband "$W/$N-lvlC.ntclevel.json"; then
        FINAL="$W/$N-lvlC"
      else
        log "lvlC is missing or outside the 1% band; not labelling it CONVERGED"
        rm -f "$W/$N-lvlC.h5ad" "$W/$N-lvlC.h5ad.meta.json"
      fi
    fi
    if [ -n "$FINAL" ]; then
      log "CONVERGED at $NEXT after $((L - 1)) levelled build(s)"
    else
      log "CONVERGED factors were not a shippable build; shipping the nearest safe build"
    fi
    break
  fi
  if ! uv run python scripts/def_day_helpers.py inbound "$NEXT"; then
    log "next factors $NEXT leave the builder bound [0.7, 1.4]; shipping the nearest safe build"
    break
  fi
  build "lvl$L" "$NEXT"
done
if [ -z "$FINAL" ]; then
  # not converged: ship the last build whose contexts are all below their cliffs, the nearest to rho* among those.
  # rho* is computed in the helper from the same q and sd as relevel_iter.py.
  FINAL="$(uv run python scripts/def_day_helpers.py best_safe "${POINTS[@]}")" || die "best_safe failed"
  [ -n "$FINAL" ] || die "no build with every rho <= 1; do not ship"
  log "NOT converged in $MAXL levels; shipping the nearest safe build: ${FINAL##*/}"
fi
FF="$(for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" = "$FINAL" ] && echo "${P##*:}"; done | head -1)"
[ -n "$FF" ] || die "final ${FINAL} has no factor string in the point list"
# Check the chosen artifact before deleting the others, so a meta mismatch still leaves them on disk.
uv run python scripts/def_day_helpers.py meta_matches "$FINAL.h5ad.meta.json" "$FF" "$SQ" \
  || die "final meta does not match factors $FF q $SQ (the h5ad is not the build this run measured)"
uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \
  || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"
for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f "${P%%.ntclevel.json:*}.h5ad" "${P%%.ntclevel.json:*}.h5ad.meta.json"; done
ln -f "$FINAL.h5ad" "$W/${N}_final.h5ad" || die "ln of the final h5ad failed"
ln -f "$FINAL.h5ad.meta.json" "$W/${N}_final.h5ad.meta.json" || die "ln of the final meta failed"
printf 'VCC_BULK_LOAD=%s\nVCC_BULK_SHARPEN_Q=%s\nVCC_BULK_CTX_SCALE=%s\nVCC_SAMELINE=%s\n' "$VCC_BULK_LOAD" "$SQ" "$FF" "${SL:-}" > "$W/def_day.final.env"
[ -n "${VCC_DE_PRIOR_GATE:-}" ] && printf 'VCC_DE_PRIOR_GATE=%s\n' "$VCC_DE_PRIOR_GATE" >> "$W/def_day.final.env"
[ -n "${VCC_DE_MAG:-}" ] && printf 'VCC_DE_MAG=%s\n' "$VCC_DE_MAG" >> "$W/def_day.final.env"
[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf 'VCC_DE_MDE_SIGN=%s\n' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"
log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"

# 5b. the rule-10 diff (481 F5), BEFORE the EXIT trap restores the builder: def_postbuild.py reads the day's operator
#     pin out of the builder source, and the allowed diffs are exactly what this run changed on purpose.
INC="${DEF_DAY_INCUMBENT_META:-submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json}"
ALLOW="whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale"
[ -n "$SL" ] && ALLOW="$ALLOW,sameline"
[ -n "${VCC_DE_PRIOR_GATE:-}" ] && ALLOW="$ALLOW,de_prior_gate"
[ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"
[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"
[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"
# 499 F2: 496's emitter pin moves config.de_sign_emitter_sha256 on every build, knob unset included. Allow that one key
# only when the builder source pin differs from the incumbent's recorded sha (today both are edb5b8dc, so nothing is
# added). A build whose recorded sha moved while the source pin still matches the incumbent still fails the gate.
_inc_em="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$INC" 2>/dev/null || true)"
_src_em="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")' scripts/make_submission_pm_cd4.py 2>/dev/null || true)"
[ -n "$_src_em" ] && [ -n "$_inc_em" ] && [ "$_src_em" != "$_inc_em" ] && ALLOW="$ALLOW,de_sign_emitter_sha256"
uv run python scripts/def_postbuild.py --meta "$W/${N}_final.h5ad.meta.json" --bundle "$B" --cd4-root "$W/cd4_compact" \
    --incumbent "$INC" --allow-config-diff "$ALLOW" --out "$W/postbuild.json" >> "$LOG" 2>&1 \
  && log "postbuild PASS vs ${INC##*/} (allowed: $ALLOW)" || die "postbuild FAILED vs ${INC##*/}: see $W/postbuild.json"

# 6. audits. Write the zero-gene audit to a file so its exit code is the audit's, not head's.
uv run python experiments/428-zero-gene-audit/audit428.py "$W/${N}_final.h5ad" --bundle "$B" > "$W/audit428.txt" 2>&1
rc=$?
summary="$(head -1 "$W/audit428.txt" 2>/dev/null || true)"
log "audit: $summary"
if [ "$rc" -ne 0 ]; then
  die "zero-gene audit failed rc=$rc (see $W/audit428.txt)"
fi
case "$summary" in
  *"ZERO block total (target excluded): 0 in "*) ;;
  *) die "zero-gene audit is not clean: $summary" ;;
esac
python3 - "$W/${N}_final.h5ad.meta.json" <<'PY' | tee -a "$LOG"
import json, sys
c = json.load(open(sys.argv[1]))["config"]
print("meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| pt", c.get("bulk_sharpen_pt"), "| sameline", json.dumps(c.get("sameline"))[:300])
if "de_prior_gate" in c:
    print("meta: de_prior_gate", json.dumps(c["de_prior_gate"])[:400])
PY
if ! "${DRV[@]}" --name "${N}_final" --stages prep > "$W/driver-prep.log" 2>&1; then
  die "prep dry run failed (see $W/driver-prep.log): $(tail -3 "$W/driver-prep.log")"
fi
log "prep dry run: $(grep -iE 'dry run|cells|verified' "$W/driver-prep.log" | head -3 | tr '\n' ' ')"
if [ "${DEF_DAY_KEEP_PINS:-0}" = "1" ]; then
  log "DONE. Next, by hand: vcc prep (real) on $W/${N}_final.h5ad, commit the pins (builder left pinned), vcc submit (-m joke, -d empty)."
else
  log "DONE. Next, by hand: vcc prep (real) on $W/${N}_final.h5ad, vcc submit (-m joke, -d empty). The exit trap restores the builder; re-run the pin stages if you need that diff to commit."
fi
