#!/usr/bin/env bash
# def_day.sh — the October D/E/F day in one gated script (2026-09-25, coordinator; generalises 448's run448.sh).
# 485 review patch: fail closed on driver rc, dump rc, audit, and prep; a bad secant step ships the
# nearest safe build instead of aborting; lvlC must sit inside the 1% band; past_cliff crashes are
# not "safe". Apply together with the def_pipeline.py dump-assert skip. The current driver exits 1
# after every successful dump, so this shell's dump rc check fails until that patch is in.
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
#   DEF_DAY_CD4_THRESHOLDS=relaxed (default unset = 047's gates, 275 per condition and 290 union): the only legal
#    non-empty value. A synthetic panel that cannot clear the gates sets this; the day leaves it unset.
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
# Rehearse on a synthetic bundle first (e.g. experiments/238-essential-free-panel/work/ef-bundle-r).
# A rerun needs a fresh work dir: this script deletes a tag's h5ad before rebuilding it, and refuses to
# treat a pre-existing file as the build it just asked for.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
B="${1:?bundle}"; W="${2:?work dir}"; N="${3:?name}"; Q="${4:-1.75}"; SLW="${5:-auto}"; SLF="${6:-pt:noprog}"; MAXL="${7:-4}"
mkdir -p "$W"; LOG="$W/def_day.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
die(){ log "STOP: $*"; exit 1; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "start: bundle $B work $W name $N q $Q sameline w $SLW flags $SLF; builder ${BSHA0:0:16} driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"
grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_SAMELINE (apply 451's patch first)"
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
if [ "${DEF_DAY_CD4_THRESHOLDS:-}" = "relaxed" ]; then
  DRV+=(--cd4-thresholds relaxed)
  log "CD4 extract thresholds relaxed (DEF_DAY_CD4_THRESHOLDS=relaxed); synthetic panels only"
elif [ -n "${DEF_DAY_CD4_THRESHOLDS:-}" ]; then
  die "DEF_DAY_CD4_THRESHOLDS must be relaxed or unset (got ${DEF_DAY_CD4_THRESHOLDS})"
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
    > "$W/driver-setup.log" 2>&1
grep -q "^--- \[repin-panel-support\].*rc=0" "$W/driver-setup.log" || die "setup did not complete (see driver-setup.log)"
log "setup: $(grep -E '^--- \[' "$W/driver-setup.log" | awk '{print $2, $3}' | tr '\n' ' ')"
[ -s "$W/identity.json" ] || die "no identity.json"
log "identity: $(python3 -c "import json; d=json.load(open('$W/identity.json')); print(d['call'], 'unanimous', d['call_unanimous_across_methods'])")"
# Context labels come from the preflight the setup stage just wrote. Hardcoding D,E,F ships a
# mismatch the moment the manifest uses other names: lvl0 has no scale and succeeds, then ntclevel
# looks for context_D.h5ad. isidentifier() rejects anything that is not a bare name.
CTXS="$(python3 -c 'import json,sys
cs=json.load(open(sys.argv[1]))["contexts"]
assert isinstance(cs, list) and len(cs)==3 and all(isinstance(c, str) and c.isidentifier() for c in cs), cs
print(",".join(cs))' "$W/preflight.json")" || die "could not read three context names from $W/preflight.json"
FACT="$(python3 -c 'import sys; print(",".join("%s=1" % c for c in sys.argv[1].split(",")))' "$CTXS")" || die "could not build the unit factor string"
log "contexts $CTXS; unit factors $FACT"

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
  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"
  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \
      --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"
  uv run python scripts/def_pin_panel_artifacts.py --whitener "$W/whitener_union_pruned_def.npz" \
      --whitener-metadata "$W/whitener_union_pruned_def.json" --whitener-id union_pruned_def --cd4-root "$W/cd4_compact" \
      --pert-scale "$LK" --pert-scale-id lodo_lo0.5_hi3_def --apply --out "$W/pin-lodo.json" >> "$LOG" 2>&1 \
      || die "LODO key pin failed"
  LODO=lodo_lo0.5_hi3_def
  log "LODO key ON for the stack dump: $LODO $(shasum -a 256 "$LK" | cut -c1-16)"
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
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$Q"
[ -n "$SL" ] && export VCC_SAMELINE="$SL"
if [ -n "${DEF_DAY_PRIOR_GATE:-}" ]; then
  [ -n "$SL" ] || die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE, but same-line mode is off: the programme never reached the builder"
  export VCC_DE_PRIOR_GATE="$DEF_DAY_PRIOR_GATE"
fi

# 4-5. lvl0, then the re-level loop
POINTS=()
build(){ local tag=$1 fact=$2 h5 rc pc
  if [ "$fact" = "$FACT" ]; then unset VCC_BULK_CTX_SCALE; else export VCC_BULK_CTX_SCALE="$fact"; fi
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
uv run python scripts/def_day_helpers.py meta_matches "$FINAL.h5ad.meta.json" "$FF" "$Q" \
  || die "final meta does not match factors $FF q $Q (the h5ad is not the build this run measured)"
for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f "${P%%.ntclevel.json:*}.h5ad" "${P%%.ntclevel.json:*}.h5ad.meta.json"; done
ln -f "$FINAL.h5ad" "$W/${N}_final.h5ad" || die "ln of the final h5ad failed"
ln -f "$FINAL.h5ad.meta.json" "$W/${N}_final.h5ad.meta.json" || die "ln of the final meta failed"
printf 'VCC_BULK_LOAD=%s\nVCC_BULK_SHARPEN_Q=%s\nVCC_BULK_CTX_SCALE=%s\nVCC_SAMELINE=%s\n' "$VCC_BULK_LOAD" "$Q" "$FF" "${SL:-}" > "$W/def_day.final.env"
[ -n "${VCC_DE_PRIOR_GATE:-}" ] && printf 'VCC_DE_PRIOR_GATE=%s\n' "$VCC_DE_PRIOR_GATE" >> "$W/def_day.final.env"
log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"

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
print("meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| sameline", json.dumps(c.get("sameline"))[:300])
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
