#!/usr/bin/env bash
# 690 R-KA2 (2026-10-18, readiness [4592ee]; preregistered in RESULTS.md §6c, c2b702e0): the whole Oct 22 chain (691 watcher AUTO_R0
# on a stub release of 180's def-bundle -> def_day -> 692) with the recipe resolved from the LIVE docs/october-decisions.json
# (OCTOBER_DECISIONS unset, no R0_ALLOW_PENDING691). PASS iff the def_day output.sha256 == R-ABC2's def_day final
# (c0a636a0..., board 0.2739849) AND def_metadiff vs R-ABC2's final is empty AND 692 writes READY AND the decisions file is
# unchanged across the run. This is run690.sh's R-KA section with only those changes; separate dirs (*_ka2, rka2_watch), so the
# real watch691 / r0ready692 / arm1015 R0 never see it. Never submits; the trial .vcc is deleted.
# Gates as run690 (HOLD, builder == HEAD, no builder, no day-path chain, AGENTS §9 memory, disk >= 35 GiB); latest start
# RKA2_LATEST (default 202610181830): after it, DEFERRED (exit 0). Terminal line, always: "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work; LOG=$W/rka2_690.log; HOLD=$W/HOLD
if [ -z "${RKA2_FROZEN:-}" ]; then   # run from a frozen per-launch copy (bash reads scripts by offset)
  FZ="$W/.frozen-$$-rka2_690.sh"; cp "$REPO/$E/rka2_690.sh" "$FZ" || exit 1
  RKA2_FROZEN=1 exec bash "$FZ" "$@"
fi
LATEST="${RKA2_LATEST:-202610181830}"
DEC=docs/october-decisions.json
REF_META=experiments/695-rabc2/work/run/rabc2_final.h5ad.meta.json; REF_RUN=experiments/695-rabc2/work/run
REF_SHA=c0a636a00437120705a8d5d8ba035bebfc84a44e1c31fcdca812fa7f98aae2fb
ZKA="$W/zips/final_defbundle.zip"; ZSHA=6d7de7718e4e1491
OPS=experiments/453-abc-dress-rehearsal/work/run
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"; rm -f "$W/.frozen-$$-rka2_690.sh"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
n_builders(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
chains(){ ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh)/' | wc -l | tr -d ' '; }
gate(){
  [ -e "$HOLD" ] && { echo "HOLD file ($(head -c 200 "$HOLD"))"; return 1; }
  local hb db; hb="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"; db="$(sha scripts/make_submission_pm_cd4.py)"
  [ "$hb" = "$db" ] || { echo "builder on disk ${db:0:16} != HEAD ${hb:0:16}"; return 1; }
  [ "$(n_builders)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  [ "$(chains)" = 0 ] || { echo "a day-path chain is running"; return 1; }
  local f d mp; f="$(fi_gib)"; d="$(disk_gib)"
  mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$f" -ge 20 ] || { echo "free+inactive ${f} GiB < 20"; return 1; }
  [ "${mp:-0}" -ge 40 ] || { echo "memory_pressure free ${mp:-?} % < 40"; return 1; }
  [ "$d" -ge 35 ] || { echo "disk ${d} GiB < 35"; return 1; }
  echo "builder ${db:0:16} == HEAD, no builder, no chain, free+inactive ${f} GiB, memory_pressure free ${mp} %, disk ${d} GiB"
}

echo $$ > "$W/rka2_690.pid"
log "rka2_690 armed (pid $$): latest start $LATEST; ref $REF_META (output ${REF_SHA:0:16}); decisions $DEC $(sha $DEC | cut -c1-16)"
[ -z "${OCTOBER_DECISIONS:-}" ] || log "note: OCTOBER_DECISIONS was set in the caller ($OCTOBER_DECISIONS); it is unset for the chain"
# preflight of the fixed inputs (any mismatch = STOP before building)
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$REF_META" 2>/dev/null)" = "$REF_SHA" ] \
  || { log "STOP: $REF_META does not record output $REF_SHA"; exit 1; }
[ "$(sha "$ZKA" | cut -c1-16)" = "$ZSHA" ] || { log "STOP: stub zip $ZKA is not $ZSHA"; exit 1; }
[ "$(sha $OPS/whitener_union_pruned_def.npz | cut -c1-8)" = 66610be7 ] || { log "STOP: 453's operator is not 66610be7"; exit 1; }
PEND="$(env -u OCTOBER_DECISIONS uv run python scripts/october_recipe.py check 2>&1 | tail -1)"; PRC=$?
log "october_recipe check on the live file: $PEND"
DEC0="$(sha $DEC)"

last=""
until r="$(gate)"; do
  [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
  [ "$(date -u +%Y%m%d%H%M)" -ge "$LATEST" ] && { log "DEFERRED: not started by $LATEST ($r); re-run after R0"; exit 0; }
  sleep 60
done
log "gate PASS: $r"
log "tree: def_day.sh $(sha scripts/def_day.sh | cut -c1-16) helpers $(sha scripts/def_day_helpers.py | cut -c1-16) builder $(sha scripts/make_submission_pm_cd4.py | cut -c1-16) emitter $(sha src/vcc2026/rc_emit_k.py | cut -c1-16) at $(git rev-parse --short HEAD)"
WW="$W/rka2_watch"
# 10-20: also the unpack scratch dir: watch691 refuses a stale "$BUNDLE.unzipped" from an earlier run (R-KA3 stopped on it)
for p in "$WW" "$WW.state" "$WW.state.dl" "$W/def_rka2" "$W/def_rka2.unzipped" "$W/dl_rka2" "$W/dr_ka2" "$W/r0ready_ka2"; do [ -e "$p" ] && mv "$p" "$p.aside-$(date -u +%H%M%S)"; done
mkdir -p "$WW"
t0=$(date +%s)
env -u OCTOBER_DECISIONS -u R0_ALLOW_PENDING691 \
  STUB_STATE="$REPO/$WW.state" STUB_ZIP="$REPO/$ZKA" VCC691="$REPO/experiments/691-release-watcher/vcc_stub691.sh" \
  W691="$WW" DL691="$W/dl_rka2" BUNDLE691="$REPO/$W/def_rka2" R0W691="$W/dr_ka2" AUTO_R0=1 \
  START691="$(date -u +%Y%m%d%H%M)" CAP691="$(date -u -v+2H +%Y%m%d%H%M)" EVERY691=2 \
  bash experiments/691-release-watcher/watch691.sh; RCK=$?
log "watch691 rc $RCK in $(( ($(date +%s) - t0) / 60 )) min; $(grep -E 'branch|AUTO_R0' "$WW/watch691.log" | tail -4 | tr '\n' ' ' | cut -c1-600)"
log "the watcher's resolved recipe env: $(python3 -c 'import json,sys; print(sorted(json.load(open(sys.argv[1]))["env"].items()))' "$WW/r0-recipe.json" 2>/dev/null)"
DEC1="$(sha $DEC)"
[ "$DEC1" = "$DEC0" ] && DECOK=1 || DECOK=0
log "decisions file across the run: ${DEC0:0:16} -> ${DEC1:0:16} ($([ $DECOK = 1 ] && echo unchanged || echo CHANGED))"
RCS=1
if [ "$RCK" = 0 ] && [ -s "$W/dr_ka2/r0def_final.h5ad" ]; then
  OUT692="$W/r0ready_ka2" PKGDIR692="$W/dr_ka2" COMMIT692=0 bash experiments/692-def-r0/r0ready692.sh "$W/def_rka2" "$W/dr_ka2" r0def \
    "$(date -u -v+3H +%Y%m%d%H%M)"
  R692=$?
  RDY="$(tail -1 "$W/r0ready_ka2/r0def.ready.log" 2>/dev/null)"
  log "692 rc $R692: ${RDY:-no ready log}"
  KAO="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$W/dr_ka2/r0def_final.h5ad.meta.json" 2>/dev/null)"
  KAMD="$(uv run python scripts/def_metadiff.py "$W/dr_ka2/r0def_final.h5ad.meta.json" "$W/dr_ka2" "$REF_META" "$REF_RUN" --expect '' \
    --json "$W/rka2_metadiff.json" 2>&1 | tail -1)"
  log "def_day output.sha256 ${KAO:-?} vs R-ABC2 $REF_SHA; $KAMD"
  RDYOK=0; case "$RDY" in *" NOT READY"*) ;; *" READY r0def "[0-9a-f]*) RDYOK=1 ;; esac
  if [ -n "$KAO" ] && [ "$KAO" = "$REF_SHA" ] && [ "${KAMD#\[metadiff\] PASS}" != "$KAMD" ] && [ "$RDYOK" = 1 ] && [ "$DECOK" = 1 ]; then
    log "R-KA2 PASS: the Oct 22 chain on the LIVE october-decisions.json rebuilt R-ABC2's def_day output byte for byte (${KAO:0:16}; board 0.2739849, nnHQfcZHsT3as9EafA3i); metadiff empty; 692 READY; decisions unchanged"
    RCS=0
  else
    log "R-KA2 FAIL: output ${KAO:-?} (want $REF_SHA); metadiff '$KAMD'; 692 '${RDY:-none}'; decisions $([ $DECOK = 1 ] && echo unchanged || echo CHANGED)"
    RCS=2
  fi
  rm -f "$W/dr_ka2/r0def.vcc"
else
  log "R-KA2 FAIL: the watcher did not build R0 (rc $RCK)"
fi
# the final h5ad (~12 GB) is redundant once its sha and meta are logged; keep logs, metas, audits
[ -s "$W/dr_ka2/r0def_final.h5ad" ] && rm -f "$W/dr_ka2/r0def_final.h5ad" "$W"/dr_ka2/r0def-lvl*.h5ad
rm -rf "$W/dl_rka2"
log "SUMMARY R-KA2 rc $RCS"
exit $RCS
