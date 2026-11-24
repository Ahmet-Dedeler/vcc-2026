#!/usr/bin/env bash
# 690 — the October dress rehearsal chain (2026-10-15, session [17adab]): R-NEW, then R-KA (v3; R-SAME retired).
#   R-NEW  def_day.sh on work/bundles/new  (foreign controls HCT116 / H1 / HEK293T, ef-bundle-r's NEW 300-target panel)
#   R-KA   the whole Oct 22 chain (691 watcher AUTO_R0 on a stub release of 180's def-bundle -> def_day -> 692) with
#          524dr2's own decisions snapshot: PASS iff it rebuilds 524dr2's BOARD-SCORED prediction (entry dVHLws6K, 0.2730889)
# R-NEW uses the resolved October recipe (`october_recipe.py command --allow-pending --json`, PENDING treated off, recorded
# per run) plus DEF_DAY_ALLOW_DEPTH=1 (Orion's controls sit at x0.80-0.84 of A/B/C's depth; real D/E/F ship at ~20k).
# After R-NEW's DONE: the REAL `vcc prep` (no upload, never `vcc submit`) to time the hand step, then the .vcc is deleted.
#
# Waits (AGENTS 6.2): the coordinator's dn671 PROCESS exits (then its log tail is recorded), then before EACH rehearsal:
#   no HOLD file ($W/HOLD: the coordinator touches it to put a slot-bound build first), the builder on disk == HEAD
#   (def_day refuses otherwise), no builder / operator fit running (the driver's own other_builders()), no other
#   day-path chain (def_day.sh / run_arm_generic.sh / run_ka_527.sh under bash), free+inactive >= 20 GiB and memory_pressure free >= 40 %, disk >= 35 GiB.
# Waiting to START stops at CAP. Never submits, never edits the builder (def_day's exit trap restores its pins).
# Terminal line, always: "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work; LOG=$W/run690.log; HOLD=$W/HOLD
# Run from a frozen per-launch copy, so an edit of this file never corrupts a live run (bash reads scripts by offset)
if [ -z "${RUN690_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-run690.sh"; cp "$REPO/$E/run690.sh" "$FZ" || exit 1
  RUN690_FROZEN=1 exec bash "$FZ" "$@"
fi
CAP="${CAP690:-202610190000}"   # v5 (10-15 21:30Z): the rehearsals now run after build705 (10-16/17); done before the 10-20 freeze
DNPID="$(cat experiments/671-dnmag-arm/work/dn671.pid 2>/dev/null || echo 0)"
OPS=experiments/453-abc-dress-rehearsal/work/run
OPSHA=66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"' EXIT
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
n_builders(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99; }
disk_gib(){ df -g /System/Volumes/Data | awk 'NR==2{print $4}'; }
chains(){ ps -axo comm=,args= | awk '$1 ~ /bash$/ && /(def_day\.sh|run_arm_generic\.sh|run_ka_527\.sh)/' | grep -v "run690" | wc -l | tr -d ' '; }
gate(){
  [ -e "$HOLD" ] && { echo "HOLD file ($(head -c 200 "$HOLD"))"; return 1; }
  local hb db; hb="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"; db="$(sha scripts/make_submission_pm_cd4.py)"
  [ "$hb" = "$db" ] || { echo "builder on disk ${db:0:16} != HEAD ${hb:0:16}"; return 1; }
  [ "$(n_builders)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  [ "$(chains)" = 0 ] || { echo "a day-path chain is running"; return 1; }
  local f d; f="$(fi_gib)"; d="$(disk_gib)"
  # 10-16: the AGENTS §9 rule (free+inactive >= 20 GiB AND memory_pressure free >= 40 %; no swap clause: since the 10-15
  # reboot there is no swap file). 26 never opened with ~7 GB held by other apps; R-ABC's level builds peaked at 20.2 GiB RSS.
  local mp; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$f" -ge "${FI_MIN690:-20}" ] || { echo "free+inactive ${f} GiB < ${FI_MIN690:-20}"; return 1; }
  [ "${mp:-0}" -ge "${MP_MIN690:-40}" ] || { echo "memory_pressure free ${mp:-?} % < ${MP_MIN690:-40}"; return 1; }
  [ "$d" -ge 35 ] || { echo "disk ${d} GiB < 35"; return 1; }
  echo "builder ${db:0:16} == HEAD, no builder, no chain, free+inactive ${f} GiB, disk ${d} GiB"
}
wait_gate(){  # $1 = label, $2 = latest start YYYYmmddHHMM (optional): rc 2 = deferred (the slot-bound R-ABC goes first)
  local last="" r
  until r="$(gate)"; do
    [ "$r" != "$last" ] && log "$1 waiting: $r"; last="$r"
    [ -n "${2:-}" ] && [ "$(date -u +%Y%m%d%H%M)" -ge "$2" ] && { log "$1 DEFERRED: not started by $2 ($r)"; return 2; }
    capped && { log "STOP: $1 not started by CAP $CAP ($r)"; return 1; }
    sleep 60
  done
  [ -n "${2:-}" ] && [ "$(date -u +%Y%m%d%H%M)" -ge "$2" ] && { log "$1 DEFERRED: its gate opened after $2"; return 2; }
  log "$1 gate PASS: $r"
}
rehearse(){  # $1 = label, $2 = bundle, $3 = work dir, $4 = name, $5 = pre-place the 453 operator (1/0)
  local L="$1" B="$2" WD="$3" N="$4" PRE="$5"
  [ -e "$WD" ] && { mv "$WD" "$WD.aside-$(date -u +%H%M%S)"; log "$L: moved an old $WD aside"; }
  mkdir -p "$WD"
  if [ "$PRE" = 1 ]; then
    [ "$(sha $OPS/whitener_union_pruned_def.npz)" = "$OPSHA" ] || { log "$L STOP: 453's operator is not 66610be7"; return 1; }
    cp -c "$OPS/whitener_union_pruned_def.npz" "$OPS/whitener_union_pruned_def.json" "$WD/" || { log "$L STOP: operator copy failed"; return 1; }
    log "$L: pre-placed 453's operator 66610be7 (+ json) in $WD (the same-panel branch)"
  fi
  local RJ="$W/recipe-$L.json"
  uv run python scripts/october_recipe.py command --bundle "$B" --work "$WD" --name "$N" --allow-pending --json > "$RJ" 2>> "$LOG" \
    || { log "$L STOP: october_recipe command failed"; return 1; }
  local ENVS ARGS
  ENVS="$(python3 -c 'import json,shlex,sys; print(" ".join(f"{k}={shlex.quote(v)}" for k,v in json.load(open(sys.argv[1]))["env"].items()))' "$RJ")"
  ARGS="$(python3 -c 'import json,shlex,sys; print(" ".join(shlex.quote(a) for a in json.load(open(sys.argv[1]))["args"]))' "$RJ")"
  log "$L: env DEF_DAY_ALLOW_DEPTH=1 $ENVS; args $ARGS; pending treated off: $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["pending_treated_off"])' "$RJ")"
  local t0; t0=$(date +%s)
  eval "env DEF_DAY_ALLOW_DEPTH=1 $ENVS bash scripts/def_day.sh $ARGS" > "$WD.nohup" 2>&1
  local rc=$?
  log "$L: def_day rc $rc in $(( ($(date +%s) - t0) / 60 )) min; $(grep -E 'DONE|STOP' "$WD/def_day.log" 2>/dev/null | tail -2 | tr '\n' ' ' | cut -c1-400)"
  log "$L: builder after: $(sha scripts/make_submission_pm_cd4.py | cut -c1-16) (HEAD $(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -c1-16))"
  return $rc
}

log "run690 armed (pid $$): wait for dn671 pid $DNPID to exit; CAP $CAP"
echo $$ > "$W/run690.pid"
# v5: a pid alone is not identity (after a reboot 20541 can be any process): the pid must still be dn671.sh under bash
while [ "$DNPID" != 0 ] && ps -o comm=,args= -p "$DNPID" 2>/dev/null | awk '$1 ~ /bash$/ && /dn671\.sh/' | /usr/bin/grep -q .; do
  capped && { log "STOP: dn671 still running at CAP $CAP"; exit 3; }
  sleep 60
done
log "dn671 exited; its tail: $(tail -3 experiments/671-dnmag-arm/work/dn671.log 2>/dev/null | tr '\n' ' ' | cut -c1-400)"
# 694 (10-15 ~05:10Z): ka693b (KAu2 then DNv2) must EXIT first. Its gate sees only the python builder / vcc prep, not
# def_day.sh's setup stage, so once the HOLD lifts, R-NEW and DNv2 could both start inside one 60 s poll (two builds).
ka693b_live(){ ps -axo comm=,args= | awk '$1 ~ /bash$/ && /experiments\/693-dnmag-speed\/ka693b\.sh/' | /usr/bin/grep -q .; }
if ka693b_live; then
  log "waiting for ka693b (KAu2 + DNv2) to exit"
  while ka693b_live; do capped && { log "STOP: ka693b still running at CAP $CAP"; exit 3; }; sleep 60; done
fi
log "ka693b not running; its verdicts: $(/usr/bin/grep -E 'KAU2 (PASS|FAIL)|DNV2 (PASS|FAIL)|STOP' experiments/693-dnmag-speed/work/ka693b.log 2>/dev/null | tr '\n' ' ' | cut -c1-300)"

# ---- R-NEW
# v4/v5: optional latest START times (RNEW_LATEST690 / RKA_LATEST690, YYYYmmddHHMM; default none): a rehearsal whose gate
# opens later is DEFERRED (logged, exit 0). ONLY690=rka runs R-KA alone. Since the 10-15 reboot the rehearsals run AFTER the
# slot builds (R-ABC, K-UP, R7a, build705), started by after690.sh, so no deadline is set by default. ONLY690=rnew runs R-NEW alone (690j).
RNEW_LATEST="${RNEW_LATEST690:-}"; RKA_LATEST="${RKA_LATEST690:-}"; RCN=skipped
if [ "${ONLY690:-}" != rka ]; then
wait_gate R-NEW "$RNEW_LATEST"; rc=$?
[ "$rc" = 2 ] && { log "SUMMARY R-NEW deferred, R-KA deferred (re-run both after R-ABC)"; exit 0; }
[ "$rc" = 0 ] || exit 3
rehearse R-NEW "$W/bundles/new" "$W/dr_new" dr690n 0; RCN=$?
if [ "$RCN" = 0 ] && [ -s "$W/dr_new/dr690n_final.h5ad" ]; then
  if r="$(gate)"; then
    t0=$(date +%s)
    /usr/bin/time -l vcc prep "$W/dr_new/dr690n_final.h5ad" -g "$W/bundles/new/gene_names.csv" --perts "$W/bundles/new/pert_counts.csv" \
      --contexts D,E,F -o "$W/dr_new/dr690n.vcc" > "$W/dr_new/prep_real.log" 2>&1
    log "R-NEW: real vcc prep rc $? in $(( $(date +%s) - t0 )) s; .vcc $(stat -f %z "$W/dr_new/dr690n.vcc" 2>/dev/null) bytes; $(grep -E 'maximum resident' "$W/dr_new/prep_real.log" | tr -s ' ')"
    rm -f "$W/dr_new/dr690n.vcc"
  else
    log "R-NEW: real prep skipped ($r)"
  fi
fi
# The final h5ad (~12 GB) has no truth to score against: keep its meta, logs and audits, free the disk.
[ "$RCN" = 0 ] && rm -f "$W/dr_new/dr690n_final.h5ad" "$W"/dr_new/dr690n-lvl*.h5ad
fi

# ---- R-KA (v3, 2026-10-15): the BOARD-VALIDATED known answer of the whole Oct 22 chain, replacing R-SAME. 691's watcher
#      (AUTO_R0=1) on a stub release of 180's def-bundle (A/B/C relabelled D/E/F: the SAME-PANEL branch, 453's operator
#      pre-placed), the recipe resolved from 524dr2's OWN decisions snapshot (OCTOBER_DECISIONS), today's patched builder
#      (631 + 690a/c/g, all inert unset); then 692 (real prep + novelty). PASS iff the novelty check recognises the
#      prediction as 524dr2's board-scored one: "ALREADY BEEN SCORED", entry dVHLws6KdMxudq0Yn5ap (board 0.2730889).
#      Any other outcome is a bug to find before Oct 22. No commit, no submit.
# 690j: ONLY690=rnew runs R-NEW alone (its rerun after def_lodo_key's RSS-guard fix; R-KA already ran on this tree)
[ "${ONLY690:-}" = rnew ] && { log "SUMMARY R-NEW rc $RCN (ONLY690=rnew: R-KA not run)"; exit 0; }
wait_gate R-KA "$RKA_LATEST"; rc=$?
[ "$rc" = 2 ] && { log "SUMMARY R-NEW rc $RCN, R-KA deferred (re-run with ONLY690=rka after R-ABC)"; exit 0; }
[ "$rc" = 0 ] || exit 3
SNAP=experiments/658-prereg-1007-builds/work/dr2/october-decisions.at-resolve.json
log "R-KA: decisions snapshot $SNAP sha256 $(sha "$SNAP" | cut -c1-16)"
ZKA="$W/zips/final_defbundle.zip"
if [ ! -s "$ZKA" ]; then
  mkdir -p "$W/zips" "$W/zipka/vcc_2026_final"; cp -c experiments/180-def-readiness/work/def-bundle/* "$W/zipka/vcc_2026_final/"
  (cd "$W/zipka" && zip -q -0 -r "$REPO/$ZKA" vcc_2026_final) || { log "R-KA STOP: zip failed"; exit 1; }
  rm -rf "$W/zipka"
fi
WW="$W/rka_watch"
for p in "$WW" "$W/def_rka" "$W/dl_rka" "$W/dr_ka" "$W/r0ready_ka"; do [ -e "$p" ] && mv "$p" "$p.aside-$(date -u +%H%M%S)"; done
mkdir -p "$WW"
t0=$(date +%s)
STUB_STATE="$REPO/$WW.state" STUB_ZIP="$REPO/$ZKA" VCC691="$REPO/experiments/691-release-watcher/vcc_stub691.sh" \
  W691="$WW" DL691="$W/dl_rka" BUNDLE691="$REPO/$W/def_rka" R0W691="$W/dr_ka" AUTO_R0=1 R0_ALLOW_PENDING691=1 \
  OCTOBER_DECISIONS="$REPO/$SNAP" START691="$(date -u +%Y%m%d%H%M)" CAP691="$(date -u -v+4H +%Y%m%d%H%M)" EVERY691=2 \
  bash experiments/691-release-watcher/watch691.sh; RCK=$?
log "R-KA: watch691 rc $RCK in $(( ($(date +%s) - t0) / 60 )) min; $(grep -E 'branch|AUTO_R0' "$WW/watch691.log" | tail -3 | tr '\n' ' ' | cut -c1-500)"
if [ "$RCK" = 0 ] && [ -s "$W/dr_ka/r0def_final.h5ad" ]; then
  OUT692="$W/r0ready_ka" PKGDIR692="$W/dr_ka" COMMIT692=0 bash experiments/692-def-r0/r0ready692.sh "$W/def_rka" "$W/dr_ka" r0def \
    "$(date -u -v+3H +%Y%m%d%H%M)"
  RR="$W/r0ready_ka/r0def.run.log"
  # 690k (10-16): judge the def_day OUTPUT, not the package. The board scored relabel453's A/B/C clone of 524dr2's def_day output,
  # and novelty hashes `context`, so a D/E/F package can never be "ALREADY BEEN SCORED" (the 10-16 09:12Z false FAIL). PASS iff
  # this run's output.sha256 == 524dr2's relabel453.source_output_sha256 (== its def_day output) and the content-aware meta
  # diff vs 524dr2's final is empty.
  KAO="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$W/dr_ka/r0def_final.h5ad.meta.json" 2>/dev/null)"
  KAS="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["relabel453"]["source_output_sha256"])' submissions/524dr2-october-recipe-abc.h5ad.meta.json 2>/dev/null)"
  KAMD="$(uv run python scripts/def_metadiff.py "$W/dr_ka/r0def_final.h5ad.meta.json" "$W/dr_ka" experiments/658-prereg-1007-builds/work/dr2/run/524dr2_final.h5ad.meta.json experiments/658-prereg-1007-builds/work/dr2/run --expect '' 2>&1 | tail -1)"
  if [ -n "$KAO" ] && [ "$KAO" = "$KAS" ] && [ "${KAMD#\[metadiff\] PASS}" != "$KAMD" ]; then
    log "R-KA PASS: the Oct 22 chain rebuilt 524dr2's def_day output byte for byte (output.sha256 ${KAO:0:16} == relabel453.source_output_sha256 of the board-scored package, entry dVHLws6KdMxudq0Yn5ap, 0.2730889); $KAMD"
    # 690i: the tree this proof covers (708 quotes it before any final build uses DEF_DAY_H1_STACK / the widened CD4 grid)
    log "R-KA tree: def_day.sh $(sha scripts/def_day.sh) helpers $(sha scripts/def_day_helpers.py) builder $(sha scripts/make_submission_pm_cd4.py | cut -c1-16) at $(git rev-parse --short HEAD); 690i $(/usr/bin/grep -q "690i: the doses come from 438's grid" scripts/def_day.sh && echo APPLIED || echo 'NOT applied')"
    RCS=0
  else
    log "R-KA FAIL: the def_day output is not 524dr2's: output.sha256 ${KAO:-?} vs ${KAS:-?}; $KAMD"
    log "R-KA tree: def_day.sh $(sha scripts/def_day.sh) at $(git rev-parse --short HEAD); 690i $(/usr/bin/grep -q "690i: the doses come from 438's grid" scripts/def_day.sh && echo APPLIED || echo 'NOT applied')"
    log "R-KA: built output sha $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$W/dr_ka/r0def_final.h5ad.meta.json" 2>/dev/null) vs 524dr2 0c394160ceaba004 (ledger content 50ec8ff127acfd99)"
    RCS=2
  fi
  rm -f "$W/dr_ka/r0def.vcc"
else
  log "R-KA FAIL: the watcher did not build R0 (rc $RCK)"; RCS=1
fi
[ -s "$W/dr_ka/r0def_final.h5ad" ] && rm -f "$W/dr_ka/r0def_final.h5ad" "$W"/dr_ka/r0def-lvl*.h5ad
rm -rf "$W/dl_rka"
log "SUMMARY R-NEW rc $RCN, R-KA rc $RCS"
