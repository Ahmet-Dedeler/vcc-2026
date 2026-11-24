#!/usr/bin/env bash
# 692 — R0 readiness adapter (2026-10-15, session [17adab], D/E/F readiness). NEVER submits.
#
# Turns a finished def_day.sh run (691's AUTO_R0, or NEXT.txt by hand) into the evidence the coordinator's slot waiter
# (experiments/670-coord-1015/arm1015.sh) gates on, so R0 goes out through the same tested submit path as every slot
# (536 retry + resume, novelty record, capture, NOGO file):
#   run log   <OUT>/<stem>.run.log   : "[<stem>] PASS …", "[<stem>] CLIFF OK …", "novelty rc=0", "package sha256: <sha>"
#   ready log <OUT>/<stem>.ready.log : "READY <stem> <sha>"  or  "NOT READY <reason>"
#
#   bash experiments/692-def-r0/r0ready692.sh <bundle> <r0 work dir> <stem> [CAP yyyymmddHHMM]
# Env: PKGDIR692 (default submissions; the .vcc goes to <PKGDIR>/<stem>.vcc), OUT692 (default experiments/692-def-r0/work),
#      COMMIT692 (default 1: commit the round's small records with an explicit pathspec once READY; 0 in rehearsals).
#
# Steps, each fail-closed (a failure writes NOT READY and END):
#   1. wait until <work>/def_day.log carries def_day's trap line "builder after restore:" AND no bash process runs
#      def_day.sh for this work dir;
#   2. def_day.log has "postbuild PASS", "DONE.", no "STOP:", and the restored builder equals the pre-pin sha;
#   3. 630dd's "uprn FAILSAFE_TOTAL <n>" line, when present, reads 0 (KA630's rule);
#   4. the final build ("final: <tag> factors …") is not past its cliff (def_day_helpers past_cliff == 1, safe);
#   5. the real `vcc prep` under the gates (no builder, no other vcc prep, free+inactive >= 26 GiB), rc 0, package written;
#   6. check_submission_novelty.py check rc 0;
#   7. provenance (warning only): 690c's builder.pinned.py sha == the final meta's hashes.implementation[builder].
# Terminal line, always: "<ts> END" in the run log.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 1
B="${1:?bundle}"; W="${2:?r0 work dir}"; N="${3:?stem}"; CAP="${4:-$(date -u -v+12H +%Y%m%d%H%M)}"
OUT="${OUT692:-experiments/692-def-r0/work}"; PKGDIR="${PKGDIR692:-submissions}"; COMMIT="${COMMIT692:-1}"
mkdir -p "$OUT" "$PKGDIR"
# Run from a frozen per-launch copy (an edit of this file never corrupts a live run; the name keeps "r0ready692.sh").
if [ -z "${R0READY692_FROZEN:-}" ]; then
  for f in "$OUT"/.frozen-*-r0ready692.sh; do [ -e "$f" ] || continue; p="${f##*/.frozen-}"; p="${p%%-*}"; kill -0 "$p" 2>/dev/null || rm -f "$f"; done
  FZ="$OUT/.frozen-$$-r0ready692.sh"; cp "$REPO/experiments/692-def-r0/r0ready692.sh" "$FZ" || exit 1
  R0READY692_FROZEN=1 exec bash "$FZ" "$@"
fi
RL="$OUT/$N.run.log"; RDY="$OUT/$N.ready.log"; PKG="$PKGDIR/$N.vcc"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$RL"; }
# arm1015.sh greps the run log for "^\[<TAG>\] PASS" and "^\[<TAG>\] CLIFF OK" ANCHORED at line start (run_arm_generic's raw
# lines), so those two evidence lines are also written raw, without the timestamp (found 10-15 while wiring R-ABC)
raw(){ printf '%s\n' "$*" >> "$RL"; }
notready(){ log "NOT READY: $*"; printf '%s NOT READY %s %s\n' "$(date -u +%FT%TZ)" "$N" "$*" >> "$RDY"; log "END"; exit 1; }
now(){ date -u +%Y%m%d%H%M; }
DL="$W/def_day.log"
# 720 R1 (same guard as 691): never two adapters on one stem
PF="$OUT/$N.r0ready692.pid"
if [ -s "$PF" ] && [ "$(cat "$PF")" != "$$" ] && ps -p "$(cat "$PF")" -o args= 2>/dev/null | /usr/bin/grep -q "r0ready692.sh"; then
  log "STOP: another r0ready692 (pid $(cat "$PF")) is alive for $N; not starting a second one"; exit 1
fi
echo $$ > "$PF"
log "r0ready692 armed (pid $$): bundle $B work $W stem $N cap $CAP pkg $PKG commit $COMMIT"
# RESUME after a restart: a package already declared READY (same sha) is done; a package without READY is an
# interrupted prep and is moved aside.
if [ -s "$PKG" ]; then
  _last="$(/usr/bin/grep -E "READY $N [0-9a-f]{64}" "$RDY" 2>/dev/null | tail -1 | awk '{print $NF}')"
  if [ -n "$_last" ] && [ "$(shasum -a 256 "$PKG" | cut -d' ' -f1)" = "$_last" ]; then log "RESUME: $PKG already READY (${_last:0:16}); nothing to do"; log "END"; exit 0; fi
  mv "$PKG" "$PKG.aside-$(date -u +%H%M%S)"; log "RESUME: $PKG existed without a READY line (interrupted prep): moved aside"
fi

# 1. def_day finished
running(){ ps -axo comm=,args= | awk -v w="$W" '$1 ~ /bash$/ && /def_day\.sh/ && index($0, w) > 0' | wc -l | tr -d ' '; }
until /usr/bin/grep -q "builder after restore:" "$DL" 2>/dev/null && [ "$(running)" = 0 ]; do
  [ "$(now)" -ge "$CAP" ] && notready "def_day did not finish by CAP $CAP"
  sleep 60
done
log "def_day finished: $(/usr/bin/grep -E 'DONE|STOP' "$DL" | tail -2 | tr '\n' ' ' | cut -c1-300)"
# 2. its own verdicts
/usr/bin/grep -q "STOP:" "$DL" && notready "def_day.log has a STOP line: $(/usr/bin/grep 'STOP:' "$DL" | tail -1 | cut -c1-200)"
/usr/bin/grep -q "postbuild PASS" "$DL" || notready "no 'postbuild PASS' in def_day.log"
/usr/bin/grep -q "DONE\. Next, by hand" "$DL" || notready "no DONE line in def_day.log"
RB="$(/usr/bin/grep 'builder after restore:' "$DL" | tail -1)"
echo "$RB" | /usr/bin/grep -qE 'restore: ([0-9a-f]{16}) \(want \1\)' || notready "builder not restored to its pre-pin sha: $RB"
log "[$N] PASS postbuild PASS + DONE + $RB"; raw "[$N] PASS postbuild PASS + DONE + builder restored"
# 3. 630dd's fail-safe count, when the UPRN pass ran
FS="$(/usr/bin/grep -oE 'uprn FAILSAFE_TOTAL -?[0-9]+' "$DL" | tail -1 | awk '{print $3}')"
if [ -n "$FS" ]; then [ "$FS" = 0 ] || notready "630dd FAILSAFE_TOTAL $FS (KA630 needs 0)"; log "630dd FAILSAFE_TOTAL 0"; fi
# 3b. 690h's DNMAG delivery report, when the pass ran: logged, never blocking (the coordinator's rule: recorded, not rebuilt)
DNL="$(/usr/bin/grep -oE 'dnmag LOW_CONTEXTS .*' "$DL" | tail -1)"
[ -n "$DNL" ] && log "690h: $DNL (report only)"
# 4. the cliff
FT="$(/usr/bin/grep -oE 'final: [^ ]+ factors' "$DL" | tail -1 | awk '{print $2}')"
[ -n "$FT" ] || notready "no 'final: <tag>' line in def_day.log"
NJ="$W/$FT.ntclevel.json"; [ -s "$NJ" ] || notready "missing $NJ"
uv run python scripts/def_day_helpers.py past_cliff "$NJ" >> "$RL" 2>&1; PC=$?
[ "$PC" = 1 ] || notready "past_cliff rc $PC on $FT (1 = safe)"
RHO="$(uv run python scripts/def_day_helpers.py rho "$NJ" 2>/dev/null)"; log "[$N] CLIFF OK final $FT: $RHO"; raw "[$N] CLIFF OK final $FT: $RHO"
H5="$W/${N}_final.h5ad"; META="$H5.meta.json"
[ -s "$H5" ] || notready "missing $H5"
# 5. the real prep, gated like run_arm_generic's
gate(){
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99)" = 0 ] || return 1
  [ "$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')" = 0 ] || return 1
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)" -ge 26 ] || return 1
}
until gate; do [ "$(now)" -ge "$CAP" ] && notready "prep gate not open by CAP"; sleep 60; done
[ -e "$PKG" ] && notready "$PKG exists (never overwrite a package)"
CTXS="$(python3 -c 'import json,sys; print(",".join(json.load(open(sys.argv[1]))["contexts"]))' "$B/manifest.json")"
/usr/bin/grep -q "UPGRADE (once" experiments/691-release-watcher/work/watch691.log 2>/dev/null && \
  log "WARNING: the vcc CLI was upgraded on the day ($(vcc --version 2>&1 | head -1)): re-run DEF-RUNBOOK §6.7's prep KA before the submit (720 R4)"
t0=$(date +%s)
/usr/bin/time -l vcc prep "$H5" -g "$B/gene_names.csv" --perts "$B/pert_counts.csv" --contexts "$CTXS" -o "$PKG" > "$OUT/$N.prep.log" 2>&1
PRC=$?; log "prep rc=$PRC in $(( $(date +%s) - t0 )) s; $(/usr/bin/grep 'maximum resident' "$OUT/$N.prep.log" | tr -s ' '); package $(stat -f %z "$PKG" 2>/dev/null) bytes"
[ "$PRC" = 0 ] && [ -s "$PKG" ] || notready "vcc prep failed rc $PRC (see $OUT/$N.prep.log)"
# 6. novelty. The meta goes next to the package FIRST: the checker identifies the prediction from <stem>.h5ad.meta.json
#    (output.sha256) and refuses an unidentifiable package (found 10-15 while designing 690's R-KA; run_arm_generic's order)
cp "$META" "$PKGDIR/$N.h5ad.meta.json"
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$RL" 2>&1; NRC=$?
log "novelty rc=$NRC"
[ "$NRC" = 0 ] || notready "novelty check refused"
SHA="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; log "package sha256: $SHA"
# 7. provenance (warning only)
if [ -s "$W/builder.pinned.py" ]; then
  PB="$(shasum -a 256 "$W/builder.pinned.py" | cut -d' ' -f1)"
  MB="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["hashes"]["implementation"].get("scripts/make_submission_pm_cd4.py",""))' "$META" 2>/dev/null)"
  [ "$PB" = "$MB" ] && log "provenance: builder.pinned.py == the meta's builder sha ${PB:0:16}" \
    || log "WARNING provenance: builder.pinned.py ${PB:0:16} != meta builder ${MB:0:16}"
else
  log "WARNING provenance: no builder.pinned.py (690c not applied?)"
fi
printf '%s READY %s %s\n' "$(date -u +%FT%TZ)" "$N" "$SHA" >> "$RDY"
log "READY $N $SHA"
if [ "$COMMIT" = 1 ]; then
  P=()
  for f in def_day.log postbuild.json identity.json census.json preflight.json pin-artifacts.json pin-lodo.json def_day.final.env \
           builder.pinned.py "${N}_final.h5ad.meta.json" timings.json sameline.npz.json whitener_union_pruned_def.json "$FT.ntclevel.json"; do
    [ -s "$W/$f" ] && P+=("$W/$f")
  done
  P+=("$PKGDIR/$N.h5ad.meta.json" "$RL" "$RDY" "$OUT/$N.prep.log")
  git add -f "${P[@]}" && git commit -q -m "692: R0 records ($N, package $SHA): def_day evidence, pinned builder, meta (r0ready692)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -- "${P[@]}" \
    && log "committed the round's records: $(git log -1 --format=%h)" || log "WARNING commit of the records failed (commit by hand)"
fi
log "END"
