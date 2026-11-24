#!/usr/bin/env bash
# 694 re-gate (10-16, readiness [4592ee]). R-ABC's build ended "NOT READY gate694 failed" (rabc694.log, 03:51:36Z) on a gate
# bug, not a build fault: the strict meta diff vs 524dr2 counted bulk_load_sha256 and sameline (sha256) although both npz
# files equal 524dr2's in every array; only the work-dir path inside a string field differs. gate694 v2 accepts those two
# keys only after that content check (test694 G6-G9). This script re-gates the FINISHED build, never rebuilds:
#   1. the package's sha256 must equal the one rabc694.log logged at build time
#   2. novelty is re-run on it; "novelty rc=<n>" logged
#   3. F1 (bulk identity vs 524dr2, IDENTICAL) is quoted from rabc694.log; nothing it read has changed since
#   4. gate694 (v2) --write; its "[rabc] ..." lines go RAW into the new log (arm1015's anchored greps)
# New ready/run log: work/rabc694.r2.log, terminal "READY rabc <sha>" | "NOT READY <reason>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/694-rabc; W=$X/work; L1=$W/rabc694.log; LOG=$W/rabc694.r2.log; STEM=submissions/694rabc-october-recipe-abc
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0; finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc)"' EXIT
[ -e "$LOG" ] && mv "$LOG" "$LOG.aside-$(date -u +%H%M%S)"
log "ready694r2 (pid $$): re-gate of the finished R-ABC build with gate694 v2 ($(shasum -a 256 $X/gate694.py | cut -c1-16)); first run: $(tail -1 "$L1" | cut -c1-120)"
WANT="$(/usr/bin/grep -oE "package $STEM.vcc sha256 [0-9a-f]{64}" "$L1" | tail -1 | awk '{print $4}')"
[ -n "$WANT" ] || finish "NOT READY no package line in $L1"
GOT="$(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
[ "$GOT" = "$WANT" ] || finish "NOT READY package sha $GOT != the build's $WANT"
log "package $STEM.vcc sha256 $GOT"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
F1="$(/usr/bin/grep -E "F1 bulk identity vs 524dr2: IDENTICAL" "$L1" | tail -1)"
[ -n "$F1" ] || finish "NOT READY no F1 IDENTICAL line in $L1"
log "F1 (quoted from rabc694.log): ${F1#* }"
uv run python $X/gate694.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate694 v2 rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate694 v2 failed (see the gate lines above)"
[ -f "$STEM.h5ad" ] && { rm -f "$STEM.h5ad"; log "removed $STEM.h5ad (the package exists, its sha is logged, every gate passed)"; }
finish "READY rabc $GOT"
