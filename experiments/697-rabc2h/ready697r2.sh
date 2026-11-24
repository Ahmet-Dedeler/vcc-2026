#!/usr/bin/env bash
# 697 re-gate (10-17, readiness [4592ee]). R-ABC2H's build ended "NOT READY gate697 failed" (15:33Z) on a false FAIL: the content-
# aware diff vs R-ABC2 listed bulk_sharpen_pt, whose only difference was the two tercile `edges` at the last float digit
# (0.20105594997123832 vs …835; 0.6530179660732803 vs …802), with every count, carried flag, key, sha and assignment identical.
# scripts/def_metadiff.py now compares floats to 1e-12 relative (test M9). This script re-gates the FINISHED build; it never rebuilds:
#   1. the package sha must equal the one rabc2h_697.log logged at build time
#   2. novelty is re-run ("novelty rc=<n>")
#   3. gate697 --write (its "[rabc2h] ..." lines go RAW into the new log)
#   4. on PASS: the new log is copied to arm1015's run log and "READY rabc2h <sha> <pkg>" is appended to READY-science.log
# New log: work/rabc2h_697.r2.log; terminal "READY rabc2h <sha>" | "NOT READY <reason>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/697-rabc2h; W=$X/work; L1=$W/rabc2h_697.log; LOG=$W/rabc2h_697.r2.log; STEM=submissions/697rabc2h-october-recipe-abc
RUNLOG_ARM=experiments/428-zero-gene-audit/work/build_rabc2h/runrabc2h.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0; finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc)"' EXIT
[ -e "$LOG" ] && mv "$LOG" "$LOG.aside-$(date -u +%H%M%S)"
log "ready697r2 (pid $$): re-gate of the finished R-ABC2H build (def_metadiff float tolerance 1e-12; gate697 $(shasum -a 256 $X/gate697.py | cut -c1-16)); first run: $(tail -1 "$L1" | cut -c1-100)"
WANT="$(/usr/bin/grep -oE "package $STEM.vcc sha256 [0-9a-f]{64}" "$L1" | tail -1 | awk '{print $4}')"
[ -n "$WANT" ] || finish "NOT READY no package line in $L1"
GOT="$(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
[ "$GOT" = "$WANT" ] || finish "NOT READY package sha $GOT != the build's $WANT"
log "package $STEM.vcc sha256 $GOT"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate697.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate697 (re-gate) rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate697 failed again (see above)"
[ -f "$STEM.h5ad" ] && { rm -f "$STEM.h5ad"; log "removed $STEM.h5ad (package exists, sha logged, every gate passed)"; }
mkdir -p "$(dirname "$RUNLOG_ARM")"; cp "$LOG" "$RUNLOG_ARM"
printf '%s READY rabc2h %s %s\n' "$(date -u +%FT%TZ)" "$GOT" "$STEM.vcc" >> "$READYLOG"
log "evidence copied to $RUNLOG_ARM; READY line appended to $READYLOG"
finish "READY rabc2h $GOT"
