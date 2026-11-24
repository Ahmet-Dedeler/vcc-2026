#!/usr/bin/env bash
# 698 re-gate (10-18, readiness [4592ee]). C1's build ended "NOT READY gate698 failed" (11:52Z) on a false FAIL: the content-aware
# diff vs R-ABC2H listed bulk_sharpen_pt, whose only difference was tercilec's own round-off diagnostic key_ratio_maxdev (8.66e-15 vs
# 1.02e-14) and the edges / key_ratio_c at ~1e-15 relative, with counts, carried key + sha, moved targets and assignments identical.
# scripts/def_metadiff.py now has an absolute float floor (atol 1e-12, test M10; f9b6e09d). This re-gates the FINISHED build; it never
# rebuilds:
#   0. carries the first run's KA and dump-check lines forward (KA-off, KA-on, stack dump == 697's, the dump def_day used == the lab's)
#   1. the package sha must equal the one c1811_698.log logged at build time
#   2. novelty is re-run ("novelty rc=<n>")
#   3. gate698 --write (its "[c1811] ..." lines go RAW into the new log)
#   4. on PASS: the h5ad is removed, the new log is copied to arm1015's run log and "READY c1811 <sha> <pkg>" goes to READY-science.log
# New log: work/c1811_698.r2.log; terminal "READY c1811 <sha>" | "NOT READY <reason>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/698-c1811; W=$X/work; L1=$W/c1811_698.log; LOG=$W/c1811_698.r2.log; STEM=submissions/698c1811-october-recipe-abc
RUNLOG_ARM=experiments/428-zero-gene-audit/work/build_c1811/runc1811.log; READYLOG=experiments/701-kup-callset/work/READY-science.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0; finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc)"' EXIT
[ -e "$LOG" ] && mv "$LOG" "$LOG.aside-$(date -u +%H%M%S)"
log "ready698r2 (pid $$): re-gate of the finished C1 build (def_metadiff rtol 1e-12 + atol 1e-12, $(shasum -a 256 scripts/def_metadiff.py | cut -c1-16); gate698 $(shasum -a 256 $X/gate698.py | cut -c1-16)); first run ended: $(tail -1 "$L1" | cut -c1-100)"
/usr/bin/grep -E "\[c1811\] (KA-off PASS|KA-on PASS|stack dump ==|dump check PASS)|GO: hook|hook committed|recipe check:" "$L1" >> "$LOG"
[ "$(/usr/bin/grep -cE "\[c1811\] (KA-off PASS|KA-on PASS|stack dump ==|dump check PASS)" "$L1")" = 4 ] \
  || finish "NOT READY the first run's log does not hold all four KA / dump-check PASS lines"
WANT="$(/usr/bin/grep -oE "package $STEM.vcc sha256 [0-9a-f]{64}" "$L1" | tail -1 | awk '{print $4}')"
[ -n "$WANT" ] || finish "NOT READY no package line in $L1"
GOT="$(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
[ "$GOT" = "$WANT" ] || finish "NOT READY package sha $GOT != the build's $WANT"
log "package $STEM.vcc sha256 $GOT"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate698.py --work "$W" --stem "$STEM" --chainlog "$LOG" --write >> "$LOG" 2>&1; G=$?
log "gate698 (re-gate) rc=$G"
[ "$G" = 0 ] || finish "NOT READY gate698 failed again (see above)"
[ -f "$STEM.h5ad" ] && { rm -f "$STEM.h5ad"; log "removed $STEM.h5ad (package exists, sha logged, every gate passed)"; }
mkdir -p "$(dirname "$RUNLOG_ARM")"; cp "$LOG" "$RUNLOG_ARM"
printf '%s READY c1811 %s %s\n' "$(date -u +%FT%TZ)" "$GOT" "$STEM.vcc" >> "$READYLOG"
log "evidence copied to $RUNLOG_ARM; READY line appended to $READYLOG"
finish "READY c1811 $GOT"
