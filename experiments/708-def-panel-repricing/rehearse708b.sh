#!/usr/bin/env bash
# 708 day-path test (after rehearse708's REHEARSAL PASS): run day708.sh end to end on 180's relabelled A/B/C bundle, with the
# day path's own inputs (453's def_day run: its day operator 66610be7 and cd4_compact) and FRESHLY FITTED held-out fold
# operators. This is exactly what 10-22 runs. Informational: its RECOMMENDATION should equal the rehearsal's (H1 x2.5,
# CD4 x1.5), but it is not bit-identical (another operator), and a different choice is reported, not failed.
#   nohup bash scripts/detach.sh bash experiments/708-def-panel-repricing/rehearse708b.sh > experiments/708-def-panel-repricing/work/rehearse708b.nohup 2>&1 &
# Terminal line, always: "<ts> DAYPATH <result>", then "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/708-def-panel-repricing; W=$X/work; LOG=$W/rehearse708b.log
T0="$(date -u +%FT%TZ)"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
trap 'log "END"' EXIT
log "armed (pid $$): waiting for rehearse708's END"
until newline "$T0" "$W/rehearse708.log" " END$" && ! alive $X/rehearse708.sh; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610201200 ] && { log "DAYPATH not run: rehearse708 never ended by the cap"; exit 0; }
  sleep 60
done
/usr/bin/grep -q " REHEARSAL PASS$" "$W/rehearse708.log" || { log "DAYPATH not run: the rehearsal did not PASS ($(/usr/bin/grep -E 'REHEARSAL' "$W/rehearse708.log" | tail -1 | cut -c1-160))"; exit 0; }
[ -e experiments/691-release-watcher/work/RELEASE-DETECTED ] && { log "DAYPATH skipped: release detected (day708 runs on the real bundle after R0 instead)"; exit 0; }
# 10-16 (AMENDMENT 5; coordinator): the day path fits 14 GB operators, so it starts only while no build / R-NEW chain is alive,
# and either before 10-17 01:30Z (done by ~02:30Z, before R-ABC2's ~03:00Z build) or after R-ABC2's chain (rabc2_695) exited
quiet(){ ! alive experiments/710-h1-dose-board/build710.sh && ! alive experiments/711-kup-dose/build711.sh && ! alive chain690j \
         && ! alive run690.sh && ! alive after690.sh && ! alive scripts/def_day.sh && ! alive run_arm_generic.sh; }
slot(){ [ "$(date -u +%Y%m%d%H%M)" -lt 202610170130 ] || ! alive rabc2_695.sh; }
until quiet && slot && sleep 60 && quiet && slot; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610201200 ] && { log "DAYPATH not run: no slot (queue / R-ABC2) by 10-20 12:00Z"; exit 0; }
  sleep 60
done
log "DAYPATH starting: queue clear and $( [ "$(date -u +%Y%m%d%H%M)" -lt 202610170130 ] && echo 'before 01:30Z' || echo 'R-ABC2 done')"
R705=NONE R707=NONE bash $X/day708.sh experiments/180-def-readiness/work/def-bundle experiments/453-abc-dress-rehearsal/work/run d708abcdef \
  > "$W/day708-d708abcdef.nohup" 2>&1
rc=$?
GOT="$(/usr/bin/grep -E '^RECOMMENDATION' "$W/dayof438-d708abcdef.txt" 2>/dev/null | tail -1 | sed -E 's/\([^)]*\)//')"
EXP="$(/usr/bin/grep -E '^RECOMMENDATION' "$W/dayof438-r708abcdef.txt" 2>/dev/null | tail -1 | sed -E 's/\([^)]*\)//')"
log "day708 rc $rc; its last line: $(tail -2 "$W/day708-d708abcdef.log" 2>/dev/null | head -1 | cut -c1-200)"
if [ "$GOT" = "$EXP" ] && [ -n "$GOT" ]; then log "DAYPATH same choice as the rehearsal: $GOT"
else log "DAYPATH different choice: day path '$GOT' vs rehearsal '$EXP' (report; not a failure by itself)"; fi
