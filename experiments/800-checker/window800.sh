#!/bin/bash
# window800.sh — is a lab reserved window open now? (charter Amendment 1 + the coordinator's 10-16 00:45Z amendment)
# Prints "WINDOW <reason>" and exits 1 when open, "CLEAR" and exits 0 otherwise. In a window: no local lab job, no Mac
# transfer > 1 GB, at most 2 lab agents (the lead counts as one).
#   1. daily 23:30-01:45Z;
#   2. the release window: from the moment D/E/F's release is detected (data/challenge-def/manifest.json exists, or the
#      coordinator's guard_release.log says "RELEASE detected") until R0 is submitted (submit-R0.log says "SUBMITTED ok");
#   3. the fixed window 10-21 20:00Z -> 10-22 14:00Z, extended until R0 is submitted if it slips.
#   4. a coordinator HOLD on local jobs only: while experiments/800-checker/HOLD exists and the log it names lacks its
#      pattern, print "WINDOW hold: ..." (so no worker starts a local job) but exit 0, because a hold does not limit the
#      number of lab agents (windowmgr, lab_subagent) or EC2/transfers (labcloud). HOLD lines: log=<path> until=<grep -F
#      pattern> reason=<text>, optional release_file=<path> and release_utc=<ISO Z>; whichever releases first lifts it. The
#      lead removes the file once it is released.
R=~/Code/virtual-cell-challenge-2026
now=$(date -u +%s); hm=$((10#$(date -u +%H%M)))
if [ "$hm" -ge 2330 ] || [ "$hm" -lt 145 ]; then echo "WINDOW daily 23:30-01:45Z"; exit 1; fi
r0_done=0
grep -q "SUBMITTED ok" "$R/experiments/670-coord-1015/work/submit-R0.log" 2>/dev/null && r0_done=1
released=0
[ -f "$R/data/challenge-def/manifest.json" ] && released=1
grep -q "RELEASE detected" "$R/experiments/670-coord-1015/work/guard_release.log" 2>/dev/null && released=1
if [ "$released" -eq 1 ] && [ "$r0_done" -eq 0 ]; then echo "WINDOW release detected, R0 not yet submitted"; exit 1; fi
t0=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-10-21T20:00:00Z" +%s); t1=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-10-22T14:00:00Z" +%s)
if [ "$now" -ge "$t0" ] && { [ "$now" -lt "$t1" ] || [ "$r0_done" -eq 0 ]; }; then
  [ "$r0_done" -eq 1 ] && [ "$now" -ge "$t1" ] || { echo "WINDOW 10-21 20:00Z -> 10-22 14:00Z (or until R0 is submitted)"; exit 1; }
fi
H="$R/experiments/800-checker/HOLD"
if [ -f "$H" ]; then
  hl=$(sed -n 's/^log=//p' "$H"); hu=$(sed -n 's/^until=//p' "$H"); hr=$(sed -n 's/^reason=//p' "$H")
  hf=$(sed -n 's/^release_file=//p' "$H"); ht=$(sed -n 's/^release_utc=//p' "$H")
  hrel=0
  [ -n "$hl" ] && [ -n "$hu" ] && /usr/bin/grep -qF -- "$hu" "$R/$hl" 2>/dev/null && hrel=1
  [ -n "$hf" ] && [ -e "$R/$hf" ] && hrel=1
  [ -n "$ht" ] && [ "$now" -ge "$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "$ht" +%s)" ] && hrel=1
  if [ "$hrel" -eq 0 ]; then
    echo "WINDOW hold (local jobs only): $hr"; exit 0
  fi
fi
echo "CLEAR"; exit 0
