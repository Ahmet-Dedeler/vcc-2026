#!/usr/bin/env bash
# chain0930b: the 09-30 rule as REVISED 2026-09-26 10:35Z, before any 09-27 entry exists (docs/NEXT-SLOT-DECISION.md):
# s1 = 458d (Arc's H1 programme from ALL public cells as every context's block prior, vs 446b) if built and gated, else 458
# (the harness-bank programme); s2 = 459dr (the H1 x3 twin of 453, vs 453) if its gate log passes. Independent of every
# earlier read; arms once the 09-29 waiters are gone. Never builds, never submits itself. Replaces chain0930.sh (killed).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/459-h1-axis-def; W=$X/work; LOG=$W/chain0930b.log; C29=experiments/457-cd4-half-def/work/chain0929.log
E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain0930b.pid
log "armed (pid $$): waiting for the 09-29 slots to finish"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609290005 ] \
      && /usr/bin/grep -qE "armed 09-29|no 09-29 slot armed|stopping" "$C29" 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /arm0929\.sh|chain0929\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609291800 ] && { log "09-29 slots not done by 18:00Z 09-29 -- arming 09-30 anyway"; break; }
  sleep 60; done
pkg(){ case "$1" in
  458d) echo $S/458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc;;
  458) echo $S/458-cap15-partial50-arch1prog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc;;
  459dr) echo $S/459dr-october-path-abc-h1x3-q1.75.vcc;;
  *) echo "/nonexistent/$1";; esac; }
name(){ case "$1" in 458d) echo "the pipette learned to meditate";; 458) echo "the autoclave read the manual";; 459dr) echo "the spectrometer took a nap";; *) echo "";; esac; }
rl(){ [ "$1" = "459dr" ] && echo $W/gates459.log || echo "$E/build_$1/run$1.log"; }
ready(){ local r p; r="$(rl "$1")"; p="$(pkg "$1")"; [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null \
         && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" && /usr/bin/grep -q "novelty rc=0" "$r" \
         && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1; }
A1=""; A2=""
if ready 458d; then A1=458d; ready 459dr && A2=459dr
elif ready 458; then A1=458; log "458d not ready -> slot 1 is 458 (the harness-bank programme)"; ready 459dr && A2=459dr
elif ready 459dr; then A1=459dr; log "458 not ready -> slot 1 is 459, slot 2 empty"
else log "none of 458d, 458, 459 is ready -- no 09-30 slot armed"; log "END"; exit 0; fi
nohup bash $X/arm0930.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 09-30" > $W/arm0930_s1.nohup 2>&1 & echo $! > $W/arm0930_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm0930.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 09-30" > $W/arm0930_s2.nohup 2>&1 & echo $! > $W/arm0930_s2.pid; }
sleep 5
log "armed 09-30: s1 $A1 (pid $(cat $W/arm0930_s1.pid)), s2 ${A2:-<none>}"
log "END"
