#!/usr/bin/env bash
# Applies "the 09-29 rule, FIXED 2026-09-26 before any 09-27 entry exists" (docs/NEXT-SLOT-DECISION.md): s1 = 457 (the CD4
# x1.5 twin of 453, vs 453) if its gate log passes, s2 = 446bd2 (DLD-1 union 2, vs 446b) if built and never submitted. The
# pair does not depend on any 09-27/28 read, so it arms once the 09-28 waiters are gone. Never builds, never submits itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/457-cd4-half-def; W=$X/work; LOG=$W/chain0929.log; C28=experiments/456-zs-riders/work/chain0928.log
E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain0929.pid
log "armed (pid $$): waiting for the 09-28 slots to finish"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609280005 ] \
      && /usr/bin/grep -qE "armed 09-28|stopping" "$C28" 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /arm0928\.sh|chain0928\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609281800 ] && { log "09-28 slots not done by 18:00Z 09-28 -- arming 09-29 anyway"; break; }
  sleep 60; done
pkg(){ case "$1" in
  457dr) echo $S/457dr-october-path-abc-cd4x15-q1.75.vcc;;
  446bd2) echo $S/446bd2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackdld1u2load-ctxlevelm-mde010coarse-zr5.vcc;;
  *) echo "/nonexistent/$1";; esac; }
name(){ case "$1" in 457dr) echo "Don't compare the outputs—compare the inputs.";; 446bd2) echo "The task is not so much to see what no one has yet seen, but to think what nobody has yet thought about that which everybody sees.";; *) echo "";; esac; }
rl(){ [ "$1" = "457dr" ] && echo $W/gates457.log || echo "$E/build_$1/run$1.log"; }
ready(){ local r p; r="$(rl "$1")"; p="$(pkg "$1")"; [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null \
         && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" && /usr/bin/grep -q "novelty rc=0" "$r" \
         && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1; }
A1=""; A2=""
if ready 457dr; then A1=457dr; ready 446bd2 && A2=446bd2
elif ready 446bd2; then A1=446bd2; log "457 not ready -> slot 1 is 446bd2, slot 2 empty"
else log "neither 457 nor 446bd2 is ready -- no 09-29 slot armed"; log "END"; exit 0; fi
nohup bash $X/arm0929.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 09-29" > $W/arm0929_s1.nohup 2>&1 & echo $! > $W/arm0929_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm0929.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 09-29" > $W/arm0929_s2.nohup 2>&1 & echo $! > $W/arm0929_s2.pid; }
sleep 5
log "armed 09-29: s1 $A1 (pid $(cat $W/arm0929_s1.pid)), s2 ${A2:-<none>}"
log "END"
