#!/usr/bin/env bash
# chain1001b: the 10-01 rule as REVISED 2026-09-26 17:20Z, before any 09-27 entry exists (docs/NEXT-SLOT-DECISION.md;
# the rule itself is rule1001b.py, which adds 09-30 s2 = 458g to the inputs and 469dr, the gated October path, to the
# candidates). Waits for the 09-30 slots, applies the rule, arms arm1001.sh s1/s2. Never builds, never submits by itself.
# Replaces chain1001.sh (killed 17:2xZ 09-26).
# 2026-09-29 11:10Z: re-armed with 446bj2zmz (527) known to pkg/name (rule1001b.py's amendment, before any 09-30 entry).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/464-slot-1001; W=$X/work; LOG=$W/chain1001b.log; P27=experiments/446-q175-level/work; P30=experiments/459-h1-axis-def/work
E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1001b.pid
log "armed (pid $$): waiting for the 09-30 second slot to finish"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609300005 ] \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /arm0930\.sh|chain0930c\.sh/' | /usr/bin/grep -q . \
      && { /usr/bin/grep -qE "CAPTURED own entry|NOT SUBMITTING" $P30/submit0930s2.log 2>/dev/null \
           || { /usr/bin/grep -qE "NOT SUBMITTING" $P30/submit0930s1.log 2>/dev/null; } \
           || { /usr/bin/grep -qE "no 09-30 slot armed" $P30/chain0930c.log 2>/dev/null; } \
           || { /usr/bin/grep -qE "armed 09-30: s1 [A-Za-z0-9]+ \(pid [0-9]+\), s2 <none>" $P30/chain0930c.log 2>/dev/null \
                && /usr/bin/grep -qE "CAPTURED own entry|NOT SUBMITTING" $P30/submit0930s1.log 2>/dev/null; }; }; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609301800 ] && { log "09-30 slot 2 never finished by 18:00Z -- proceeding with what exists"; break; }
  sleep 60; done
sleep 30
ID27="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P27/submit0927s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
ID1="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P30/submit0930s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
TAG1="$(/usr/bin/grep -o 'ARMED slot 1: [A-Za-z0-9]*' $P30/submit0930s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
ID2="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P30/submit0930s2.log 2>/dev/null | tail -1 | awk '{print $4}')"
TAG2="$(/usr/bin/grep -o 'ARMED slot 2: [A-Za-z0-9]*' $P30/submit0930s2.log 2>/dev/null | tail -1 | awk '{print $4}')"
log "inputs: 09-27 s1 (446bj2) ${ID27:-<none>}; 09-30 s1 (${TAG1:-?}) ${ID1:-<none>}; 09-30 s2 (${TAG2:-?}) ${ID2:-<none>}"
RULE="$(uv run python $X/rule1001b.py "${ID27:-none}" "${ID1:-none}" "${TAG1:-none}" "${ID2:-none}" "${TAG2:-none}" 2>&1 | tail -1)"
log "rule: $RULE"
pkg(){ case "$1" in
  461dr) echo $S/461dr-october-path-arcprog-abc-q1.75.vcc;;
  460dr) echo $S/460dr-october-path-lodo-abc-q1.75.vcc;;
  469dr) echo $S/469dr-october-path-arcprog-gate33-abc-q1.75.vcc;;
  458e) echo $S/458e-cap15-partial50-arch1eblend-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc;;
  446bj2zmz) echo $S/446bj2zmz-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-hkmz05le1.vcc;;
  *) echo "/nonexistent/$1";; esac; }
name(){ case "$1" in
  461dr) echo "Microsoft has had clear competitors in the past. It's a good thing we have museums to document that. (Bill Gates)";; 460dr) echo "Hate the chains that hold you back. Hate the doubts that whisper lies. Hate the part of you that still seeks comfort in mediocrity. Use it. Let it drive you. Become the force that reshapes the world, or be forgotten by it.";;
  469dr) echo "When something truly matters and you possess the ability to create change, it becomes your moral duty to do so.";; 458e) echo "I would like to own all the land in this valley, as far as I can see.";; 446bj2zmz) echo "Sun Tzu, The Art of Being A Decent Human Being";; *) echo "";; esac; }
rl(){ case "$1" in
  461dr) echo experiments/461-arcprog-october-path/work/gates461.log;; 460dr) echo experiments/460-lodo-dayof/work/gates460.log;;
  469dr) echo experiments/469-gated-october-path/work/gates469.log;;
  *) echo "$E/build_$1/run$1.log";; esac; }
ready(){ local r p; r="$(rl "$1")"; p="$(pkg "$1")"; [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null \
         && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" && /usr/bin/grep -q "novelty rc=0" "$r" \
         && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1; }
A1=""; A2=""
for c in ${RULE%%|*}; do
  [ "$c" = "-" ] && continue
  if ready "$c"; then [ -z "$A1" ] && A1=$c || { [ -z "$A2" ] && A2=$c; }; else log "$c not ready (gate/package/novelty) -- skipped"; fi
done
[ -n "$A1" ] || { log "no eligible 10-01 arm -- nothing armed (the next session decides)"; log "END"; exit 0; }
nohup bash $X/arm1001.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-01" > $W/arm1001_s1.nohup 2>&1 & echo $! > $W/arm1001_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1001.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-01" > $W/arm1001_s2.nohup 2>&1 & echo $! > $W/arm1001_s2.pid; }
sleep 5
log "armed 10-01: s1 $A1 (pid $(cat $W/arm1001_s1.pid)), s2 ${A2:-<none>}"
log "END"
