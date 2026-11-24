#!/usr/bin/env bash
# chain1007: the 10-07 rule (653's v4 of 537's draft), FIXED 2026-10-02 before any 10-06 entry exists
# (docs/NEXT-SLOT-DECISION.md, "the 10-07 / 10-08 rule, FIXED before any 10-06 entry exists").
# At 10-06 20:00Z it reads the captured entries of 446bj2zcp (10-01 s2), 495 (10-03 s1), 496a (10-05 s1), 446bj2zk16 and,
# if it was scored, 446bj2zk20 (10-06). It runs the sha-pinned v4 rule WITHOUT --arc-hook (652 F1: MZ / CPMZ never enter),
# then arms arm1007.sh with the first READY candidate of each slot's list. READY = package present; "[<tag>] PASS",
# "[<tag>] CLIFF OK" and "novelty rc=0" in its gate log; novelty re-checked now.
# R6 (653): the coordinator reads 524dr (10-06 s1) before 20:00Z. On LAG or a falsifier it writes work/HOLD-1007-s1 and
# slot 1 is not armed; slot 2's pick then goes alone. Never builds, never submits by itself.
# DRY=1: no wait, no arming; it prints the plan to the log (for testing on any day).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/660-slot-1007; W=$X/work; LOG=$W/chain1007.log; E=experiments/428-zero-gene-audit/work; S=submissions
RULE=experiments/653-rule1007-review/rule1007_v4_653.py
RULE_SHA=e00af6c7d126d4768dfc9f8c5d1541f44d80d98bba40f10fefe0383fd18128ce
DRY="${DRY:-0}"
P=""; [ "$DRY" = 1 ] && P="[DRY] "
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$P$*" >> "$LOG"; }
if [ "$DRY" != 1 ]; then
  echo $$ > $W/chain1007.pid
  log "armed (pid $$): waiting for 10-06 20:00Z"
  until [ "$(date -u +%Y%m%d%H%M)" -ge 202610062000 ]; do sleep 60; done
  # a machine that slept through the evening must not arm a day it can no longer reach in time
  if [ "$(date -u +%Y%m%d%H%M)" -ge 202610070600 ]; then log "woke after 10-07 06:00Z: too late to arm 10-07 -- nothing armed"; log "END"; exit 0; fi
fi
GOT="$(shasum -a 256 "$RULE" | cut -d' ' -f1)"
if [ "$GOT" != "$RULE_SHA" ]; then log "rule sha mismatch: $GOT != $RULE_SHA -- nothing armed"; log "END"; exit 1; fi
# the captured entry id of a tag (scripts/slot_entry_id.py: the capture after the tag's LAST ARMED line in its slot log)
eid(){ uv run python scripts/slot_entry_id.py "$1" experiments/*/work/submit100[1-6]s[12].log 2>/dev/null || echo none; }
ICP="$(eid 446bj2zcp)"; IC="$(eid 495)"; IH="$(eid 496a)"; IK16="$(eid 446bj2zk16)"; IK20="$(eid 446bj2zk20)"
log "inputs: cp(446bj2zcp) ${ICP:-none}; c(495) ${IC:-none}; h(496a) ${IH:-none}; k16 ${IK16:-none}; k20 ${IK20:-none}"
OUT="$(uv run python "$RULE" "${ICP:-none}" "${IH:-none}" "${IK16:-none}" "${IK20:-none}" "${IC:-none}" 2>>"$LOG" | tail -1)"
log "rule: $OUT"
L1="${OUT%%|*}"; REST="${OUT#*|}"; L2="${REST%%|*}"
if [ -f "$W/HOLD-1007-s1" ]; then log "HOLD-1007-s1 present (R6: 524dr read) -- slot 1 list dropped: $L1"; L1=""; fi
# rule tags -> build tags
btag(){ case "$1" in CPH) echo 446bj2zcph;; K20) echo 446bj2zk20;; 459dr) echo 459dr;; *) echo "";; esac; }
pkg(){ case "$1" in
  446bj2zcph) echo $S/446bj2zcph-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph.vcc;;
  446bj2zk20) echo $S/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc;;
  459dr) echo $S/459dr-october-path-abc-h1x3-q1.75.vcc;;
  *) echo "/nonexistent/$1";; esac; }
rl(){ case "$1" in
  459dr) echo experiments/459-h1-axis-def/work/gates459.log;;
  *) echo "$E/build_$1/run$1.log";; esac; }
# -m names from docs/MODEL-NAMES.md's queue (docs/model-names-claimed.tsv): 446bj2zcph claimed 10-02 20:01Z,
# 446bj2zk20 09-29 01:04Z, 459dr 09-28 15:18Z. -d stays empty (AGENTS §14).
name(){ case "$1" in
  446bj2zcph) echo "The reasonable man adapts himself to the world: the unreasonable one persists in trying to adapt the world to himself. Therefore all progress depends on the unreasonable man.";;
  446bj2zk20) echo "You're behind—very behind.";;
  459dr) echo "The future is here; it is just not evenly distributed.";;
  *) echo "";; esac; }
ready(){ local t="$1" r p; r="$(rl "$t")"; p="$(pkg "$t")"
  [ -s "$p" ] && /usr/bin/grep -q "^\[$t\] PASS" "$r" 2>/dev/null && /usr/bin/grep -q "^\[$t\] CLIFF OK" "$r" \
    && /usr/bin/grep -q "novelty rc=0" "$r" && [ -n "$(name "$t")" ] \
    && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1; }
A1=""; A2=""
for c in $L1; do b="$(btag "$c")"; [ -n "$b" ] || { log "s1 candidate $c has no built package -- skipped"; continue; }
  if ready "$b"; then A1=$b; break; else log "s1 candidate $c ($b) not ready -- skipped"; fi; done
for c in $L2; do b="$(btag "$c")"; [ -n "$b" ] || { log "s2 candidate $c has no built package -- skipped"; continue; }
  [ "$b" = "$A1" ] && continue
  if ready "$b"; then A2=$b; break; else log "s2 candidate $c ($b) not ready -- skipped"; fi; done
[ -z "$A1" ] && [ -n "$A2" ] && { log "slot 1 empty -- slot 2's pick $A2 goes as slot 1"; A1=$A2; A2=""; }
if [ -z "$A1" ]; then log "no eligible 10-07 arm -- nothing armed (the next session decides)"; log "END"; exit 0; fi
if [ "$DRY" = 1 ]; then log "would arm 10-07: s1 $A1 ($(pkg $A1)), s2 ${A2:-<none>}"; log "END"; exit 0; fi
nohup bash $X/arm1007.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-07" > $W/arm1007_s1.nohup 2>&1 & echo $! > $W/arm1007_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1007.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-07" > $W/arm1007_s2.nohup 2>&1 & echo $! > $W/arm1007_s2.pid; }
sleep 5
log "armed 10-07: s1 $A1 (pid $(cat $W/arm1007_s1.pid)), s2 ${A2:-<none>}"
log "END"
