#!/usr/bin/env bash
# chain1008 (659): the 10-08 slot chain. A DRAFT for the coordinator, UNARMED until the coordinator arms it.
# Modelled on 660/chain1007.sh. At 10-07 20:00Z it:
#   1. checks the frozen rule's sha (653's v4; the same pin as chain1007);
#   2. reads the captured entries (scripts/slot_entry_id.py over experiments/*/work/submit100[1-7]s[12].log):
#        cp 446bj2zcp (10-01), h 496a (10-05), k16 446bj2zk16, k20 446bj2zk20 (scored 10-04 or later), c 495 (10-03),
#        cph 446bj2zcph (10-07);
#   3. derives --scored=<rule tags scored on 10-07> from 660's submit1007s1/s2.log. A slot's LAST "ARMED slot <n>: <tag>"
#      run that reached "SUBMITTED ok" or "CAPTURED own entry" scored <tag> (a submitted but uncaptured entry counts);
#   4. runs the rule WITHOUT --arc-hook and WITH --scored=...: its "10-08:" list mapped to build tags, then the queue tail
#      from tail1008.tsv (630arm -> 524dr2 -> allsw658, read now, at decision time), is each slot's candidate list;
#   5. arms arm1008.sh with the first READY candidate of each slot (slot 2 never takes slot 1's pick).
# READY:
#   - a built rule-list arm (446bj2zcph, 446bj2zk20, 459dr): package present; its gate log has "[tag] PASS",
#     "[tag] CLIFF OK" and "novelty rc=0"; a claimed name (docs/model-names-claimed.tsv); the novelty check passes now;
#   - a tail arm (a tail1008.tsv row): its DRIVER log's last line is "<ts> READY <tag> <sha>"; the driver's
#     "package <path> sha256 <sha>" line names the package and shasum(package) == <sha>; its gate log has
#     "[gate_tag] PASS", "[gate_tag] CLIFF OK", "novelty rc=0" and "package sha256: <sha>"; a claimed name; novelty passes.
# HOLD-1008-s1 / HOLD-1008-s2 in work/ (written by the coordinator before 20:00Z 10-07) drop that slot's list (R6 style).
# QUEUEFIRST-1008 in work/ (opt-in): the queue tail goes before the rule's own fallbacks (459dr, ALLSW, S1); see step 4.
# NOGO-1008-s<n> stops a slot inside arm1008.sh. Every path ends with "<ts> END". bash 3.2 / BSD grep safe.
# DRY=1: no wait and no arming; the plan goes to the log. Test-only overrides, honoured ONLY with DRY=1:
#   CH1008_W (work dir), CH1008_S1007 (dir of the submit1007 logs), CH1008_TAIL (spec), CH1008_CLAIMS (names),
#   CH1008_NOVELTY (the novelty command, e.g. "true").
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/659-chain1008; W=$X/work; E=experiments/428-zero-gene-audit/work; S=submissions
S1007=experiments/660-slot-1007/work; TAIL=$X/tail1008.tsv; CLAIMS=docs/model-names-claimed.tsv
NOVELTY="uv run python scripts/check_submission_novelty.py check"
RULE=experiments/653-rule1007-review/rule1007_v4_653.py
RULE_SHA=e00af6c7d126d4768dfc9f8c5d1541f44d80d98bba40f10fefe0383fd18128ce
DRY="${DRY:-0}"
if [ "$DRY" = 1 ]; then
  W="${CH1008_W:-$W}"; S1007="${CH1008_S1007:-$S1007}"; TAIL="${CH1008_TAIL:-$TAIL}"; CLAIMS="${CH1008_CLAIMS:-$CLAIMS}"
  NOVELTY="${CH1008_NOVELTY:-$NOVELTY}"
fi
mkdir -p "$W"; LOG=$W/chain1008.log
P=""; [ "$DRY" = 1 ] && { P="[DRY] "; LOG=$W/dry-chain1008.log; }   # a DRY run never writes the real chain's log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$P$*" >> "$LOG"; }
if [ "$DRY" != 1 ]; then
  echo $$ > $W/chain1008.pid
  log "armed (pid $$): waiting for 10-07 20:00Z"
  until [ "$(date -u +%Y%m%d%H%M)" -ge 202610072000 ]; do sleep 60; done
  # a machine that slept through the evening must not arm a day it can no longer reach in time
  if [ "$(date -u +%Y%m%d%H%M)" -ge 202610080600 ]; then log "woke after 10-08 06:00Z: too late to arm 10-08 -- nothing armed"; log "END"; exit 0; fi
fi
GOT="$(shasum -a 256 "$RULE" | cut -d' ' -f1)"
if [ "$GOT" != "$RULE_SHA" ]; then log "rule sha mismatch: $GOT != $RULE_SHA -- nothing armed"; log "END"; exit 1; fi

# ---- inputs: the captured entry of each tag (the capture after the tag's LAST ARMED line in its slot log)
eid(){ uv run python scripts/slot_entry_id.py "$1" experiments/*/work/submit100[1-7]s[12].log 2>/dev/null || echo none; }
ICP="$(eid 446bj2zcp)"; IC="$(eid 495)"; IH="$(eid 496a)"; IK16="$(eid 446bj2zk16)"; IK20="$(eid 446bj2zk20)"; ICPH="$(eid 446bj2zcph)"
log "inputs: cp(446bj2zcp) ${ICP:-none}; h(496a) ${IH:-none}; k16 ${IK16:-none}; k20 ${IK20:-none}; c(495) ${IC:-none}; cph(446bj2zcph) ${ICPH:-none}"

# ---- --scored: build tags whose LAST armed run on 10-07 submitted, mapped to rule tags
rtag(){ case "$1" in 446bj2zcph) echo CPH;; 446bj2zk20) echo K20;; 459dr) echo 459dr;; allsw658) echo ALLSW;;
                     446bj2zcphk20) echo CPHK20;; *) echo "";; esac; }
scored_tag(){ # $1 = a submit1007 log: print the build tag of its last ARMED run if that run submitted
  [ -s "$1" ] || return 0
  awk '/ARMED slot [12]: / { t = $0; sub(/.*ARMED slot [12]: /, "", t); sub(/ .*/, "", t); ok = 0 }
       /SUBMITTED ok|CAPTURED own entry/ { if (t != "") ok = 1 }
       END { if (t != "" && ok) print t }' "$1"; }
SCORED=""; SCORED_B=""
for f in "$S1007/submit1007s1.log" "$S1007/submit1007s2.log"; do
  b="$(scored_tag "$f")"; [ -n "$b" ] || continue
  SCORED_B="$SCORED_B $b"; r="$(rtag "$b")"
  [ -n "$r" ] && SCORED="${SCORED:+$SCORED,}$r"
done
log "scored on 10-07 (build tags):${SCORED_B:- none}; --scored=$SCORED"

OUT="$(uv run python "$RULE" "${ICP:-none}" "${IH:-none}" "${IK16:-none}" "${IK20:-none}" "${IC:-none}" "${ICPH:-none}" "--scored=$SCORED" 2>>"$LOG" | tail -1)"
log "rule: $OUT"
D8="$(printf '%s\n' "$OUT" | awk -F'|' '{print $3}' | sed -E 's/^10-08: ?//')"
case "$OUT" in *"|10-08:"*) ;; *) log "the rule printed no 10-08 list -- nothing armed"; log "END"; exit 1;; esac

# ---- candidates: the rule's 10-08 list (rule tags -> build tags), then the tail in tail1008.tsv's row order
btag(){ case "$1" in CPH) echo 446bj2zcph;; K20) echo 446bj2zk20;; 459dr) echo 459dr;; ALLSW) echo allsw658;;
                     CPHK20) echo 446bj2zcphk20;; *) echo "";; esac; }
tail_tags(){ [ -s "$TAIL" ] && awk -F'\t' '$0 !~ /^#/ && NF >= 5 { print $1 }' "$TAIL"; }
tail_get(){ [ -s "$TAIL" ] && awk -F'\t' -v t="$1" -v k="$2" '$0 !~ /^#/ && $1 == t { v = $k } END { print v }' "$TAIL"; }
is_tail(){ tail_tags | /usr/bin/grep -qx "$1"; }
CAND=""
add(){ case " $CAND " in *" $1 "*) ;; *) CAND="$CAND $1";; esac; }
if [ -f "$W/QUEUEFIRST-1008" ]; then
  # opt-in (the coordinator creates work/QUEUEFIRST-1008 before 20:00Z 10-07): the rule's DE follow-ups (CPHK20, CPMZ, K20,
  # MZ) first, then the queue tail in tail1008.tsv's order, then the rule's own fallbacks (459dr, ALLSW, S1). Without the
  # file the literal order holds: the rule's whole 10-08 list (which ends with its fallbacks 459dr, ALLSW), then the tail.
  log "QUEUEFIRST-1008 present: rule follow-ups, then the queue tail, then the rule's fallbacks"
  for c in $D8; do case "$c" in 459dr|ALLSW|S1) ;; *) b="$(btag "$c")"; if [ -n "$b" ]; then add "$b"; else log "rule candidate $c has no build tag -- skipped"; fi;; esac; done
  for t in $(tail_tags); do add "$t"; done
  for c in $D8; do case "$c" in 459dr|ALLSW|S1) b="$(btag "$c")"; if [ -n "$b" ]; then add "$b"; else log "rule candidate $c has no build tag -- skipped"; fi;; esac; done
else
  for c in $D8; do b="$(btag "$c")"; if [ -n "$b" ]; then add "$b"; else log "rule candidate $c has no build tag -- skipped"; fi; done
  for t in $(tail_tags); do add "$t"; done
fi
log "candidates (in order):$CAND"

# ---- readiness
name(){ [ -s "$CLAIMS" ] && awk -F'\t' -v t="$1" '$2 == t { n = $3 } END { print n }' "$CLAIMS"; }
pkg(){ case "$1" in
  446bj2zcph) echo $S/446bj2zcph-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph.vcc;;
  446bj2zk20) echo $S/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc;;
  459dr) echo $S/459dr-october-path-abc-h1x3-q1.75.vcc;;
  *) if is_tail "$1"; then tail_pkg "$1"; else echo "/nonexistent/$1"; fi;; esac; }
rl(){ case "$1" in
  459dr) echo experiments/459-h1-axis-def/work/gates459.log;;
  446bj2zcph|446bj2zk20) echo "$E/build_$1/run$1.log";;
  *) if is_tail "$1"; then tail_get "$1" 4; else echo "/nonexistent/$1.log"; fi;; esac; }
gtag(){ if is_tail "$1"; then tail_get "$1" 5; else echo "$1"; fi; }
tail_sha(){ # the sha of the driver's terminal READY line, or nothing
  local d l; d="$(tail_get "$1" 3)"; [ -n "$d" ] && [ "$d" != TBD ] && [ -s "$d" ] || return 0
  l="$(awk 'NF { x = $0 } END { print x }' "$d")"
  printf '%s\n' "$l" | sed -nE "s/^[0-9T:Z.-]+ READY $1 ([0-9a-f]{64})$/\1/p"; }
tail_pkg(){ # the package the driver logged with the READY sha
  local d s; d="$(tail_get "$1" 3)"; s="$(tail_sha "$1")"; [ -n "$s" ] || { echo "/nonexistent/$1"; return 0; }
  /usr/bin/grep -E "package [^ ]+ sha256 $s" "$d" | tail -1 | sed -E "s/.*package ([^ ]+) sha256 $s.*/\1/"; }
ready(){ local t="$1" r p g n s
  p="$(pkg "$t")"; r="$(rl "$t")"; g="$(gtag "$t")"; n="$(name "$t")"
  if is_tail "$t"; then
    s="$(tail_sha "$t")"
    [ -n "$s" ] || { log "  $t: driver log ($(tail_get "$t" 3)) does not end with 'READY $t <sha>'"; return 1; }
    [ -s "$p" ] || { log "  $t: no package named by the driver for sha ${s:0:16}"; return 1; }
    [ "$(shasum -a 256 "$p" | cut -d' ' -f1)" = "$s" ] || { log "  $t: package sha != the driver's READY sha"; return 1; }
    /usr/bin/grep -q "package sha256: $s" "$r" 2>/dev/null || { log "  $t: gate log $r has no 'package sha256: ${s:0:16}...'"; return 1; }
  else
    [ -s "$p" ] || { log "  $t: no package at $p"; return 1; }
  fi
  /usr/bin/grep -q "^\[$g\] PASS" "$r" 2>/dev/null || { log "  $t: no '[$g] PASS' in $r"; return 1; }
  /usr/bin/grep -q "^\[$g\] CLIFF OK" "$r" || { log "  $t: no '[$g] CLIFF OK' in $r"; return 1; }
  /usr/bin/grep -q "novelty rc=0" "$r" || { log "  $t: no 'novelty rc=0' in $r"; return 1; }
  [ -n "$n" ] || { log "  $t: no claimed name in $CLAIMS (scripts/next_model_name.py --claim $t)"; return 1; }
  $NOVELTY "$p" > /dev/null 2>&1 || { log "  $t: the novelty check refuses $p (already scored?)"; return 1; }
  return 0; }

A1=""; A2=""
if [ -f "$W/HOLD-1008-s1" ]; then log "HOLD-1008-s1 present -- slot 1 list dropped"
else for c in $CAND; do if ready "$c"; then A1=$c; break; else log "s1 candidate $c not ready -- skipped"; fi; done; fi
if [ -f "$W/HOLD-1008-s2" ]; then log "HOLD-1008-s2 present -- slot 2 list dropped"
else for c in $CAND; do [ "$c" = "$A1" ] && continue
  if ready "$c"; then A2=$c; break; else log "s2 candidate $c not ready -- skipped"; fi; done; fi
[ -z "$A1" ] && [ -n "$A2" ] && { log "slot 1 empty -- slot 2's pick $A2 goes as slot 1"; A1=$A2; A2=""; }
if [ -z "$A1" ]; then log "no eligible 10-08 arm -- nothing armed (the next session decides)"; log "END"; exit 0; fi
if [ "$DRY" = 1 ]; then
  log "would arm 10-08: s1 $A1 ($(pkg $A1); gate log $(rl $A1) [$(gtag $A1)]; name '$(name $A1)'), s2 ${A2:-<none>}${A2:+ ($(pkg $A2); gate log $(rl $A2) [$(gtag $A2)]; name '$(name $A2)')}"
  log "END"; exit 0; fi
nohup bash $X/arm1008.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-08" "$(gtag $A1)" > $W/arm1008_s1.nohup 2>&1 & echo $! > $W/arm1008_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1008.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-08" "$(gtag $A2)" > $W/arm1008_s2.nohup 2>&1 & echo $! > $W/arm1008_s2.pid; }
sleep 5
log "armed 10-08: s1 $A1 (pid $(cat $W/arm1008_s1.pid)), s2 ${A2:-<none>}"
log "END"
