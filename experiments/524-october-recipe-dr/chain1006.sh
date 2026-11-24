#!/usr/bin/env bash
# chain1006: the 10-06 rule, FIXED 2026-09-29 01:10Z before any 09-29 entry exists (docs/NEXT-SLOT-DECISION.md; bins in
# 524/PREREGISTRATION.md and 523/PREREGISTRATION.md).
#   s1: 524dr (the October recipe as resolved on 10-05, on relabelled A/B/C; read vs 453dr) if READY, else 446bj2zk16;
#   s2: 446bj2zk16 (523's K 1,600, read vs 446bj2z) when s1 is 524dr, else 446bj2zk20 (the dose partner, K 2,000), so a
#       missing 524dr puts the whole preregistered K pair on one day.
# READY = the package exists; its gate log has "[tag] PASS", "[tag] CLIFF OK", "novelty rc=0" and a "package sha256:"
# equal to the file's; the novelty check passes again now. Runs at 10-05 20:00Z, after recipe_dr524.sh's process exits
# (or at 23:30Z regardless). Never builds, never submits by itself: it arms arm1006.sh.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/524-october-recipe-dr; W=$X/work; LOG=$W/chain1006.log; E=experiments/428-zero-gene-audit/work; S=submissions
mkdir -p $W
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1006.pid
log "armed (pid $$): waiting for 10-05 20:00Z and recipe_dr524.sh's process to exit (23:30Z at the latest)"
rdr_alive(){ local p; p="$(cat $W/recipe_dr524.pid 2>/dev/null)"; [ -n "$p" ] && ps -p "$p" -o args= 2>/dev/null | /usr/bin/grep -q recipe_dr524; }
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610052000 ] && ! rdr_alive; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610052330 ] && { log "recipe_dr524 still running at 23:30Z -- deciding without it"; break; }
  sleep 60; done
log "recipe_dr524's last lines: $(tail -3 $W/recipe_dr524.log 2>/dev/null | tr '\n' ' ' | cut -c1-400)"
pkg(){ case "$1" in
  524dr) echo $S/524dr-october-recipe-abc.vcc;;
  446bj2zk16) echo $S/446bj2zk16-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1200.vcc;;
  446bj2zk20) echo $S/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc;;
  *) echo "/nonexistent/$1";; esac; }
rl(){ case "$1" in
  524dr) echo $W/gates524.log;;
  *) echo "$E/build_$1/run$1.log";; esac; }
# -m names from the MODEL-NAMES.md queue (claimed 2026-09-29 01:04Z, docs/model-names-claimed.tsv)
name(){ case "$1" in
  524dr) echo "The most dangerous worldviews are the worldviews of those who have never viewed the world.";;
  446bj2zk16) echo "I wish to see my children eating off golden plates.";;
  446bj2zk20) echo "You're behind—very behind.";;
  *) echo "";; esac; }
ready(){ local r p want got; r="$(rl "$1")"; p="$(pkg "$1")"
  [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" \
    && /usr/bin/grep -q "novelty rc=0" "$r" || return 1
  want="$(/usr/bin/grep 'package sha256:' "$r" | tail -1 | awk '{print $NF}')"; got="$(shasum -a 256 "$p" | cut -d' ' -f1)"
  [ -n "$want" ] && [ "$want" = "$got" ] || { log "$1: package sha ${got:0:16} != gate log ${want:0:16}"; return 1; }
  uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1; }
A1=""; A2=""
for c in 524dr 446bj2zk16; do if ready "$c"; then A1=$c; break; else log "s1 candidate $c not ready -- skipped"; fi; done
if [ "$A1" = 524dr ]; then L2="446bj2zk16 446bj2zk20"; else L2="446bj2zk20"; fi
for c in $L2; do [ "$c" = "$A1" ] && continue; if ready "$c"; then A2=$c; break; else log "s2 candidate $c not ready -- skipped"; fi; done
[ -z "$A1" ] && [ -n "$A2" ] && { A1=$A2; A2=""; }
[ -n "$A1" ] || { log "no eligible 10-06 arm -- nothing armed (the next session decides)"; log "END"; exit 0; }
nohup bash $X/arm1006.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-06" > $W/arm1006_s1.nohup 2>&1 & echo $! > $W/arm1006_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1006.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-06" > $W/arm1006_s2.nohup 2>&1 & echo $! > $W/arm1006_s2.pid; }
sleep 5
log "armed 10-06: s1 $A1 (pid $(cat $W/arm1006_s1.pid)), s2 ${A2:-<none>}"
log "END"
