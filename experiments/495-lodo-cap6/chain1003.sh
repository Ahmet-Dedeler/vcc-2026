#!/usr/bin/env bash
# chain1003: the 10-03 rule, FIXED 2026-09-27 before any 09-28 entry exists (docs/NEXT-SLOT-DECISION.md). Two bulk
# questions vs 446bj2 that do not depend on 10-02's reads: s1 495 (the LODO key at cap 6; PREREGISTRATION 706ad5f0),
# s2 471c (the h1_2025_full union swap, re-levelled; 471's PREREGISTRATION addenda). Fallback if either is not READY:
# 497dr (the October path + percentile cap + S_P2, vs 460dr; gate497.py). Applies at 10-02 20:00Z and arms arm1003.sh.
# READY = package present, config diff PASS, CLIFF OK, novelty rc 0 in its log, novelty re-checked now, and for the
# production arms every context's ntclevel rho within 0.5 % of 0.9793. Never builds, never submits by itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/495-lodo-cap6; W=$X/work; LOG=$W/chain1003.log; E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1003.pid
log "armed (pid $$): waiting for 10-02 20:00Z"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610022000 ]; do sleep 60; done
pkg(){ case "$1" in
  495) echo $S/495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc;;
  471c) echo $S/471c-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevel-mde010coarse-zr5.vcc;;
  497dr) echo $S/497dr-october-path-lodopct-spt-abc.vcc;;
  *) echo "/nonexistent/$1";; esac; }
# -m names: the MODEL-NAMES.md queue, claimed 2026-09-27 06:22Z (docs/model-names-claimed.tsv)
name(){ case "$1" in
  495) echo "taking action is the best action—anytime. dont let it stay in your mind, go grind!!";;
  471c) echo "If you don't believe it or don't get it, I don't have the time to try to convince you, sorry. (Satoshi)";;
  497dr) echo "You're not respawning. What are you waiting for?";;
  *) echo "";; esac; }
rl(){ case "$1" in 497dr) echo experiments/497-october-switches/work/gates497.log;; *) echo "$E/build_$1/run$1.log";; esac; }
nlj(){ case "$1" in
  495) echo $E/build_495/the-water-bath-filed-for-divorce.ntclevel.json;;
  471c) echo $E/build_471c/the-ph-meter-joined-a-union.ntclevel.json;; esac; }
rho_ok(){ [ "$1" = 497dr ] && return 0; uv run python - "$(nlj "$1")" <<'PY'
import json, sys
per = json.load(open(sys.argv[1]))
per = per.get("per_context", per)
r = {c: float(v["rho"]) for c, v in per.items() if isinstance(v, dict) and "rho" in v}
ok = len(r) == 3 and all(abs(x / 0.9793 - 1) <= 0.005 for x in r.values())
print(r, "OK" if ok else "OUT OF BAND"); sys.exit(0 if ok else 1)
PY
}
ready(){ local r p; r="$(rl "$1")"; p="$(pkg "$1")"
  if [ "$1" = 497dr ]; then uv run python experiments/497-october-switches/gate497.py --write >> "$LOG" 2>&1 || return 1; fi
  [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null \
         && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" && /usr/bin/grep -q "novelty rc=0" "$r" \
         && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1 && rho_ok "$1" >> "$LOG" 2>&1; }
A1=""; A2=""
for c in 495 471c 497dr; do
  if ready "$c"; then [ -z "$A1" ] && A1=$c || { [ -z "$A2" ] && A2=$c; }; else log "$c not ready (gate/package/novelty/rho) -- skipped"; fi
done
[ -n "$A1" ] || { log "no eligible 10-03 arm -- nothing armed (the next session decides)"; log "END"; exit 0; }
nohup bash $X/arm1003.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-03" > $W/arm1003_s1.nohup 2>&1 & echo $! > $W/arm1003_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1003.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-03" > $W/arm1003_s2.nohup 2>&1 & echo $! > $W/arm1003_s2.pid; }
sleep 5
log "armed 10-03: s1 $A1 (pid $(cat $W/arm1003_s1.pid)), s2 ${A2:-<none>}"
log "END"
