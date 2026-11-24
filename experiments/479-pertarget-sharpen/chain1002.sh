#!/usr/bin/env bash
# chain1002: the 10-02 rule, FIXED 2026-09-27 before any 09-28 entry exists (docs/NEXT-SLOT-DECISION.md; 479's
# PREREGISTRATION: s1 P2g, s2 S_P2, a same-day pair vs 446bj2). The rule reads no 10-01 entry, so this arms at
# 10-01 20:00Z (after the 10-01 slots). Candidates in order: 479p2g, 479sp2, then the fallbacks 495 (cap-6 key; its own
# PREREGISTRATION) and 471c (the re-levelled union swap; 471's addendum). The first two READY arms are armed.
# READY = package present, config diff PASS, CLIFF OK, novelty rc 0 in the build log, novelty re-checked now, and for the
# re-levelled arms every context's ntclevel rho within 0.5 % of 0.9793. Never builds, never submits by itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/479-pertarget-sharpen; W=$X/work; LOG=$W/chain1002.log; E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1002.pid
log "armed (pid $$): waiting for 10-01 20:00Z"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610012000 ]; do sleep 60; done
pkg(){ case "$1" in
  479p2g) echo $S/479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5.vcc;;
  479sp2) echo $S/479sp2-cap15-partial50-clean4q80-pertscale3-ballast-sharpensp2-stacklodo05load-ctxlevelr-mde010coarse-zr5.vcc;;
  495) echo $S/495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc;;
  471c) echo $S/471c-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevel-mde010coarse-zr5.vcc;;
  *) echo "/nonexistent/$1";; esac; }
# -m names: the MODEL-NAMES.md queue, claimed 2026-09-27 06:22Z (docs/model-names-claimed.tsv)
name(){ case "$1" in
  479p2g) echo "Be the change you wanna see in the world.";;
  479sp2) echo "You can also commit injustice by doing nothing. (Marcus Aurelius)";;
  495) echo "taking action is the best action—anytime. dont let it stay in your mind, go grind!!";;
  471c) echo "If you don't believe it or don't get it, I don't have the time to try to convince you, sorry. (Satoshi)";;
  *) echo "";; esac; }
nlj(){ case "$1" in
  479p2g) echo $E/build_479p2g/the-lyophilizer-wrote-a-cookbook.ntclevel.json;;
  479sp2) echo $E/build_479sp2/the-magnetic-stirrer-ran-for-office.ntclevel.json;;
  495) echo $E/build_495/the-water-bath-filed-for-divorce.ntclevel.json;;
  471c) echo $E/build_471c/the-ph-meter-joined-a-union.ntclevel.json;; esac; }
rho_ok(){ uv run python - "$(nlj "$1")" <<'PY'
import json, sys
per = json.load(open(sys.argv[1]))
per = per.get("per_context", per)
r = {c: float(v["rho"]) for c, v in per.items() if isinstance(v, dict) and "rho" in v}
ok = len(r) == 3 and all(abs(x / 0.9793 - 1) <= 0.005 for x in r.values())
print(r, "OK" if ok else "OUT OF BAND"); sys.exit(0 if ok else 1)
PY
}
ready(){ local r p; r="$E/build_$1/run$1.log"; p="$(pkg "$1")"; [ -s "$p" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null \
         && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" && /usr/bin/grep -q "novelty rc=0" "$r" \
         && uv run python scripts/check_submission_novelty.py check "$p" > /dev/null 2>&1 && rho_ok "$1" >> "$LOG" 2>&1; }
A1=""; A2=""
for c in 479p2g 479sp2 495 471c; do
  if ready "$c"; then [ -z "$A1" ] && A1=$c || { [ -z "$A2" ] && A2=$c; }; else log "$c not ready (gate/package/novelty/rho) -- skipped"; fi
done
[ -n "$A1" ] || { log "no eligible 10-02 arm -- nothing armed (the next session decides)"; log "END"; exit 0; }
nohup bash $X/arm1002.sh 1 "$A1" "$E/build_$A1/run$A1.log" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-02" > $W/arm1002_s1.nohup 2>&1 & echo $! > $W/arm1002_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1002.sh 2 "$A2" "$E/build_$A2/run$A2.log" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-02" > $W/arm1002_s2.nohup 2>&1 & echo $! > $W/arm1002_s2.pid; }
sleep 5
log "armed 10-02: s1 $A1 (pid $(cat $W/arm1002_s1.pid)), s2 ${A2:-<none>}"
log "END"
