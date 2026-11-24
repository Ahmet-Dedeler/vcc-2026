#!/usr/bin/env bash
# chain1005: the 10-05 rule, FIXED 2026-09-27 before any 09-28 entry exists (rule1005.py; docs/NEXT-SLOT-DECISION.md).
# At 10-04 20:00Z it reads the captured entries of 453dr (09-28), the 479 pair (10-02) and 495 (10-03), applies the rule,
# then arms arm1005.sh with the first READY candidate of each slot's list. READY = package present, gates PASS (config
# diff / CLIFF OK or the October-path gates), novelty rc 0 in its log, novelty re-checked now, and for the production
# re-levelled arms every context's ntclevel rho within 0.5 % of 0.9793. Never builds, never submits by itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/496-uncalled-sign-hook; W=$X/work; LOG=$W/chain1005.log; E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1005.pid
log "armed (pid $$): waiting for 10-04 20:00Z"
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610042000 ]; do sleep 60; done
# the captured entry id of a tag (scripts/slot_entry_id.py: the capture that follows the tag's LAST ARMED line in its
# slot log, so a re-armed slot cannot lend its entry to another tag; a fallback arm can run on another day)
eid(){ uv run python scripts/slot_entry_id.py "$1" experiments/456-zs-riders/work/submit0928s[12].log experiments/479-pertarget-sharpen/work/submit1002s[12].log experiments/495-lodo-cap6/work/submit1003s[12].log experiments/498-combined-bulk/work/submit1004s[12].log 2>/dev/null || echo none; }
I453="$(eid 453dr)"; IP="$(eid 479p2g)"; IS="$(eid 479sp2)"; I495="$(eid 495)"
log "inputs: 453dr ${I453:-none}; 479p2g ${IP:-none}; 479sp2 ${IS:-none}; 495 ${I495:-none}"
RULE="$(uv run python $X/rule1005.py "${I453:-none}" "${IP:-none}" "${IS:-none}" "${I495:-none}" 2>&1 | tail -1)"
log "rule: $RULE"
pkg(){ case "$1" in
  496a) echo $S/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05.vcc;;
  498p2g) echo $S/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc;;
  498sp2) echo $S/498sp2-cap15-partial50-clean4q80-lodo6-ballast-sharpensp2-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc;;
  497dr) echo $S/497dr-october-path-lodopct-spt-abc.vcc;;
  460dr) echo $S/460dr-october-path-lodo-abc-q1.75.vcc;;
  461dr) echo $S/461dr-october-path-arcprog-abc-q1.75.vcc;;
  469dr) echo $S/469dr-october-path-arcprog-gate33-abc-q1.75.vcc;;
  458e) echo $S/458e-cap15-partial50-arch1eblend-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc;;
  *) echo "/nonexistent/$1";; esac; }
# -m names: 498p2g / 498sp2 / 497dr from the MODEL-NAMES.md queue (claimed 09-27, docs/model-names-claimed.tsv); the
# four October-path/programme fallbacks keep the names chain1001b baked for them on 09-26.
name(){ case "$1" in
  496a) echo "The world is bottlenecked by the ability of organizations to turn ideas into software. (John Collison)";;
  498p2g) echo "Even if I knew the world would end tomorrow, I'd still plant an apple tree.";;
  498sp2) echo "Perhaps, we don't need new ideas, but the same ideas to be known by everyone.";;
  497dr) echo "You're not respawning. What are you waiting for?";;
  461dr) echo "Microsoft has had clear competitors in the past. It's a good thing we have museums to document that. (Bill Gates)";; 460dr) echo "Hate the chains that hold you back. Hate the doubts that whisper lies. Hate the part of you that still seeks comfort in mediocrity. Use it. Let it drive you. Become the force that reshapes the world, or be forgotten by it.";;
  469dr) echo "When something truly matters and you possess the ability to create change, it becomes your moral duty to do so.";; 458e) echo "I would like to own all the land in this valley, as far as I can see.";; *) echo "";; esac; }
rl(){ case "$1" in
  497dr) echo experiments/497-october-switches/work/gates497.log;;
  460dr) echo experiments/460-lodo-dayof/work/gates460.log;;
  461dr) echo experiments/461-arcprog-october-path/work/gates461.log;;
  469dr) echo experiments/469-gated-october-path/work/gates469.log;;
  *) echo "$E/build_$1/run$1.log";; esac; }
nlj(){ case "$1" in
  498p2g) echo $E/build_498p2g/the-parafilm-joined-a-band.ntclevel.json;;
  498sp2) echo $E/build_498sp2/the-cuvette-learned-to-juggle.ntclevel.json;; esac; }
rho_ok(){ case "$1" in 498p2g|498sp2) ;; *) return 0;; esac
  uv run python - "$(nlj "$1")" <<'PY'
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
L1="${RULE%%|*}"; REST="${RULE#*|}"; L2="${REST%%|*}"
A1=""; A2=""
for c in $L1; do [ "$c" = "-" ] && continue; if ready "$c"; then A1=$c; break; else log "s1 candidate $c not ready -- skipped"; fi; done
for c in $L2; do [ "$c" = "-" ] && continue; [ "$c" = "$A1" ] && continue; if ready "$c"; then A2=$c; break; else log "s2 candidate $c not ready -- skipped"; fi; done
[ -z "$A1" ] && [ -n "$A2" ] && { A1=$A2; A2=""; }
[ -n "$A1" ] || { log "no eligible 10-05 arm -- nothing armed (the next session decides)"; log "END"; exit 0; }
nohup bash $X/arm1005.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 10-05" > $W/arm1005_s1.nohup 2>&1 & echo $! > $W/arm1005_s1.pid
[ -n "$A2" ] && { nohup bash $X/arm1005.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 10-05" > $W/arm1005_s2.nohup 2>&1 & echo $! > $W/arm1005_s2.pid; }
sleep 5
log "armed 10-05: s1 $A1 (pid $(cat $W/arm1005_s1.pid)), s2 ${A2:-<none>}"
log "END"
