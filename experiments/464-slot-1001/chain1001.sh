#!/usr/bin/env bash
# Applies "the 10-01 rule, FIXED 2026-09-26 before any 09-27 entry exists" (docs/NEXT-SLOT-DECISION.md) once the 09-30 slots
# finish, then arms arm1001.sh s1/s2. Never builds, never submits by itself.
#   key_ok : 09-27 s1 (446bj2) vs 446b: overall >= +0.001 and pds >= 0 (433's PREREGISTRATION-0926 bins)
#   P      : 09-30 s1 (458d, or the fallback 458) vs 446b by fid raw: positive >= +0.0015 (positive_reach if also reach
#            <= -0.005 raw and overall <= 0), partial +0.0005..+0.0015, else null/harm/missing (458's bins + addendum 11:31Z)
#   candidates, first two that are ready:
#     positive       -> 461dr, 458e, 460dr (if key_ok)
#     positive_reach -> 458e, 460dr (if key_ok)
#     partial        -> 458e, 460dr (if key_ok)
#     null/harm/none -> 460dr (if key_ok)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/464-slot-1001; W=$X/work; LOG=$W/chain1001.log; P27=experiments/446-q175-level/work; P30=experiments/459-h1-axis-def/work
E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain1001.pid
log "armed (pid $$): waiting for the 09-30 second slot to finish"
until /usr/bin/grep -qE "CAPTURED own entry|NOT SUBMITTING|exiting| END$" $P30/submit0930s2.log 2>/dev/null \
      || { /usr/bin/grep -qE "NOT SUBMITTING" $P30/submit0930s1.log 2>/dev/null; } \
      || { [ -f $P30/chain0930b.log ] && /usr/bin/grep -qE "no 09-30 slot armed" $P30/chain0930b.log; }; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609301800 ] && { log "09-30 slot 2 never finished by 18:00Z -- proceeding with what exists"; break; }
  sleep 60; done
sleep 30
ID27="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P27/submit0927s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
ID30="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P30/submit0930s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
TAG30="$(/usr/bin/grep -o 'ARMED slot 1: [A-Za-z0-9]*' $P30/submit0930s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
log "inputs: 09-27 s1 (446bj2) ${ID27:-<none>}; 09-30 s1 (${TAG30:-?}) ${ID30:-<none>}"
RULE="$(uv run python - "${ID27:-none}" "${ID30:-none}" "${TAG30:-none}" <<'PY'
import json, sys
from pathlib import Path
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
def load(i):
    p = Path(f"data/board-entries/{i}.json")
    return find(json.loads(p.read_text())) if i != "none" and p.exists() else None
b = load("UDgjOKL9gI1udAexq6ue"); k = load(sys.argv[1]); a = load(sys.argv[2]); tag = sys.argv[3]
key_ok = bool(k) and (k["score_avg"] - b["score_avg"]) >= 0.001 and (k["pds_cosine"] - b["pds_cosine"]) >= 0
P = "none"
if a and tag in ("458d", "458"):
    dfid = a["de_wilcoxon_direction_fidelity_yield_raw"] - b["de_wilcoxon_direction_fidelity_yield_raw"]
    dreach = a["de_wilcoxon_direction_reach_raw"] - b["de_wilcoxon_direction_reach_raw"]
    dov = a["score_avg"] - b["score_avg"]
    if dfid >= 0.0015:
        P = "positive_reach" if (dreach <= -0.005 and dov <= 0) else "positive"
    elif dfid >= 0.0005:
        P = "partial"
    else:
        P = "null"
    detail = f"dfid {dfid:+.5f} dreach {dreach:+.5f} dov {dov:+.7f}"
else:
    detail = "no programme entry"
cands = {"positive": ["461dr", "458e"], "positive_reach": ["458e"], "partial": ["458e"]}.get(P, [])
if key_ok: cands.append("460dr")
print(" ".join(cands) if cands else "-", f"| key_ok {key_ok} P {P} ({tag}: {detail})")
PY
)"
log "rule: $RULE"
pkg(){ case "$1" in
  461dr) echo $S/461dr-october-path-arcprog-abc-q1.75.vcc;;
  460dr) echo $S/460dr-october-path-lodo-abc-q1.75.vcc;;
  458e) echo $S/458e-cap15-partial50-arch1eblend-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc;;
  *) echo "/nonexistent/$1";; esac; }
name(){ case "$1" in
  461dr) echo "the bunsen burner took a gap year";; 460dr) echo "the fridge magnet learned calculus";;
  458e) echo "the lab coat filed a complaint";; *) echo "";; esac; }
rl(){ case "$1" in
  461dr) echo experiments/461-arcprog-october-path/work/gates461.log;; 460dr) echo experiments/460-lodo-dayof/work/gates460.log;;
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
