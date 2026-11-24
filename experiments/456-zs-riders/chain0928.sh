#!/usr/bin/env bash
# Applies "the 09-28 rule, FIXED before any 09-27 entry exists" (docs/NEXT-SLOT-DECISION.md, 2026-09-26 01:40Z) once
# the 09-27 slots finish, then arms arm0928.sh s1/s2. Never builds, never submits by itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/456-zs-riders; W=$X/work; LOG=$W/chain0928.log; P=experiments/446-q175-level/work
E=experiments/428-zero-gene-audit/work; S=submissions
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain0928.pid
log "armed (pid $$): waiting for the 09-27 second slot to finish"
until /usr/bin/grep -qE "CAPTURED own entry|NOT SUBMITTING|exiting| END$" $P/submit0927s2.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271400 ] && { log "09-27 slot 2 never finished by 14:00Z -- stopping (PI's call)"; exit 1; }
  sleep 60; done
sleep 30
ID1="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P/submit0927s1.log 2>/dev/null | tail -1 | awk '{print $4}')"
ID2="$(/usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $P/submit0927s2.log 2>/dev/null | tail -1 | awk '{print $4}')"
log "inputs: 446bj2 ${ID1:-<none>}; 446br2 ${ID2:-<none>}"
BR="$(uv run python - "${ID1:-none}" "${ID2:-none}" <<'PY'
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
b = load("UDgjOKL9gI1udAexq6ue"); k = load(sys.argv[1]); f = load(sys.argv[2])
dk = (k["score_avg"] - b["score_avg"]) if k else None
key_ok = bool(k) and dk >= 0.001 and (k["pds_cosine"] - b["pds_cosine"]) >= 0
df = (f["score_avg"] - b["score_avg"]) if f else None
floor_ok = bool(f) and (f["de_wilcoxon_lfc_nmae"] - b["de_wilcoxon_lfc_nmae"]) <= -0.002
if key_ok and floor_ok: arm = "446bj2z" if dk >= df else "446br2z"
elif key_ok: arm = "446bj2z"
elif floor_ok: arm = "446br2z"
else: arm = "446bz"
print(arm, f"key_ok {key_ok} (d {dk}) floor_ok {floor_ok} (d {df})")
PY
)"
log "rule: $BR"
ZS="${BR%% *}"
# bash 3.2 (macOS /bin/bash) has no associative arrays: the -A flag fails and a keyed assignment silently indexes 0
pkg(){ case "$1" in
  446bz) echo $S/446bz-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5-zs250.vcc;;
  446bj2z) echo $S/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.vcc;;
  446br2z) echo $S/446br2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2-zs250.vcc;;
  446bd2) echo $S/446bd2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackdld1u2load-ctxlevelm-mde010coarse-zr5.vcc;;
  453dr) echo $S/453dr-october-path-abc-q1.75.vcc;;
  *) echo "/nonexistent/$1";; esac; }
name(){ case "$1" in
  446bz) echo "the pipette tip went on strike";; 446bj2z) echo "the fume hood took up yoga";; 446br2z) echo "the gel box wrote a haiku";;
  446bd2) echo "the microscope joined a book club";; 453dr) echo "the lab notebook ran a marathon";; *) echo "";; esac; }
rl(){ [ "$1" = "453dr" ] && echo experiments/453-abc-dress-rehearsal/work/gates453.log || echo "$E/build_$1/run$1.log"; }
ready(){ local r; r="$(rl "$1")"; [ -s "$(pkg "$1")" ] && /usr/bin/grep -q "^\[$1\] PASS" "$r" 2>/dev/null && /usr/bin/grep -q "^\[$1\] CLIFF OK" "$r" \
         && /usr/bin/grep -q "novelty rc=0" "$r" && /usr/bin/grep -q "package sha256:" "$r"; }
ready "$ZS" || { log "$ZS not ready (build failed or missing) -> 446bz"; ZS=446bz; }
if ready 453dr; then A1=453dr; A2=$ZS
else log "453 not ready -> slot 1 takes the zs arm"; A1=$ZS; A2=""; [ "$ZS" = "446bz" ] && ready 446bd2 && A2=446bd2; fi
nohup bash $X/arm0928.sh 1 "$A1" "$(rl $A1)" "$(pkg $A1)" "$(name $A1)" "$A1 slot 1 on 09-28" > $W/arm0928_s1.nohup 2>&1 & echo $! > $W/arm0928_s1.pid
if [ -n "$A2" ]; then nohup bash $X/arm0928.sh 2 "$A2" "$(rl $A2)" "$(pkg $A2)" "$(name $A2)" "$A2 slot 2 on 09-28" > $W/arm0928_s2.nohup 2>&1 & echo $! > $W/arm0928_s2.pid; fi
sleep 5
log "armed 09-28: s1 $A1 (pid $(cat $W/arm0928_s1.pid)), s2 ${A2:-<none>} ($( [ -n "$A2" ] && cat $W/arm0928_s2.pid ))"
log "END"
