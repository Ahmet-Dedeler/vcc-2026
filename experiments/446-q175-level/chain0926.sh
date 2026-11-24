#!/usr/bin/env bash
# After the 09-26 reads: apply the 09-27 rule fixed in docs/NEXT-SLOT-DECISION.md (2026-09-25 18:40Z, before any 09-26
# entry existed), build what the branch needs (serial, builder 1dcb29a3 = b6c3f131 when VCC_SAMELINE is unset, KA-1), and
# arm both 09-27 slots with arm0927.sh. Reads only captured board records. Everything it decides is logged.
#   branch 1a  446b - 433k >= +0.001 with pds >= -0.0005, and 446bja - 446b >= +0.0010  -> build 446bjaj2, 446bjar2
#   branch 1b  446b - 433k >= +0.001 with pds >= -0.0005, otherwise                     -> 446bj2, 446br2 (built)
#   branch 2   446b - 433k <  +0.001 (or pds < -0.0005)                                  -> 433kj (built), build 433kr2
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/446-q175-level/work; LOG=$W/chain0926.log; X=experiments/450-same-line-programme
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "armed (pid $$): waiting for the 09-26 second slot to finish"
until grep -qE "CAPTURED own entry|NO SLOT 2|NOT SUBMITTING|exiting" $W/submit0926s2.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609261400 ] && { log "09-26 slot 2 never finished by 14:00Z -- stopping (PI's call)"; exit 1; }
  sleep 60; done
ID1="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $W/submit446b.log | tail -1 | awk '{print $4}')"
ID2="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' $W/submit0926s2.log | tail -1 | awk '{print $4}')"
TAG2="$(grep -o 'rule: [a-z0-9]*' $W/submit0926s2.log | tail -1 | awk '{print $2}')"
[ -n "$ID1" ] || { log "446b was not captured -- stopping (PI's call)"; exit 1; }
log "inputs: 446b $ID1; slot 2 ${TAG2:-<none>} ${ID2:-<none>}"
BR="$(uv run python - "$ID1" "${ID2:-none}" "${TAG2:-none}" <<'PY'
import json, sys
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
load = lambda i: find(json.load(open(f"data/board-entries/{i}.json")))
b, k = load(sys.argv[1]), load("T76JWD3tOajk5IYvzwUW")
d1, dp = b["score_avg"] - k["score_avg"], b["pds_cosine"] - k["pds_cosine"]
d2 = None
if sys.argv[2] != "none" and sys.argv[3] == "446bja":
    d2 = load(sys.argv[2])["score_avg"] - b["score_avg"]
q175 = d1 >= 0.001 and dp >= -0.0005
br = ("1a" if (d2 is not None and d2 >= 0.001) else "1b") if q175 else "2"
print(br, f"446b-433k {d1:+.7f} pds {dp:+.6f}; 446bja-446b " + ("n/a" if d2 is None else f"{d2:+.7f}"))
PY
)"
log "branch: $BR"
B=4feae27d37ffa436; BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
build(){  # build <arm> <artifact name> <package stem> <ref meta> <expect keys> ; env exported by the caller
  local arm=$1 nm=$2 stem=$3 ref=$4 exp=$5
  [ "${BF:0:16}" = "$B" ] || { log "builder is ${BF:0:16}, not $B -- stopping"; exit 1; }
  log "build $arm start"
  bash $X/run_arm_450.sh "$arm" "$ARM_Q" "$ref" "$nm" "$stem" "$BF" "$EM" "$exp" >> "$W/build_$arm.nohup" 2>&1
  local rl=$X/work/build_$arm/run$arm.log
  log "build $arm: $(grep -E "^\[$arm\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping" "$rl" | tr '\n' ' ' | cut -c1-500)"
  rm -f "$X/work/build_$arm/$nm.h5ad"
  grep -q "^\[$arm\] PASS" "$rl" && grep -q "^\[$arm\] CLIFF OK" "$rl" && grep -q "novelty rc=0" "$rl" || { log "build $arm failed its gates -- stopping"; exit 1; }; }
S=submissions; L446=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz
LODO=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz
case "${BR%% *}" in
  1a)
    ( export ARM_DE_SIGN_PRIOR=clean4q80ja VCC_BULK_LOAD="$LODO:$(shasum -a 256 $LODO | cut -d' ' -f1)" VCC_DE_MDE=0.1,0.1,0.3,coarse \
             VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" ARM_Q=1.75
      build 446bjaj2 the-lab-coat-applied-for-tenure 446bjaj2-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5 \
        "$REPO/$S/446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5.h5ad.meta.json" \
        de_sign_prior,de_sign_prior_context_A,de_sign_prior_sha256 ) || exit 1
    ( export ARM_DE_SIGN_PRIOR=clean4q80ja VCC_BULK_LOAD="$L446:$(shasum -a 256 $L446 | cut -d' ' -f1)" VCC_DE_MDE=0.1,0.1,0.3,coarse \
             VCC_DE_ZERO_REPAIR=-2 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" ARM_Q=1.75
      build 446bjar2 the-centrifuge-joined-a-choir 446bjar2-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2 \
        "$REPO/$S/446br2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2.h5ad.meta.json" \
        de_sign_prior,de_sign_prior_context_A,de_sign_prior_sha256 ) || exit 1
    A1=(446bjaj2 $X/work/build_446bjaj2/run446bjaj2.log $S/446bjaj2-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5.vcc
        "the lab coat applied for tenure" "446bjaj2 slot 1 on 09-27 (446bja + LODO key)")
    A2=(446bjar2 $X/work/build_446bjar2/run446bjar2.log $S/446bjar2-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2.vcc
        "the centrifuge joined a choir" "446bjar2 slot 2 on 09-27 (446bja + repair floor -2)");;
  1b)
    A1=(446bj2 experiments/428-zero-gene-audit/work/build_446bj2/run446bj2.log $S/446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5.vcc
        "the petri dish hired an agent" "446bj2 slot 1 on 09-27 (446b + LODO key at matched rho)")
    A2=(446br2 experiments/428-zero-gene-audit/work/build_446br2/run446br2.log $S/446br2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2.vcc
        "the pipette rack learned origami" "446br2 slot 2 on 09-27 (446b + repair floor -2)");;
  2)
    ( export VCC_BULK_LOAD="$L446:$(shasum -a 256 $L446 | cut -d' ' -f1)" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-2 \
             VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114" ARM_Q=1.5
      build 433kr2 the-ice-bucket-started-journaling 433kr2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr2 \
        "$REPO/$S/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json" de_zero_repair ) || exit 1
    A1=(433kj experiments/428-zero-gene-audit/work/build_433kj/run433kj.log $S/433kj-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stacklodo05load-ctxlevel-mde010coarse-zr5.vcc
        "the petri dish hired an agent" "433kj slot 1 on 09-27 (433k + LODO key)")
    A2=(433kr2 $X/work/build_433kr2/run433kr2.log $S/433kr2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr2.vcc
        "the ice bucket started journaling" "433kr2 slot 2 on 09-27 (433k + repair floor -2)");;
  *) log "branch did not resolve: $BR -- stopping"; exit 1;;
esac
nohup bash experiments/446-q175-level/arm0927.sh 1 "${A1[@]}" > $W/arm0927_s1.nohup 2>&1 & echo $! > $W/arm0927_s1.pid
nohup bash experiments/446-q175-level/arm0927.sh 2 "${A2[@]}" > $W/arm0927_s2.nohup 2>&1 & echo $! > $W/arm0927_s2.pid
sleep 5
log "armed 09-27: s1 ${A1[0]} (pid $(cat $W/arm0927_s1.pid)), s2 ${A2[0]} (pid $(cat $W/arm0927_s2.pid))"
log "END"
