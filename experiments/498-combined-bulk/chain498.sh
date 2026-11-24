#!/usr/bin/env bash
# 498: the combined bulk arms for 10-04's candidates (540c: 479 x 495 stack cleanly): 446bj2 + the cap-6 dump (495) + each
# per-target sharpening (P2g tercile:2.5,2.0,1.75 and S_P2 stercile:2.5,2.0,1.75), each built at 495's final factors,
# ntclevel'd, rebuilt at factors x EV-scale (rho within 0.5 % of 0.9793). Config diff vs 495: bulk_sharpen_q,
# bulk_sharpen_pt (+ bulk_ctx_scale in round 2). Serial after chain497 (process exit, no other build). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/498-combined-bulk; LOG=$X/work/chain498.log; E=experiments/428-zero-gene-audit; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain498.pid
log "armed (pid $$): waiting for chain497's process to exit (and no other build), builder == HEAD, 495 package present"
P497="$(cat experiments/497-october-switches/work/chain497.pid 2>/dev/null)"
until ! ps -p "${P497:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain497 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain47[69]b?\.sh|chain49[57]\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610011200 ] && { log "497 never ended by 12:00Z 10-01 -- stopping"; log "END"; exit 1; }; sleep 60; done
R="$REPO/submissions/495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json"
N2="$REPO/experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz"
[ -s "$R" ] && [ -s "$N2" ] || { log "495's meta or cap-6 dump missing -- stopping"; log "END"; exit 1; }
DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1)
[ "$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"]["bulk_load_sha256"])' "$R")" = "$DS2" ] || { log "495's meta does not carry this dump -- stopping"; log "END"; exit 1; }
F0="$(uv run python -c 'import json,sys; s=json.load(open(sys.argv[1]))["config"]["bulk_ctx_scale"]; print(",".join(f"{k}={s[k]:.4f}" for k in "ABC"))' "$R")"
B="$(shasum -a 256 $BU | cut -d' ' -f1)"; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
log "released: builder ${B:0:16}; cap-6 dump ${DS2:0:16}; 495's factors $F0"
build(){ local tag=$1 qv=$2 fact=$3 name=$4 stem=$5 expect=$6
  ( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
    export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$fact"
    bash $E/run_arm_generic.sh "$tag" "$qv" "$R" "$name" "$stem" "$B" "$EM" "$expect" ) >> $X/work/build_$tag.nohup 2>&1
  log "$tag: $(/usr/bin/grep -E "^\[$tag\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping" $E/work/build_$tag/run$tag.log | tr '\n' ' ' | cut -c1-500)"
  /usr/bin/grep -q "^\[$tag\] PASS" $E/work/build_$tag/run$tag.log; }
level(){ uv run python experiments/417-bulk-sharpen/ntclevel417.py "$1" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1; }
for arm in "p2g:tercile:2.5,2.0,1.75:the-parafilm-joined-a-band" "sp2:stercile:2.5,2.0,1.75:the-cuvette-learned-to-juggle"; do
  a="${arm%%:*}"; rest="${arm#*:}"; mode="${rest%%:*}"; rest="${rest#*:}"; qs="${rest%%:*}"; nm="${rest#*:}"; QV="$mode:$qs"
  t1="498${a}r1"; t2="498${a}"
  S1="$t1-cap15-partial50-clean4q80-lodo6-ballast-sharpen${a}-stacklodo05hi6load-ctxlevel-mde010coarse-zr5"
  build "$t1" "$QV" "$F0" "$nm-r1" "$S1" "bulk_sharpen_pt,bulk_sharpen_q" || { log "$t1 build failed -- next arm"; continue; }
  H1="$E/work/build_$t1/$nm-r1.h5ad"; level "$H1" || { log "$t1 ntclevel failed"; continue; }
  EV="$(/usr/bin/grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' $LOG | tail -1 | awk -F'= ' '{print $2}')"
  F1="$(uv run python -c "import sys; f0={k:float(v) for k,v in (p.split('=') for p in sys.argv[1].split(','))}; ev={k:float(v) for k,v in (p.split('=') for p in sys.argv[2].split(','))}; print(','.join(f'{k}={f0[k]*ev[k]:.4f}' for k in 'ABC'))" "$F0" "$EV")"
  log "$t1 rho: $(uv run python scripts/def_day_helpers.py rho "${H1%.h5ad}.ntclevel.json" 2>&1 | tail -1); EV $EV -> round-2 factors $F1"
  rm -f "$H1" "$REPO/submissions/$S1.vcc"
  build "$t2" "$QV" "$F1" "$nm" "$t2-cap15-partial50-clean4q80-lodo6-ballast-sharpen${a}-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5" \
        "bulk_ctx_scale,bulk_sharpen_pt,bulk_sharpen_q" || { log "$t2 build failed"; continue; }
  H2="$E/work/build_$t2/$nm.h5ad"; level "$H2"
  log "$t2 rho: $(uv run python scripts/def_day_helpers.py rho "${H2%.h5ad}.ntclevel.json" 2>&1 | tail -1) (must be within 0.5 % of 0.9793)"
  rm -f "$H2"
done
log "END"
