#!/usr/bin/env bash
# 479b: the arms P2g (tercile:2.5,2.0,1.75) and S_P2 (stercile:2.5,2.0,1.75) on 446bj2. Each arm: round 1 with 446bj2's
# factors -> ntclevel -> round 2 at factors x the EV-scale -> ntclevel check (rho within 0.5 % of 0.9793).
# Re-armed 2026-09-27 06:40Z after 499 F1: waits for chain495's process to exit (no other build), the builder at HEAD
# (495's b4975024 or later) and the driver with the per-target row (44714b9a+). chain476 waits for this chain (its
# regex names chain479b.sh). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/479-pertarget-sharpen; LOG=$X/work/chain479b.log; E=experiments/428-zero-gene-audit; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain479b.pid
log "armed (pid $$): waiting for chain495's process to exit (no other build), builder == HEAD, driver with the per-target row"
P495="$(cat experiments/495-lodo-cap6/work/chain495.pid 2>/dev/null)"
until ! ps -p "${P495:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain495 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain479\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] \
      && /usr/bin/grep -q '"tercile", "stercile"' scripts/def_pipeline.py \
      && [ "$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)" = "$(git show HEAD:scripts/def_pipeline.py | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609290800 ] && { log "not released by 08:00Z 09-29 -- stopping"; log "END"; exit 1; }; sleep 60; done
B="$(shasum -a 256 $BU | cut -d' ' -f1)"; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R="$REPO/$(ls submissions/446bj2-*.h5ad.meta.json | head -1)"
F0="A=1.0962,B=1.1461,C=1.1681"
log "released: builder ${B:0:16}; dump ${DS:0:16}"
build(){ local tag=$1 qv=$2 fact=$3 name=$4 stem=$5 expect=$6
  ( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG
    export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$fact"
    bash $E/run_arm_generic.sh "$tag" "$qv" "$R" "$name" "$stem" "$B" "$EM" "$expect" ) >> $X/work/build_$tag.nohup 2>&1
  log "$tag: $(/usr/bin/grep -E "^\[$tag\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping" $E/work/build_$tag/run$tag.log | tr '\n' ' ' | cut -c1-500)"
  /usr/bin/grep -q "^\[$tag\] PASS" $E/work/build_$tag/run$tag.log; }
level(){ uv run python experiments/417-bulk-sharpen/ntclevel417.py "$1" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1; }
for arm in "p2g:tercile:2.5,2.0,1.75:the-lyophilizer-wrote-a-cookbook" "sp2:stercile:2.5,2.0,1.75:the-magnetic-stirrer-ran-for-office"; do
  a="${arm%%:*}"; rest="${arm#*:}"; mode="${rest%%:*}"; rest="${rest#*:}"; qs="${rest%%:*}"; nm="${rest#*:}"; QV="$mode:$qs"
  t1="479${a}r1"; t2="479${a}"
  build "$t1" "$QV" "$F0" "$nm-r1" "$t1-cap15-partial50-clean4q80-pertscale3-ballast-sharpen${a}-stacklodo05load-ctxlevel-mde010coarse-zr5" \
        "bulk_sharpen_pt,bulk_sharpen_q" || { log "$t1 build failed -- next arm"; continue; }
  H1="$E/work/build_$t1/$nm-r1.h5ad"; level "$H1" || { log "$t1 ntclevel failed"; continue; }
  EV="$(/usr/bin/grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' $LOG | tail -1 | awk -F'= ' '{print $2}')"
  F1="$(uv run python -c "import sys; f0={k:float(v) for k,v in (p.split('=') for p in sys.argv[1].split(','))}; ev={k:float(v) for k,v in (p.split('=') for p in sys.argv[2].split(','))}; print(','.join(f'{k}={f0[k]*ev[k]:.4f}' for k in 'ABC'))" "$F0" "$EV")"
  log "$t1 rho: $(uv run python scripts/def_day_helpers.py rho "${H1%.h5ad}.ntclevel.json" 2>&1 | tail -1); EV $EV -> round-2 factors $F1"
  rm -f "$H1" "$REPO/submissions/$t1-cap15-partial50-clean4q80-pertscale3-ballast-sharpen${a}-stacklodo05load-ctxlevel-mde010coarse-zr5.vcc"  # round 1 is never submitted; keep its meta
  build "$t2" "$QV" "$F1" "$nm" "$t2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen${a}-stacklodo05load-ctxlevelr-mde010coarse-zr5" \
        "bulk_ctx_scale,bulk_sharpen_pt,bulk_sharpen_q" || { log "$t2 build failed"; continue; }
  H2="$E/work/build_$t2/$nm.h5ad"; level "$H2"
  log "$t2 rho: $(uv run python scripts/def_day_helpers.py rho "${H2%.h5ad}.ntclevel.json" 2>&1 | tail -1) (must be within 0.5 % of 0.9793)"
  rm -f "$H2"
done
log "END"
