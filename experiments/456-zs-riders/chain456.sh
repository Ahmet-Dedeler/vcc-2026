#!/usr/bin/env bash
# 456: pre-build the 09-28 slot-2 rider (the zs pass, VCC_DE_MAG=dn_k250_zs) on each base the 09-27 reads can produce, so
# 09-28 can be armed with no build on the day. 446bz (on 446b) already exists (8020a8d5...). This builds:
#   446bj2z = 446bj2 + zs   (config diff vs 446bj2 exactly de_mag)
#   446br2z = 446br2 + zs   (config diff vs 446br2 exactly de_mag)
# Waits for 455's chain to END (it runs the driver). run_arm_generic accepts the driver before or after 455.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/456-zs-riders; LOG="$X/work/chain456.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain456.pid
log "armed (pid $$): waiting for 455's chain to END"
# chain455.log quotes 453's "... END" line early on, so wait for 455's OWN last line ("driver now ...") and its exit
until /usr/bin/grep -q "driver now " experiments/455-driver-q175/work/chain455.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain455\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262000 ] && { log "455 not done by 20:00Z -- stopping"; exit 1; }; sleep 60; done
B=4feae27d37ffa43687a2eaa96d4ca835d85900fbaac63504dc1f5987fe6cd2e9; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; [ "$BF" = "$B" ] || { log "builder ${BF:0:16} != 4feae27d -- stopping"; exit 1; }
S=submissions
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "$DLS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DLS unexpected -- stopping"; exit 1; }
build(){ local arm=$1 ref=$2 nm=$3 stem=$4
  local scale; scale="$(uv run python -c "import json,sys; s=json.load(open(sys.argv[1]))['config']['bulk_ctx_scale']; print(','.join(f'{k}={v}' for k,v in s.items()))" "$ref")"
  log "$arm start (scale $scale)"
  bash $E/run_arm_generic.sh "$arm" 1.75 "$ref" "$nm" "$stem" "$B" "$EM" de_mag >> "$X/work/build_$arm.nohup" 2>&1
  log "$arm: $(/usr/bin/grep -E "^\[$arm\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping" $E/work/build_$arm/run$arm.log | tr '\n' ' ' | cut -c1-700)"
  rm -f "$E/work/build_$arm/$nm.h5ad"; }
( export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_DE_MAG=dn_k250_zs
  export VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681"
  build 446bj2z "$REPO/$S/446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5.h5ad.meta.json" \
    the-fume-hood-took-up-yoga 446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250 )
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-2 VCC_DE_MAG=dn_k250_zs
  export VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  build 446br2z "$REPO/$S/446br2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2.h5ad.meta.json" \
    the-gel-box-wrote-a-haiku 446br2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2-zs250 )
log "END"
