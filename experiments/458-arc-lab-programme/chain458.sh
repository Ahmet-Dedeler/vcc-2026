#!/usr/bin/env bash
# 458: after 457's chain ends, build 446b + VCC_SAMELINE (Arc's H1 programme for A/B/C, w 0) through run_arm_generic.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/458-arc-lab-programme; LOG=$X/work/chain458.log; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain458.pid
log "armed (pid $$): waiting for 457's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END$|stopping|did not finish|failed" experiments/457-cd4-half-def/work/chain457.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain457\.sh|def_day_cd4x15\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262300 ] && { log "457 never ended by 23:00Z -- stopping"; exit 1; }; sleep 60; done
B=4feae27d37ffa43687a2eaa96d4ca835d85900fbaac63504dc1f5987fe6cd2e9; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
A=$X/work/arch1_programme_abc.npz; AS=$(shasum -a 256 $A | cut -d' ' -f1)
[ "$AS" = "32b828501cd883180e65345bfe0205eb67e650ad2372b540c70ced8b022f62e5" ] || { log "artifact sha $AS unexpected -- stopping"; exit 1; }
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
log "458 start"
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$REPO/$A:$AS:0"
  bash $E/run_arm_generic.sh 458 1.75 "$REF" the-autoclave-read-the-manual \
    458-cap15-partial50-arch1prog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" sameline ) >> $X/work/build458.nohup 2>&1
log "458: $(/usr/bin/grep -E '^\[458\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping' $E/work/build_458/run458.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_458/the-autoclave-read-the-manual.h5ad
log "END"
