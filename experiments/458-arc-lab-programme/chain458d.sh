#!/usr/bin/env bash
# 458d: after 460's chain ends (the builder is restored), build 446b + VCC_SAMELINE (Arc's H1 programme from ALL public cells,
# w 0) through run_arm_generic, exactly as chain458.sh built 458. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/458-arc-lab-programme; LOG=$X/work/chain458d.log; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain458d.pid
log "armed (pid $$): waiting for 460's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|def_day did not finish.*|relabel failed|prep failed|prep dry run failed)$" experiments/460-lodo-dayof/work/chain460.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain460\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262300 ] && { log "460 never ended by 23:00Z -- stopping"; exit 1; }; sleep 60; done
B=4feae27d37ffa43687a2eaa96d4ca835d85900fbaac63504dc1f5987fe6cd2e9; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$BF" = "$B" ] || { log "builder $BF is not $B (still pinned?) -- stopping"; exit 1; }
A=$X/work/arch1d_programme_abc.npz; AS=$(shasum -a 256 $A | cut -d' ' -f1)
[ "$AS" = "6951dbc7498e1f2aeabf405d5cbe94e7b695a7d3c0848a3a9460daa52d87a5ba" ] || { log "artifact sha $AS unexpected -- stopping"; exit 1; }
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
log "458d start (460 ended: $(tail -1 experiments/460-lodo-dayof/work/chain460.log | cut -c1-80))"
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$REPO/$A:$AS:0"
  bash $E/run_arm_generic.sh 458d 1.75 "$REF" the-pipette-learned-to-meditate \
    458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" sameline ) >> $X/work/build458d.nohup 2>&1
log "458d: $(/usr/bin/grep -E '^\[458d\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping' $E/work/build_458d/run458d.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_458d/the-pipette-learned-to-meditate.h5ad
log "END"
